module med_diag_mod

  !----------------------------------------------------------------------------
  ! Compute spatial and time averages of fluxed quatities for water and
  ! energy balance
  !
  ! Sign convention for fluxes is positive downward with hierarchy being
  !    atm/glc/lnd/rof/ice/ocn
  ! Sign convention:
  !    positive value <=> the model is gaining water, heat, momentum, etc.
  ! Unit convention:
  !    heat flux     ~ W/m^2
  !    momentum flux ~ N/m^2
  !    water flux    ~ (kg/s)/m^2
  !    salt  flux    ~ (kg/s)/m^2
  !----------------------------------------------------------------------------

  use ESMF
  use shr_sys_mod           , only : shr_sys_abort
  use shr_const_mod         , only : shr_const_rearth, shr_const_pi, shr_const_latice
  use shr_const_mod         , only : shr_const_ice_ref_sal, shr_const_ocn_ref_sal, shr_const_isspval
  use med_constants_mod     , only : R8, CS, CL
  use med_internalstate_mod , only : InternalState, logunit, mastertask
  use shr_nuopc_methods_mod , only : fldchk => shr_nuopc_methods_FB_FldChk
  use shr_nuopc_utils_mod   , only : chkerr => shr_nuopc_utils_ChkErr
  use shr_nuopc_methods_mod , only : FB_GetFldPtr => shr_nuopc_methods_FB_GetFldPtr
  use esmFlds               , only : compatm, compocn, compice, complnd, comprof, compwav, compglc

  implicit none
  private

  public  :: med_diag_zero
  public  :: med_phases_diag_atm
  public  :: med_phases_diag_lnd
  public  :: med_phases_diag_rof
  public  :: med_phases_diag_glc
  public  :: med_phases_diag_ocn
  public  :: med_phases_diag_ice_ice2med
  public  :: med_phases_diag_ice_med2ice
  public  :: med_phases_diag_print

  private :: med_diag_accum
  private :: med_diag_sum_master
  private :: med_diag_print_atm
  private :: med_diag_print_lnd_ice_ocn
  private :: med_diag_print_summary

  integer, parameter :: budget_nmax = 50
  type, public :: budget_diag_type
     character(CS) :: name
     integer       :: index
  end type budg_diag_type
  type(budget_diag_type) :: budget_fields(budget_namx)
  type(budget_diag_type) :: budget_comps(budget_nmax)
  type(budget_diag_type) :: budget_periods(budget_nmax)

  ! ---------------------------------
  ! print options (obtained from mediator config input)
  ! ---------------------------------

  ! sets the diagnotics level of the annual budgets. [0,1,2,3],
  !  0 = none,
  !  1 = net summary budgets
  !  2 = 1 + detailed lnd/ocn/ice component budgets
  !  3 = 2 + detailed atm budgets

  integer :: budget_print_inst  ! default is 0
  integer :: budget_print_daily ! default is 0
  integer :: budget_print_month ! default is 1
  integer :: budget_print_ann   ! default is 1
  integer :: budget_print_ltann ! default is 1
  integer :: budget_print_ltend ! default is 0

  ! formats for output tables
  character(*), parameter :: F00  = "('(med_phases_diag_print) ',4a)"
  character(*), parameter :: FAH  = "(4a,i9,i6)"
  character(*), parameter :: FA0  = "('    ',12x,6(6x,a8,1x))"
  character(*), parameter :: FA1  = "('    ',a12,6f15.8)"
  character(*), parameter :: FA0r = "('    ',12x,8(6x,a8,1x))"
  character(*), parameter :: FA1r = "('    ',a12,8f15.8)"

  ! ---------------------------------
  ! local constants
  ! ---------------------------------
  real(r8), parameter :: HFLXtoWFLX = & ! water flux implied by latent heat of fusion
       &  - (shr_const_ocn_ref_sal-shr_const_ice_ref_sal) / &
       &    (shr_const_ocn_ref_sal*shr_const_latice)
  real(r8), parameter :: SFLXtoWFLX = & ! water flux implied by salt flux (kg/m^2s)
       -1._r8/(shr_const_ocn_ref_sal*1.e-3_r8)

  ! ---------------------------------
  ! C for component
  ! ---------------------------------

  ! "r" is receive in the mediator
  ! "s" is send from the mediator

  integer  :: c_atm_asend  ! model index: atm
  integer  :: c_atm_arecv  ! model index: atm
  integer  :: c_inh_isend  ! model index: ice, northern
  integer  :: c_inh_irecv  ! model index: ice, northern
  integer  :: c_ish_isend  ! model index: ice, southern
  integer  :: c_ish_irecv  ! model index: ice, southern
  integer  :: c_lnd_lsend  ! model index: lnd
  integer  :: c_lnd_lrecv  ! model index: lnd
  integer  :: c_ocn_osend  ! model index: ocn
  integer  :: c_ocn_orecv  ! model index: ocn
  integer  :: c_rof_rsend  ! model index: rof
  integer  :: c_rof_rrecv  ! model index: rof
  integer  :: c_glc_gsend  ! model index: glc
  integer  :: c_glc_grecv  ! model index: glc

  ! The folowing is needed for detailing the atm budgets and breakdown into components
  integer  :: c_inh_asend  ! model index: ice, northern, on atm grid
  integer  :: c_inh_arecv  ! model index: ice, northern, on atm grid
  integer  :: c_ish_asend  ! model index: ice, southern, on atm grid
  integer  :: c_ish_arecv  ! model index: ice, southern, on atm grid
  integer  :: c_lnd_asend  ! model index: lnd, on atm grid
  integer  :: c_lnd_arecv  ! model index: lnd, on atm grid
  integer  :: c_ocn_asend  ! model index: ocn, on atm grid
  integer  :: c_ocn_arecv  ! model index: ocn, on atm grid
  integer  :: c_size

  ! ---------------------------------
  ! F for field
  ! ---------------------------------

  integer :: f_area          ! area (wrt to unit sphere)
  integer :: f_heat_frz      ! heat : latent, freezing
  integer :: f_heat_melt     ! heat : latent, melting
  integer :: f_heat_swnet    ! heat : short wave, net
  integer :: f_heat_lwdn     ! heat : longwave down
  integer :: f_heat_lwup     ! heat : longwave up
  integer :: f_heat_latv     ! heat : latent, vaporization
  integer :: f_heat_latf     ! heat : latent, fusion, snow
  integer :: f_heat_ioff     ! heat : latent, fusion, frozen runoff
  integer :: f_heat_sen      ! heat : sensible
  integer :: f_watr_frz      ! water: freezing
  integer :: f_watr_melt     ! water: melting
  integer :: f_watr_rain     ! water: precip, liquid
  integer :: f_watr_snow     ! water: precip, frozen
  integer :: f_watr_evap     ! water: evaporation
  integer :: f_watr_salt     ! water: water equivalent of salt flux
  integer :: f_watr_roff     ! water: runoff/flood
  integer :: f_watr_ioff     ! water: frozen runoff
  integer :: f_watr_frz_16O  ! water isotope: freezing
  integer :: f_watr_melt_16O ! water isotope: melting
  integer :: f_watr_rain_16O ! water isotope: precip, liquid
  integer :: f_watr_snow_16O ! water isotope: prcip, frozen
  integer :: f_watr_evap_16O ! water isotope: evaporation
  integer :: f_watr_roff_16O ! water isotope: runoff/flood
  integer :: f_watr_ioff_16O ! water isotope: frozen runoff
  integer :: f_watr_frz_18O  ! water isotope: freezing
  integer :: f_watr_melt_18O ! water isotope: melting
  integer :: f_watr_rain_18O ! water isotope: precip, liquid
  integer :: f_watr_snow_18O ! water isotope: precip, frozen
  integer :: f_watr_evap_18O ! water isotope: evaporation
  integer :: f_watr_roff_18O ! water isotope: runoff/flood
  integer :: f_watr_ioff_18O ! water isotope: frozen runoff
  integer :: f_watr_frz_HDO  ! water isotope: freezing
  integer :: f_watr_melt_HDO ! water isotope: melting
  integer :: f_watr_rain_HDO ! water isotope: precip, liquid
  integer :: f_watr_snow_HDO ! water isotope: precip, frozen
  integer :: f_watr_evap_HDO ! water isotope: evaporation
  integer :: f_watr_roff_HDO ! water isotope: runoff/flood
  integer :: f_watr_ioff_HDO ! water isotope: frozen runoff
  integer :: f_size          ! Total array size of all elements

  integer :: f_heat_beg      ! 1st index  for heat
  integer :: f_heat_end      ! Last index for heat
  integer :: f_watr_beg      ! 1st index  for water
  integer :: f_watr_end      ! Last index for water

  integer :: f_16O_beg       ! 1st index  for 16O water isotope
  integer :: f_16O_end       ! Last index for 16O water isotope
  integer :: f_18O_beg       ! 1st index  for 18O water isotope
  integer :: f_18O_end       ! Last index for 18O water isotope
  integer :: f_HDO_beg       ! 1st index  for HDO water isotope
  integer :: f_HDO_end       ! Last index for HDO water isotope

  ! water isotopes names and indices
  logical :: flds_wiso  = .false.! If water isotope fields are active -
  ! TODO: for now set to .false. - but this needs to be set in an initialization phase

  integer, parameter          :: nisotopes = 3
  integer, parameter          :: iso0(nisotopes)    = (/ f_16O_beg, f_18O_beg, f_hdO_beg /)
  integer, parameter          :: isof(nisotopes)    = (/ f_16O_end, f_18O_end, f_hdO_end /)
  character(len=5), parameter :: isoname(nisotopes) = (/ 'H216O',   'H218O',   '  HDO'   /)

  ! ---------------------------------
  ! P for period
  ! ---------------------------------

  integer :: period_inst
  integer :: period_day
  integer :: period_mon
  integer :: period_ann
  integer :: period_inf
  integer :: p_size

  ! ---------------------------------
  ! public data members
  ! ---------------------------------

  ! note: call med_diag_sum_master then save budget_global and budget_counter on restart from/to root pe ---

  real(r8),public, allocatable :: budget_local  (:,:,:) ! local sum, valid on all pes
  real(r8),public, allocatable :: budget_global (:,:,:) ! global sum, valid only on root pe
  real(r8),public, allocatable :: budget_counter(:,:,:) ! counter, valid only on root pe

  character(len=*), parameter :: modName   = "(med_diag) "
  character(len=*), parameter :: u_FILE_u  = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_diag_init(gcomp, rc)

    c_size = 1
    call add_to_budget_type(budget_comps , c_size, c_atm_asend, ' c2a_atm' ) ! model index: atm
    call add_to_budget_type(budget_comps , c_size, c_atm_arecv, ' a2c_atm' ) ! model index: atm
    call add_to_budget_type(budget_comps , c_size, c_inh_isend, ' c2i_inh' ) ! model index: ice, northern
    call add_to_budget_type(budget_comps , c_size, c_inh_irecv, ' i2c_inh' ) ! model index: ice, northern
    call add_to_budget_type(budget_comps , c_size, c_ish_isend, ' c2i_ish' ) ! model index: ice, southern
    call add_to_budget_type(budget_comps , c_size, c_ish_irecv, ' i2c_ish' ) ! model index: ice, southern
    call add_to_budget_type(budget_comps , c_size, c_lnd_lsend, ' c2l_lnd' ) ! model index: lnd
    call add_to_budget_type(budget_comps , c_size, c_lnd_lrecv, ' l2c_lnd' ) ! model index: lnd
    call add_to_budget_type(budget_comps , c_size, c_ocn_osend, ' c2o_ocn' ) ! model index: ocn
    call add_to_budget_type(budget_comps , c_size, c_ocn_orecv, ' o2c_ocn' ) ! model index: ocn
    call add_to_budget_type(budget_comps , c_size, c_rof_rsend, ' c2r_rof' ) ! model index: rof
    call add_to_budget_type(budget_comps , c_size, c_rof_rrecv, ' r2c_rof' ) ! model index: rof
    call add_to_budget_type(budget_comps , c_size, c_glc_gsend, ' c2g_glc' ) ! model index: glc
    call add_to_budget_type(budget_comps , c_size, c_glc_grecv, ' g2c_glc' ) ! model index: glc
    call add_to_budget_type(budget_comps , c_size, c_inh_asend, ' c2a_inh' ) ! model index: ice, northern, on atm grid
    call add_to_budget_type(budget_comps , c_size, c_inh_arecv, ' a2c_inh' ) ! model index: ice, northern, on atm grid
    call add_to_budget_type(budget_comps , c_size, c_ish_asend, ' c2a_ish' ) ! model index: ice, southern, on atm grid
    call add_to_budget_type(budget_comps , c_size, c_ish_arecv, ' a2c_ish' ) ! model index: ice, southern, on atm grid
    call add_to_budget_type(budget_comps , c_size, c_lnd_asend, ' c2a_lnd' ) ! model index: lnd, on atm grid
    call add_to_budget_type(budget_comps , c_size, c_lnd_arecv, ' a2c_lnd' ) ! model index: lnd, on atm grid
    call add_to_budget_type(budget_comps , c_size, c_ocn_asend, ' c2a_ocn' ) ! model index: ocn, on atm grid
    call add_to_budget_type(budget_comps , c_size, c_ocn_arecv, ' a2c_ocn' ) ! model index: ocn, on atm grid

    f_size = 1
    call add_to_budget_type(budget_fields , f_size, f_area          ,'area'        ) ! area (wrt to unit sphere)
    call add_to_budget_type(budget_fields , f_size, f_heat_frz      ,'hfreeze'     ) ! heat : latent, freezing
    call add_to_budget_type(budget_fields , f_size, f_heat_melt     ,'hmelt'       ) ! heat : latent, melting
    call add_to_budget_type(budget_fields , f_size, f_heat_swnet    ,'hnetsw'      ) ! heat : short wave, net
    call add_to_budget_type(budget_fields , f_size, f_heat_lwdn     ,'hlwdn'       ) ! heat : longwave down
    call add_to_budget_type(budget_fields , f_size, f_heat_lwup     ,'hlwup'       ) ! heat : longwave up
    call add_to_budget_type(budget_fields , f_size, f_heat_latv     ,'hlatvap'     ) ! heat : latent, vaporization
    call add_to_budget_type(budget_fields , f_size, f_heat_latf     ,'hlatfus'     ) ! heat : latent, fusion, snow
    call add_to_budget_type(budget_fields , f_size, f_heat_ioff     ,'hiroff'      ) ! heat : latent, fusion, frozen runoff
    call add_to_budget_type(budget_fields , f_size, f_heat_sen      ,'hsen'        ) ! heat : sensible
    call add_to_budget_type(budget_fields , f_size, f_watr_frz      ,'wfreeze'     ) ! water: freezing
    call add_to_budget_type(budget_fields , f_size, f_watr_melt     ,'wmelt'       ) ! water: melting
    call add_to_budget_type(budget_fields , f_size, f_watr_rain     ,'wrain'       ) ! water: precip, liquid
    call add_to_budget_type(budget_fields , f_size, f_watr_snow     ,'wsnow'       ) ! water: precip, frozen
    call add_to_budget_type(budget_fields , f_size, f_watr_evap     ,'wevap'       ) ! water: evaporation
    call add_to_budget_type(budget_fields , f_size, f_watr_salt     ,'weqsaltf'    ) ! water: water equivalent of salt flux
    call add_to_budget_type(budget_fields , f_size, f_watr_roff     ,'wrunoff'     ) ! water: runoff/flood
    call add_to_budget_type(budget_fields , f_size, f_watr_ioff     ,'wfrzrof'     ) ! water: frozen runoff
    call add_to_budget_type(budget_fields , f_size, f_watr_frz_16O  ,'wfreeze_16O' ) ! water isotope: freezing
    call add_to_budget_type(budget_fields , f_size, f_watr_melt_16O ,'wmelt_16O'   ) ! water isotope: melting
    call add_to_budget_type(budget_fields , f_size, f_watr_rain_16O ,'wrain_16O'   ) ! water isotope: precip, liquid
    call add_to_budget_type(budget_fields , f_size, f_watr_snow_16O ,'wsnow_16O'   ) ! water isotope: prcip, frozen
    call add_to_budget_type(budget_fields , f_size, f_watr_evap_16O ,'wevap_16O'   ) ! water isotope: evaporation
    call add_to_budget_type(budget_fields , f_size, f_watr_roff_16O ,'wrunoff_16O' ) ! water isotope: runoff/flood
    call add_to_budget_type(budget_fields , f_size, f_watr_ioff_16O ,'wfrzrof_16O' ) ! water isotope: frozen runoff
    call add_to_budget_type(budget_fields , f_size, f_watr_frz_18O  ,'wfreeze_18O' ) ! water isotope: freezing
    call add_to_budget_type(budget_fields , f_size, f_watr_melt_18O ,'wmelt_18O'   ) ! water isotope: melting
    call add_to_budget_type(budget_fields , f_size, f_watr_rain_18O ,'wrain_18O'   ) ! water isotope: precip, liquid
    call add_to_budget_type(budget_fields , f_size, f_watr_snow_18O ,'wsnow_18O'   ) ! water isotope: precip, frozen
    call add_to_budget_type(budget_fields , f_size, f_watr_evap_18O ,'wevap_18O'   ) ! water isotope: evaporation
    call add_to_budget_type(budget_fields , f_size, f_watr_roff_18O ,'wrunoff_18O' ) ! water isotope: runoff/flood
    call add_to_budget_type(budget_fields , f_size, f_watr_ioff_18O ,'wfrzrof_18O' ) ! water isotope: frozen runoff
    call add_to_budget_type(budget_fields , f_size, f_watr_frz_HDO  ,'wfreeze_HDO' ) ! water isotope: freezing
    call add_to_budget_type(budget_fields , f_size, f_watr_melt_HDO ,'wmelt_HDO'   ) ! water isotope: melting
    call add_to_budget_type(budget_fields , f_size, f_watr_rain_HDO ,'wrain_HDO'   ) ! water isotope: precip, liquid
    call add_to_budget_type(budget_fields , f_size, f_watr_snow_HDO ,'wsnow_HDO'   ) ! water isotope: precip, frozen
    call add_to_budget_type(budget_fields , f_size, f_watr_evap_HDO ,'wevap_HDO'   ) ! water isotope: evaporation
    call add_to_budget_type(budget_fields , f_size, f_watr_roff_HDO ,'wrunoff_HDO' ) ! water isotope: runoff/flood
    call add_to_budget_type(budget_fields , f_size, f_watr_ioff_HDO ,'wfrzrof_HDO' ) ! water isotope: frozen runoff

    f_heat_beg = f_heat_frz      ! 1st index  for heat
    f_heat_end = f_heat_sen      ! Last index for heat
    f_watr_beg = f_watr_frz      ! 1st index  for water
    f_watr_end = f_watr_ioff     ! Last index for water

    f_16O_beg  = f_watr_frz_16O  ! 1st index  for 16O water isotope
    f_16O_end  = f_watr_ioff_16O ! Last index for 16O water isotope
    f_18O_beg  = f_watr_frz_18O  ! 1st index  for 18O water isotope
    f_18O_end  = f_watr_ioff_18O ! Last index for 18O water isotope
    f_HDO_beg  = f_watr_frz_HDO  ! 1st index  for HDO water isotope
    f_HDO_end  = f_watr_ioff_HDO ! Last index for HDO water isotope

    p_size = 1
    call add_to_budget_type(budget_fields , p_size, period_inst,'    inst')
    call add_to_budget_type(budget_fields , p_size, period_day ,'   daily')
    call add_to_budget_type(budget_fields , p_size, period_mon ,' monthly')
    call add_to_budget_type(budget_fields , p_size, period_ann ,'  annual')
    call add_to_budget_type(budget_fields , p_size, period_inf ,'all_time')

    ! allocate module budget arrays
    allocate(budget_local  (f_size, c_size, p_size)) ! local sum, valid on all pes
    allocate(budget_global (f_size, c_size, p_size)) ! global sum, valid only on root pe
    allocate(budget_counter(f_size, c_size, p_size)) ! counter, valid only on root pe

  end subroutine med_diag_init

  !===============================================================================

  subroutine med_diag_zero( gcomp, mode, rc)

    ! Zero out global budget diagnostic data.

    ! input/output variables
    type(ESMF_GridComp)                   :: gcomp
    character(len=*), intent(in),optional :: mode
    integer, intent(out)                  :: rc

    ! local variables
    type(ESMF_Clock) :: clock
    type(ESMF_Time)  :: currTime
    integer          :: ip
    integer          :: curr_year, curr_mon, curr_day, curr_tod
    character(*), parameter :: subName = '(med_diag_zero) '
    ! ------------------------------------------------------------------

    if(present(mode)) then

       if (trim(mode) == 'inst') then
          budget_local(:,:,period_inst) = 0.0_r8
          budget_global(:,:,period_inst) = 0.0_r8
          budget_counter(:,:,period_inst) = 0.0_r8
       elseif (trim(mode) == 'day') then
          budget_local(:,:,period_day) = 0.0_r8
          budget_global(:,:,period_day) = 0.0_r8
          budget_counter(:,:,period_day) = 0.0_r8
       elseif (trim(mode) == 'mon') then
          budget_local(:,:,period_mon) = 0.0_r8
          budget_global(:,:,period_mon) = 0.0_r8
          budget_counter(:,:,period_mon) = 0.0_r8
       elseif (trim(mode) == 'ann') then
          budget_local(:,:,period_ann) = 0.0_r8
          budget_global(:,:,period_ann) = 0.0_r8
          budget_counter(:,:,period_ann) = 0.0_r8
       elseif (trim(mode) == 'inf') then
          budget_local(:,:,period_inf) = 0.0_r8
          budget_global(:,:,period_inf) = 0.0_r8
          budget_counter(:,:,period_inf) = 0.0_r8
       elseif (trim(mode) == 'all') then
          budget_local(:,:,:) = 0.0_r8
          budget_global(:,:,:) = 0.0_r8
          budget_counter(:,:,:) = 0.0_r8
       else
          call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
       endif

    else
       call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet( currTime, yy=curr_year, mm=curr_mon, dd=curr_day, s=curr_tod, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       do ip = 1,period_ntypes
          if (ip == period_inst) then
             budget_local(:,:,ip) = 0.0_r8
             budget_global(:,:,ip) = 0.0_r8
             budget_counter(:,:,ip) = 0.0_r8
          endif
          if (ip==period_day .and. curr_tod==0) then
             budget_local(:,:,ip) = 0.0_r8
             budget_global(:,:,ip) = 0.0_r8
             budget_counter(:,:,ip) = 0.0_r8
          endif
          if (ip==period_mon .and. curr_day==1 .and. curr_tod==0) then
             budget_local(:,:,ip) = 0.0_r8
             budget_global(:,:,ip) = 0.0_r8
             budget_counter(:,:,ip) = 0.0_r8
          endif
          if (ip==period_ann .and. curr_mon==1 .and. curr_day==1 .and. curr_tod==0) then
             budget_local(:,:,ip) = 0.0_r8
             budget_global(:,:,ip) = 0.0_r8
             budget_counter(:,:,ip) = 0.0_r8
          endif
       enddo
    end if

  end subroutine med_diag_zero

  !===============================================================================

  subroutine med_diag_accum()

    ! ------------------------------------------------------------------
    ! Accumulate out global budget diagnostic data.
    ! ------------------------------------------------------------------

    ! local variables
    integer     :: ip
    character(*), parameter :: subName = '(med_diag_accum) '
    ! ------------------------------------------------------------------

    do ip = period_inst+1,period_ntypes
       budget_local(:,:,ip) = budget_local(:,:,ip) + budget_local(:,:,period_inst)
    enddo
    budget_counter(:,:,:) = budget_counter(:,:,:) + 1.0_r8

  end subroutine med_diag_accum

  !===============================================================================

  subroutine med_diag_sum_master(gcomp, rc)

    ! ------------------------------------------------------------------
    ! Sum local values to global on root
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM) :: vm
    integer       :: count
    real(r8), allocatable   :: budget_global_tmp(:,:,:) ! temporary sum
    character(*), parameter :: subName = '(med_diag_sum_master) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp, vm=vm, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(budget_global_tmp(f_size,c_size,period_ntypes)

    count = size(budget_global_tmp)
    budget_global_tmp(:,:,:) = 0.0_r8

    call ESMF_VMReduce(vm, budget_local, budget_global_tmp, count, ESMF_REDUCE_SUM, 0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    budget_global(:,:,:) = budget_global(:,:,:) + budget_global_tmp(:,:,:)
    budget_local(:,:,:) = 0.0_r8

    deallocate(budget_global_tmp)

  end subroutine med_diag_sum_master

  !===============================================================================

  subroutine med_phases_diag_atm(gcomp,  rc)

    ! ------------------------------------------------------------------
    ! Compute global atm input/output flux diagnostics
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(in)           :: gcomp
    integer             , intent(out)          :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n,k,ic,ip
    real(r8)            :: wgt, area
    logical             :: atm_wiso_recv = .false.
    logical             :: atm_wiso_send = .false.
    real(r8), pointer   :: afrac(:)
    real(r8), pointer   :: lfrac(:)
    real(r8), pointer   :: ifrac(:)
    real(r8), pointer   :: ofrac(:)
    real(r8), pointer   :: Faxa_swnet(:)        ! from atm to mediator
    real(r8), pointer   :: Faxa_lwdn(:)         ! from atm to mediator
    real(r8), pointer   :: Faxa_rainc(:)        ! from atm to mediator
    real(r8), pointer   :: Faxa_rainl(:)        ! from atm to mediator
    real(r8), pointer   :: Faxa_snowc(:)        ! from atm to mediator
    real(r8), pointer   :: Faxa_snowl(:)        ! from atm to mediator
    real(r8), pointer   :: Faxa_rainl_wiso(:,:) ! from atm to mediator
    real(r8), pointer   :: Faxa_rainc_wiso(:,:) ! from atm to mediator
    real(r8), pointer   :: Faxx_lwup(:)         ! from  mediator to atm
    real(r8), pointer   :: Faxx_lat(:)          ! from  mediator to atm
    real(r8), pointer   :: Faxx_sen(:)          ! from  mediator to atm
    real(r8), pointer   :: Faxx_evap(:)         ! from  mediator to atm
    real(r8), pointer   :: Faxx_evap_wiso(:,:)  ! from  mediator to atm
    character(*), parameter :: subName = '(med_phases_diag_atm) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get fractions on atm mesh
    call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'afrac', fldptr1=afrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', fldptr1=lfrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ifrac', fldptr1=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', fldptr1=ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------
    ! from atm to mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swnet', fldptr1=Faxa_swnet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , fldptr1=Faxa_lwdn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc', fldptr1=Faxa_rainc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl', fldptr1=Faxa_rainl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc', fldptr1=Faxa_snowc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl', fldptr1=Faxa_snowl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc_wiso', rc=rc)) then
       atm_wiso_recv = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl_wiso', fldptr2=Faxa_rainl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc_wiso', fldptr2=Faxa_rainc_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    do n = 1,size(Faxa_swnet)
       do k = 1,4  ! loop over components that atm sends data to

          area = is_local%wrap%mesh_info(compatm)%areas(n)
          if (k == 1) then
             wgt = -area * afrac(n)
          elseif (k == 2) then
             wgt =  area * lfrac(n)
          elseif (k == 3) then
             wgt =  area * ofrac(n)
          elseif (k == 4) then
             wgt =  area * ifrac(n)
          endif

          if (k == 1) then
             ic = c_atm_arecv
          elseif (k == 2) then
             ic = c_lnd_arecv
          elseif (k == 3) then
             ic = c_ocn_arecv
          elseif (k == 4) then
             if (is_local%wrap%mesh_info(compatm)%lats(n) > 0.0_r8) then
                ic = c_inh_arecv
             else
                ic = c_ish_arecv
             endif
          endif

          budget_local(f_area      ,ic,ip) = budget_local(f_area      ,ic,ip) + wgt
          budget_local(f_heat_swnet,ic,ip) = budget_local(f_heat_swnet,ic,ip) + wgt*Faxa_swnet(n)
          budget_local(f_heat_lwdn ,ic,ip) = budget_local(f_heat_lwdn ,ic,ip) + wgt*Faxa_lwdn(n)
          budget_local(f_watr_rain ,ic,ip) = budget_local(f_watr_rain ,ic,ip) + wgt*(Faxa_rainc(n) + Faxa_rainl(n))
          budget_local(f_watr_snow ,ic,ip) = budget_local(f_watr_snow ,ic,ip) + wgt*(Faxa_snowc(n) + Faxa_snowl(n))
          if ( atm_wiso_recv )then
             budget_local(f_watr_rain_16O,ic,ip) = budget_local(f_watr_rain_16O,ic,ip) + &
                  wgt*(Faxa_rainc_wiso(1,n) + Faxa_rainl_wiso(1,n))
             budget_local(f_watr_rain_18O,ic,ip) = budget_local(f_watr_rain_18O,ic,ip) + &
                  wgt*(Faxa_rainc_wiso(2,n) + Faxa_rainl_wiso(2,n))
             budget_local(f_watr_rain_HDO,ic,ip) = budget_local(f_watr_rain_HDO,ic,ip) + &
                  wgt*(Faxa_rainc_wiso(3,n) + Faxa_rainl_wiso(3,n))
          end if
       enddo
    enddo

    ! heat implied by snow flux
    budget_local(f_heat_latf,c_atm_arecv,ip) = -budget_local(f_watr_snow,c_atm_arecv,ip)*shr_const_latice
    budget_local(f_heat_latf,c_lnd_arecv,ip) = -budget_local(f_watr_snow,c_lnd_arecv,ip)*shr_const_latice
    budget_local(f_heat_latf,c_ocn_arecv,ip) = -budget_local(f_watr_snow,c_ocn_arecv,ip)*shr_const_latice
    budget_local(f_heat_latf,c_inh_arecv,ip) = -budget_local(f_watr_snow,c_inh_arecv,ip)*shr_const_latice
    budget_local(f_heat_latf,c_ish_arecv,ip) = -budget_local(f_watr_snow,c_ish_arecv,ip)*shr_const_latice

    !-------------------------------
    ! from mediator to atm
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_lwup', fldptr1=Faxx_lwup, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_lat' , fldptr1=Faxx_lat , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_sen' , fldptr1=Faxx_sen , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_evap', fldptr1=Faxx_evap, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBExp(compatm), 'Faxx_evap_wiso', rc=rc)) then
       atm_wiso_send = .true.
       call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_evap_wiso', fldptr2=Faxx_evap_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    do n = 1,size(Faxx_lwup)

       ! TODO: this loop must be over the number of components that are being merged to obtain the
       ! atm input plus 1
       do k = 1,4  ! sum over merged, lfrac, ifrac, ofrac (3 inputs to the atm, 1 merged input)

          ! Determine weight
          area = is_local%wrap%mesh_info(compatm)%areas(n)
          if (k == 1) then
             wgt = -area * afrac(n)
          else if (k == 2) then
             wgt =  area * lfrac(n)
          elseif (k == 3) then
             wgt =  area * ifrac(n)
          elseif (k == 4) then
             wgt =  area * ofrac(n)
          endif

          ! Determine budget_local index
          if (k == 1) then
             ic = c_atm_asend
          elseif (k == 2) then
             ic = c_lnd_asend
          elseif (k == 3) then
             ic = c_ocn_asend
          elseif (k == 4) then
             if (is_local%wrap%mesh_info(compatm)%lats(n) > 0.0_r8) then
                ic = c_inh_asend
             else
                ic = c_ish_asend
             endif
          endif

          budget_local(f_area     ,ic,ip) = budget_local(f_area     ,ic,ip) + wgt
          budget_local(f_heat_lwup,ic,ip) = budget_local(f_heat_lwup,ic,ip) + wgt*Faxx_lwup(n)
          budget_local(f_heat_latv,ic,ip) = budget_local(f_heat_latv,ic,ip) + wgt*Faxx_lat(n)
          budget_local(f_heat_sen ,ic,ip) = budget_local(f_heat_sen ,ic,ip) + wgt*Faxx_sen(n)
          budget_local(f_watr_evap,ic,ip) = budget_local(f_watr_evap,ic,ip) + wgt*Faxx_evap(n)
          if ( atm_wiso_send )then
             budget_local(f_watr_evap_16O,ic,ip) = budget_local(f_watr_evap_16O,ic,ip) + wgt*Faxx_evap_wiso(1,n)
             budget_local(f_watr_evap_18O,ic,ip) = budget_local(f_watr_evap_18O,ic,ip) + wgt*Faxx_evap_wiso(2,n)
             budget_local(f_watr_evap_HDO,ic,ip) = budget_local(f_watr_evap_HDO,ic,ip) + wgt*Faxx_evap_wiso(3,n)
          end if

       enddo
    enddo

  end subroutine med_phases_diag_atm

  !===============================================================================

  subroutine med_phases_diag_lnd( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global lnd input/output flux diagnostics
    ! ------------------------------------------------------------------

    ! intput/output variables
    type(ESMF_GridComp) , intent(in)  :: gcomp
    integer             , intent(out) :: rc

    ! local variables
    type(InternalState)    :: is_local
    real(r8), pointer      :: lfrac(:)
    integer                :: n,ic,ip
    real(r8)               :: wgt
    logical                :: lnd_wiso_recv = .false.
    logical                :: lnd_wiso_send = .false.
    logical                :: lnd_irrig = .false.
    real(r8), pointer      :: Fall_swnet(:)        ! to mediator from land
    real(r8), pointer      :: Fall_lwup(:)         ! to mediator from land
    real(r8), pointer      :: Fall_lat(:)          ! to mediator from land
    real(r8), pointer      :: Fall_sen(:)          ! to mediator from land
    real(r8), pointer      :: Fall_evap(:)         ! to mediator from land
    real(r8), pointer      :: Fall_evap_wiso(:,:)  ! to mediator from land
    real(r8), pointer      :: Flrl_rofsur(:)       ! to mediator from land
    real(r8), pointer      :: Flrl_rofgwl(:)       ! to mediator from land
    real(r8), pointer      :: Flrl_rofsub(:)       ! to mediator from land
    real(r8), pointer      :: Flrl_rofdto(:)       ! to mediator from land
    real(r8), pointer      :: Flrl_rofi(:)         ! to mediator from land
    real(r8), pointer      :: Flrl_irrig(:)        ! to mediator from land
    real(r8), pointer      :: Flrl_rofl_wiso(:,:)  ! to mediator from land
    real(r8), pointer      :: Flrl_rofi_wiso(:,:)  ! to mediator from land
    real(r8), pointer      :: Faxa_lwdn(:)         ! from mediator to land
    real(r8), pointer      :: Faxa_rainc(:)        ! from mediator to land
    real(r8), pointer      :: Faxa_rainl(:)        ! from mediator to land
    real(r8), pointer      :: Faxa_snowc(:)        ! from mediator to land
    real(r8), pointer      :: Faxa_snowl(:)        ! from mediator to land
    real(r8), pointer      :: Flrl_flood(:)        ! from mediator to land
    real(r8), pointer      :: Faxa_rainl_wiso(:,:) ! from mediator to land
    real(r8), pointer      :: Faxa_rainc_wiso(:,:) ! from mediator to land
    real(r8), pointer      :: Faxa_snowl_wiso(:,:) ! from mediator to land
    real(r8), pointer      :: Faxa_snowc_wiso(:,:) ! from mediator to land
    real(r8), pointer      :: Flrl_flood_wiso(:,:) ! from mediator to land
    character(*), parameter :: subName = '(med_phases_diag_lnd) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get fractions on lnd mesh
    call FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrac', lfrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------
    ! from land to mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_swnet', fldptr1=Fall_swnet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_lwup', fldptr1=Fall_lwup, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_lat', fldptr1=Fall_lat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_sen', fldptr1=Fall_sen, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap', fldptr1=Fall_evap, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofsur', fldptr1=Flrl_rofsur, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofgwl', fldptr1=Flrl_rofgwl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofsub', fldptr1=Flrl_rofsub, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_roftdo', fldptr1=Flrl_roftdo, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofi', fldptr1=Flrl_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Flrl_irrig', rc=rc)) then
       lnd_irrig = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_irrig', fldptr1=Flrl_irrig, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       lnd_irrig = .false.
    end if

    if ( fldchk(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofl_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofi_wiso', rc=rc)) then
       lnd_wiso_recv = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap_wiso', fldptr2=Fall_evap_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofl_wiso', fldptr2=Flrl_rofl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofi_wiso', fldptr2=Flrl_rofi_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    ic = c_lnd_lrecv
    do n = 1, size(Fall_swnet)
       wgt = is_local%wrap%mesh_info(complnd)%areas(n) * lfrac(n)

       budget_local(f_area      ,ic,ip) = budget_local(f_area      ,ic,ip) + wgt
       budget_local(f_heat_swnet,ic,ip) = budget_local(f_heat_swnet,ic,ip) + wgt*Fall_swnet(n)
       budget_local(f_heat_lwup ,ic,ip) = budget_local(f_heat_lwup ,ic,ip) + wgt*Fall_lwup(n)
       budget_local(f_heat_latv ,ic,ip) = budget_local(f_heat_latv ,ic,ip) + wgt*Fall_lat(n)
       budget_local(f_heat_sen  ,ic,ip) = budget_local(f_heat_sen  ,ic,ip) + wgt*Fall_sen(n)
       budget_local(f_watr_evap ,ic,ip) = budget_local(f_watr_evap ,ic,ip) + wgt*Fall_evap(n)

       budget_local(f_watr_roff ,ic,ip) = budget_local(f_watr_roff ,ic,ip) &
            - wgt*Flrl_rofsur(n) - wgt*Flrl_rofgwl(n) - wgt*Flrl_rofsub(n) - wgt*Flrl_rofdto(n)
       if (lnd_irrig) then
          budget_local(f_watr_roff,ic,ip) = budget_local(f_watr_roff,ic,ip) - wgt*Flrl_irrig(n)
       end if
       budget_local(f_watr_ioff ,ic,ip) = budget_local(f_watr_ioff,ic,ip)  - wgt*Flrl_rofi(n)

       if ( lnd_wiso_recv ) then
          budget_local(f_watr_evap_16O,ic,ip) = budget_local(f_watr_evap_16O,ic,ip) + wgt*Fall_evap_wiso(1,n)
          budget_local(f_watr_evap_18O,ic,ip) = budget_local(f_watr_evap_18O,ic,ip) + wgt*Fall_evap_wiso(2,n)
          budget_local(f_watr_evap_HDO,ic,ip) = budget_local(f_watr_evap_HDO,ic,ip) + wgt*Fall_evap_wiso(3,n)

          budget_local(f_watr_roff_16O,ic,ip) = budget_local(f_watr_roff_16O,ic,ip) - wgt*Flrl_rofl_wiso(1,n)
          budget_local(f_watr_roff_18O,ic,ip) = budget_local(f_watr_roff_18O,ic,ip) - wgt*Flrl_rofl_wiso(2,n)
          budget_local(f_watr_roff_HDO,ic,ip) = budget_local(f_watr_roff_HDO,ic,ip) - wgt*Flrl_rofl_wiso(3,n)

          budget_local(f_watr_ioff_16O,ic,ip) = budget_local(f_watr_ioff_16O,ic,ip) - wgt*Flrl_rofi_wiso(1,n)
          budget_local(f_watr_ioff_18O,ic,ip) = budget_local(f_watr_ioff_18O,ic,ip) - wgt*Flrl_rofi_wiso(2,n)
          budget_local(f_watr_ioff_HDO,ic,ip) = budget_local(f_watr_ioff_HDO,ic,ip) - wgt*Flrl_rofi_wiso(3,n)
       end if
    end do
    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice

    !-------------------------------
    ! to land from mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_lwdn', fldptr1=Faxa_lwdn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainc', fldptr1=Faxa_rainc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainl', fldptr1=Faxa_rainl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainc', fldptr1=Faxa_snowc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainl', fldptr1=Faxa_snowl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Flrl_flood', fldptr1=Flrl_flood, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBExp(complnd), 'Faxa_rainl_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBExp(complnd), 'Faxa_snowl_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBExp(complnd), 'Flrl_flood_wiso', rc=rc)) then
       lnd_wiso_send = .true.
       call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainl_wiso', fldptr2=Faxa_rainl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainc_wiso', fldptr2=Faxa_rainc_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_snowl_wiso', fldptr2=Faxa_snowl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_snowc_wiso', fldptr2=Faxa_snowc_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Flrl_flood_wiso', fldptr2=Flrl_flood_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    ic = c_lnd_lsend
    do n = 1,size(Faxa_lwdn)
       wgt = is_local%wrap%mesh_info(complnd)%areas(n) * lfrac(n)

       budget_local(f_area     ,ic,ip) = budget_local(f_area     ,ic,ip) + wgt
       budget_local(f_heat_lwdn,ic,ip) = budget_local(f_heat_lwdn,ic,ip) + wgt*Faxa_lwdn(n)
       budget_local(f_watr_rain,ic,ip) = budget_local(f_watr_rain,ic,ip) + wgt*Faxa_rainc(n) + wgt*Faxa_rainl(n)
       budget_local(f_watr_snow,ic,ip) = budget_local(f_watr_snow,ic,ip) + wgt*Faxa_snowc(n) + wgt*Faxa_snowl(n)
       budget_local(f_watr_roff,ic,ip) = budget_local(f_watr_roff,ic,ip) - wgt*Flrr_flood(n)

       if ( lnd_wiso_send )then
          budget_local(f_watr_rain_16O,ic,ip) = budget_local(f_watr_rain_16O,ic,ip) + wgt*Faxa_rainc_wiso(1,n) + wgt*Faxa_rainl_wiso(1,n)
          budget_local(f_watr_rain_18O,ic,ip) = budget_local(f_watr_rain_18O,ic,ip) + wgt*Faxa_rainc_wiso(2,n) + wgt*Faxa_rainl_wiso(2,n)
          budget_local(f_watr_rain_HDO,ic,ip) = budget_local(f_watr_rain_HDO,ic,ip) + wgt*Faxa_rainc_wiso(3,n) + wgt*Faxa_rainl_wiso(3,n)

          budget_local(f_watr_snow_16O,ic,ip) = budget_local(f_watr_snow_16O,ic,ip) + wgt*Faxa_snowc_wiso(1,n) + wgt*Faxa_snowl_wiso(1,n)
          budget_local(f_watr_snow_18O,ic,ip) = budget_local(f_watr_snow_18O,ic,ip) + wgt*Faxa_snowc_wiso(2,n) + wgt*Faxa_snowl_wiso(2,n)
          budget_local(f_watr_snow_HDO,ic,ip) = budget_local(f_watr_snow_HDO,ic,ip) + wgt*Faxa_snowc_wiso(3,n) + wgt*Faxa_snowl_wiso(3,n)

          budget_local(f_watr_roff_16O,ic,ip) = budget_local(f_watr_roff_16O,ic,ip) - wgt*Flrr_flood_wiso(1,n)
          budget_local(f_watr_roff_18O,ic,ip) = budget_local(f_watr_roff_18O,ic,ip) - wgt*Flrr_flood_wiso(2,n)
          budget_local(f_watr_roff_HDO,ic,ip) = budget_local(f_watr_roff_HDO,ic,ip) - wgt*Flrr_flood_wiso(3,n)
       end if

    end do
    budget_local(f_heat_latf,ic,ip) = -budget_local(f_watr_snow,ic,ip)*shr_const_latice

  end subroutine med_phases_diag_lnd

  !===============================================================================

  subroutine med_phases_diag_rof( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global river input/output
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(in)           :: gcomp
    integer             , intent(out) :: rc

    ! local variables
    type(InternalState)    :: is_local
    real(r8), pointer      :: lfrac(:)
    integer                :: n,ic,ip
    real(r8)               :: wgt
    logical                :: rof_wiso = .false.
    logical                :: rof_irrig = .false.
    real(r8), pointer      :: Flrl_rofsur(:)       ! to river from mediator
    real(r8), pointer      :: Flrl_rofgwl(:)       ! to river from mediator
    real(r8), pointer      :: Flrl_rofsub(:)       ! to river from mediator
    real(r8), pointer      :: Flrl_roftdo(:)       ! to river from mediator
    real(r8), pointer      :: Flrl_rofi(:)         ! to river from mediator
    real(r8), pointer      :: Flrl_irrig(:)        ! to river from mediator
    real(r8), pointer      :: Flrl_rofl_wiso(:,:)  ! to river from mediator
    real(r8), pointer      :: Flrl_rofi_wiso(:,:)  ! to river from mediator
    real(r8), pointer      :: Forr_rofl(:)         ! from river to mediator
    real(r8), pointer      :: Forr_rofi(:)         ! from river to mediator
    real(r8), pointer      :: Flrr_flood(:)        ! from river to mediator
    real(r8), pointer      :: Firr_rofi(:)         ! from river to mediator
    real(r8), pointer      :: Forr_rofl_wiso(:,:)  ! from river to mediator
    real(r8), pointer      :: Forr_rofi_wiso(:,:)  ! from river to mediator
    real(r8), pointer      :: Flrr_flood_wiso(:,:) ! from river to mediator
    character(*), parameter :: subName = '(med_phases_diag_rof) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------
    ! to river from mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofsur', fldptr1=Flrl_rofsur, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofgwl', fldptr1=Flrl_rofgwl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofsub', fldptr1=Flrl_rofsub, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofdto', fldptr1=Flrl_rofdto, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofi', fldptr1=Flrl_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (fldchk(is_local%wrap%FBImp(comprof,comprof), 'Flrl_irrig', rc=rc)) then
       rof_irrig = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_irrig', fldptr1=Flrl_irrig, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if ( fldchk(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofl_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofi_wiso', rc=rc)) then
       rof_wiso = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofl_wiso', fldptr2=Flrl_rofl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofi_wiso', fldptr2=Flrl_rofi_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    ic = c_rof_rrecv
    do n = 1,size(Flrl_rofdto)
       wgt = is_local%wrap%mesh_info(comprof)%areas(n)

       budget_local(f_watr_roff,ic,ip) = budget_local(f_watr_roff,ic,ip) &
            + wgt*Flrl_rofsur(n) + wgt*Flrl_rofgwl(n) + wgt*Flrl_rofsub(n) + wgt*Flrl_rofdto(n)
       if (rof_irrig /= 0) then
          budget_local(f_watr_roff,ic,ip) = budget_local(f_watr_roff,ic,ip) + wgt*Flrl_irrig(n)
       end if
       budget_local(f_watr_ioff,ic,ip) = budget_local(f_watr_ioff,ic,ip) + wgt*Flrl_rofi(n)

       if ( rof_wiso )then
          budget_local(f_watr_roff_16O,ic,ip) = budget_local(f_watr_roff_16O,ic,ip) + wgt*Flrl_rofl_wiso(1,n)
          budget_local(f_watr_roff_18O,ic,ip) = budget_local(f_watr_roff_18O,ic,ip) + wgt*Flrl_rofl_wiso(2,n)
          budget_local(f_watr_roff_HDO,ic,ip) = budget_local(f_watr_roff_HDO,ic,ip) + wgt*Flrl_rofl_wiso(3,n)

          budget_local(f_watr_ioff_16O,ic,ip) = budget_local(f_watr_ioff_16O,ic,ip) + wgt*Flrl_rofi_wiso(1,n)
          budget_local(f_watr_ioff_18O,ic,ip) = budget_local(f_watr_ioff_18O,ic,ip) + wgt*Flrl_rofi_wiso(2,n)
          budget_local(f_watr_ioff_HDO,ic,ip) = budget_local(f_watr_ioff_HDO,ic,ip) + wgt*Flrl_rofi_wiso(3,n)
       end if
    end do
    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice

    !-------------------------------
    ! from river to mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Forr_rofl', fldptr1=Forr_rofl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Flrl_flood', fldptr1=Flrl_flood, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Forr_rofi', fldptr1=Forr_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Firr_rofi', fldptr1=Firr_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBExp(comprof), 'Forr_rofl_wiso' , rc=rc) .and. &
         fldchk(is_local%wrap%FBExp(comprof), 'Forr_rofi_wiso' , rc=rc) .and. &
         fldchk(is_local%wrap%FBExp(comprof), 'Flrl_flood_wiso', rc=rc)) then
       rof_wiso = .true.
       call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Forr_rofl_wiso', fldptr2=Forr_rofl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Forr_rofi_wiso', fldptr2=Forr_rofi_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Flrl_flood_wiso', fldptr2=Flrl_flood_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    ic = c_rof_rsend
    do n = 1,size(Flrr_flood)
       wgt = is_local%wrap%mesh_info(comprof)%areas(n)

       budget_local(f_watr_roff,ic,ip) = budget_local(f_watr_roff,ic,ip) - wgt*Forr_rofl(n) + wgt*Flrr_flood(n)
       budget_local(f_watr_ioff,ic,ip) = budget_local(f_watr_ioff,ic,ip) - wgt*Forr_rofi(n) - wgt*Firr_rofi(n)

       if ( rof_wiso )then
          budget_local(f_watr_roff_16O,ic,ip) = budget_local(f_watr_roff_16O,ic,ip) - wgt*Forr_rofl_wiso(1,n)
          budget_local(f_watr_roff_18O,ic,ip) = budget_local(f_watr_roff_18O,ic,ip) - wgt*Forr_rofl_wiso(2,n)
          budget_local(f_watr_roff_HDO,ic,ip) = budget_local(f_watr_roff_HDO,ic,ip) - wgt*Forr_rofl_wiso(3,n)

          budget_local(f_watr_ioff_16O,ic,ip) = budget_local(f_watr_ioff_16O,ic,ip) - wgt*Forr_rofi_wiso(1,n)
          budget_local(f_watr_ioff_18O,ic,ip) = budget_local(f_watr_ioff_18O,ic,ip) - wgt*Forr_rofi_wiso(2,n)
          budget_local(f_watr_ioff_HDO,ic,ip) = budget_local(f_watr_ioff_HDO,ic,ip) - wgt*Forr_rofi_wiso(3,n)

          budget_local(f_watr_roff_16O,ic,ip) = budget_local(f_watr_roff_16O,ic,ip) + wgt*Flrr_flood_wiso(1,n)
          budget_local(f_watr_roff_18O,ic,ip) = budget_local(f_watr_roff_18O,ic,ip) + wgt*Flrr_flood_wiso(2,n)
          budget_local(f_watr_roff_HDO,ic,ip) = budget_local(f_watr_roff_HDO,ic,ip) + wgt*Flrr_flood_wiso(3,n)
       end if
    end do
    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice

  end subroutine med_phases_diag_rof

  !===============================================================================

  subroutine med_phases_diag_glc( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global glc output
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(in)  :: gcomp
    integer             , intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n,ic,ip
    real(r8)            :: wgt
    real(r8), pointer   :: Fogg_rofl(:) ! from glc to mediator
    real(r8), pointer   :: Fogg_rofi(:) ! from glc to mediator
    real(r8), pointer   :: Figg_rofi(:) ! from glc to mediator
    character(*), parameter :: subName = '(med_phases_diag_glc) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------
    ! from glc to mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBImp(compglc,compglc), 'Fogg_rofl', fldptr1=Fogg_rofl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compglc,compglc), 'Fogg_rofi', fldptr1=Fogg_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compglc,compglc), 'Figg_rofi', fldptr1=Figg_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ip = period_inst
    ic = c_glc_gsend
    do n = 1,size(Fogg_rofl)
       wgt = is_local%wrap%mesh_info(compglc)%areas(n)
       budget_local(f_watr_roff,ic,ip) = budget_local(f_watr_roff,ic,ip) - wgt*Fogg_rofl(n)
       budget_local(f_watr_ioff,ic,ip) = budget_local(f_watr_ioff,ic,ip) - wgt*Fogg_rofi(n) - wgt*Figg_rofi(n)
    end do
    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice

  end subroutine med_phases_diag_glc

  !===============================================================================

  subroutine med_phases_diag_ocn( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global ocn input from mediator
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(in)  :: gcomp
    integer             , intent(out) :: rc

    ! local variables
    type(InternalState)    :: is_local
    integer                :: n,ic,ip
    real(r8)               :: wgt_i,wgt_o
    logical                :: ocn_wiso = .false.
    character(*), parameter :: subName = '(med_phases_diag_ocn) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------
    ! from ocn to mediator
    !-------------------------------

    call FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac', fldptr1=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac', fldptr1=ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ip = period_inst
    ic = c_ocn_orecv
    do n = 1,size(Fioo_q)
       wgt_o = is_local%wrap%mesh_info(compocn)%areas(n) * ofrac(n)
       wgt_i = is_local%wrap%mesh_info(compocn)%areas(n) * ifrac(n)

       budget_local(f_area,ic,ip) = budget_local(f_area,ic,ip) + wgt_o
       budget_local(f_heat_frz,ic,ip) = budget_local(f_heat_frz,ic,ip) + (wgt_o + wgt_i)*max(0.0_r8,Fioo_q(n))
    end do
    budget_local(f_watr_frz,ic,ip) = budget_local(f_heat_frz,ic,ip) * HFLXtoWFLX

    if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup', rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat' , rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_sen' , rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc)) then

       call FB_getFldPtr(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup', fldptr1=Faox_lwup, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBMed_aoflux_o, 'Faox_lat', fldptr1=Faox_lat, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBMed_aoflux_o, 'Faox_sen', fldptr1=Faox_sen, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', fldptr1=Faox_evap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ip = period_inst
       ic = c_ocn_orecv
       do n = 1,size(Faox_lwup)
          wgt_o = is_local%wrap%mesh_info(compocn)%areas(n) * ofrac(n)

          budget_local(f_heat_lwup,ic,ip) = budget_local(f_heat_lwup,ic,ip) + wgt_o*Faox_lwup(n)
          budget_local(f_heat_latv,ic,ip) = budget_local(f_heat_latv,ic,ip) + wgt_o*Faox_lat(n)
          budget_local(f_heat_sen ,ic,ip) = budget_local(f_heat_sen ,ic,ip) + wgt_o*Faox_sen(n)
          budget_local(f_watr_evap,ic,ip) = budget_local(f_watr_evap,ic,ip) + wgt_o*Faox_evap(n)

          if ( ocn_wiso) then
             budget_local(f_watr_evap_16O,ic,ip) = budget_local(f_watr_evap_16O,ic,ip) + wgt_o*Faox_evap_wiso(1,n)
             budget_local(f_watr_evap_18O,ic,ip) = budget_local(f_watr_evap_18O,ic,ip) + wgt_o*Faox_evap_wiso(2,n)
             budget_local(f_watr_evap_HDO,ic,ip) = budget_local(f_watr_evap_HDO,ic,ip) + wgt_o*Faox_evap_wiso(3,n)
          end if
       end do
    end if

    !-------------------------------
    ! from mediator to ocn
    !-------------------------------

    call FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac', fldptr1=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac', fldptr1=ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Foxx_lwup', rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Foxx_lat' , rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Foxx_sen' , rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Foxx_evap', rc=rc)) then

       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_lwup', fldptr1=Foxx_lwup, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_lat', fldptr1=Foxx_lat, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_sen', fldptr1=Foxx_sen, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_evap', fldptr1=Foxx_evap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ip = period_inst
       ic = c_ocn_orecv
       do n = 1,size(Foxx_lwup)
          wgt_o = is_local%wrap%mesh_info(compocn)%areas(n) * ofrac(n)
          wgt_i = is_local%wrap%mesh_info(compocn)%areas(n) * ifrac(n)

          budget_local(f_heat_lwup,ic,ip) = budget_local(f_heat_lwup,ic,ip) + (wgt_o+wgt_i)*Foxx_lwup(n)
          budget_local(f_heat_latv,ic,ip) = budget_local(f_heat_latv,ic,ip) + (wgt_o+wgt_i)*Foxx_lat(n)
          budget_local(f_heat_sen ,ic,ip) = budget_local(f_heat_sen ,ic,ip) + (wgt_o+wgt_i)*Foxx_sen(n)
          budget_local(f_watr_evap,ic,ip) = budget_local(f_watr_evap,ic,ip) + (wgt_o+wgt_i)*Foxx_evap(n)
       end do
    endif

    call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Fioi_meltw', fldptr1=Fioi_meltw, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Fioi_melth', fldptr1=Fioi_melth, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBExp(compocn), 'Fioi_bergw', rc=rc)) then
       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Fioi_bergw', fldptr1=Fioi_bergw, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ice_bergw = .true.
    end if
    if ( fldchk(is_local%wrap%FBExp(compocn), 'Fioi_bergh', rc=rc)) then
       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Fioi_bergh', fldptr1=Fioi_bergh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ice_bergh = .true.
    end if

    ip = period_inst
    ic = c_ocn_osend
    do n = 1,size(Fioi_meltw)
       wgt_o = is_local%wrap%mesh_info(compocn)%areas(n) * ofrac(n)
       wgt_i = is_local%wrap%mesh_info(compocn)%areas(n) * ifrac(n)

       budget_local(f_area,ic,ip) = budget_local(f_area,ic,ip) + wgt_o

       if (.not. ice_bergw) then
          budget_local(f_watr_melt,ic,ip) = budget_local(nf,ic,ip) + (wgt_o+wgt_i)*Fioi_meltw(n)
       else
          budget_local(f_watr_melt,ic,ip) = budget_local(nf,ic,ip) + (wgt_o+wgt_i)*(Fioi_meltw(n)+Fioi_bergw(n))
       endif

       if (.not. ice_bergh) then
          budget_local(f_heat_melt,ic,ip) = budget_local(f_heat_melt,ic,ip) + (wgt_o+wgt_i)*Fioi_melth(n)
       else
          budget_local(f_heat_melt,ic,ip) = budget_local(f_heat_melt,ic,ip) + (wgt_o+wgt_i)*(Fioi_melth(n) + Fioi_bergh(n))
       endif

       budget_local(f_watr_salt ,ic,ip) = budget_local(f_watr_salt ,ic,ip) + (wgt_o+wgt_i)*Fioi_salt(n) * SFLXtoWFLX
       budget_local(f_heat_swnet,ic,ip) = budget_local(f_heat_swnet,ic,ip) + (wgt_o+wgt_i)*Foxx_swnet(n)
       budget_local(f_heat_lwdn ,ic,ip) = budget_local(f_heat_lwdn ,ic,ip) + (wgt_o+wgt_i)*Faxa_lwdn(n)
       budget_local(f_watr_rain ,ic,ip) = budget_local(f_watr_rain ,ic,ip) + (wgt_o+wgt_i)*Faxa_rain(n)
       budget_local(f_watr_snow ,ic,ip) = budget_local(f_watr_snow ,ic,ip) + (wgt_o+wgt_i)*Faxa_snow(n)
       budget_local(f_watr_roff ,ic,ip) = budget_local(f_watr_roff ,ic,ip) + (wgt_o+wgt_i)*Foxx_rofl(n)
       budget_local(f_watr_ioff ,ic,ip) = budget_local(f_watr_ioff ,ic,ip) + (wgt_o+wgt_i)*Foxx_rofi(n)

       if ( ocn_wiso )then
          budget_local(f_watr_melt_16O,ic,ip) = budget_local(f_watr_melt_16O,ic,ip) + (wgt_o+wgt_i)*Fioi_meltw_wiso(1,n)
          budget_local(f_watr_melt_18O,ic,ip) = budget_local(f_watr_melt_18O,ic,ip) + (wgt_o+wgt_i)*Fioi_meltw_wiso(2,n)
          budget_local(f_watr_melt_HDO,ic,ip) = budget_local(f_watr_melt_HDO,ic,ip) + (wgt_o+wgt_i)*Fioi_meltw_wiso(3,n)

          budget_local(f_watr_rain_16O,ic,ip) = budget_local(f_watr_rain_16O,ic,ip) + (wgt_o+wgt_i)*Faxa_rain_wiso(1,n)
          budget_local(f_watr_rain_18O,ic,ip) = budget_local(f_watr_rain_18O,ic,ip) + (wgt_o+wgt_i)*Faxa_rain_wiso(2,n)
          budget_local(f_watr_rain_HDO,ic,ip) = budget_local(f_watr_rain_HDO,ic,ip) + (wgt_o+wgt_i)*Faxa_rain_wiso(3,n)

          budget_local(f_watr_snow_16O,ic,ip) = budget_local(f_watr_snow_16O,ic,ip) + (wgt_o+wgt_i)*Faxa_snow_wiso(1,n)
          budget_local(f_watr_snow_18O,ic,ip) = budget_local(f_watr_snow_18O,ic,ip) + (wgt_o+wgt_i)*Faxa_snow_wiso(2,n)
          budget_local(f_watr_snow_HDO,ic,ip) = budget_local(f_watr_snow_HDO,ic,ip) + (wgt_o+wgt_i)*Faxa_snow_wiso(3,n)

          budget_local(f_watr_roff_16O,ic,ip) = budget_local(f_watr_roff_16O,ic,ip) + (wgt_o+wgt_i)*Foxx_rofl_wiso(1,n)
          budget_local(f_watr_roff_18O,ic,ip) = budget_local(f_watr_roff_18O,ic,ip) + (wgt_o+wgt_i)*Foxx_rofl_wiso(2,n)
          budget_local(f_watr_roff_HDO,ic,ip) = budget_local(f_watr_roff_HDO,ic,ip) + (wgt_o+wgt_i)*Foxx_rofl_wiso(3,n)

          budget_local(f_watr_ioff_16O,ic,ip) = budget_local(f_watr_ioff_16O,ic,ip) + (wgt_o+wgt_i)*Foxx_rofi_wiso(1,n)
          budget_local(f_watr_ioff_18O,ic,ip) = budget_local(f_watr_ioff_18O,ic,ip) + (wgt_o+wgt_i)*Foxx_rofi_wiso(2,n)
          budget_local(f_watr_ioff_HDO,ic,ip) = budget_local(f_watr_ioff_HDO,ic,ip) + (wgt_o+wgt_i)*Foxx_rofi_wiso(3,n)
       end if
    end do
    budget_local(f_heat_latf,ic,ip) = -budget_local(f_watr_snow,ic,ip)*shr_const_latice
    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice

    ! EBK -- isotope r2x_Forr_rofl/i?

  end subroutine med_phases_diag_ocn_med2ocn

  !===============================================================================

  subroutine med_phases_diag_ice_ice2med( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global ice input/output flux diagnostics
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(in)  :: gcomp
    integer             , intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n,ic,ip
    real(r8)            :: wgt_i, wgt_o
    real(r8), pointer   :: ofrac(:)
    real(r8), pointer   :: ifrac(:)
    real(r8), pointer   :: Fioi_melth(:)
    real(r8), pointer   :: Fioi_meltw(:)
    real(r8), pointer   :: Fioi_salt(:)
    real(r8), pointer   :: Fioi_swpen(:)
    real(r8), pointer   :: Faii_swnet(:)
    real(r8), pointer   :: Faii_lwup(:)
    real(r8), pointer   :: Faii_sen(:)
    real(r8), pointer   :: Faii_lat(:)
    real(r8), pointer   :: Faii_evap(:)
    real(r8), pointer   :: Fioi_meltw_wiso(:,:)
    real(r8), pointer   :: Faii_evap_wiso(:,:)
    logical             :: ice_wiso = .false.
    character(*), parameter :: subName = '(med_phases_diag_ice_ice2med) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ifrac', fldptr1=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ofrac', fldptr1=ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Fioi_melth', fldptr1=Fioi_melth, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Fioi_meltw', fldptr1=Fioi_meltw, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Fioi_salt', fldptr1=Fioi_salt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Fioi_swpen', fldptr1=Fioi_swpen, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_swnet', fldptr1=Faii_swnet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_lwup', fldptr1=Faii_lwup, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_sen', fldptr1=Faii_sen, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_lat', fldptr1=Faii_lat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_evap', fldptr1=Faii_evap, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_meltw_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap_wiso' , rc=rc)) then
       ice_wiso = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Fioi_meltw_wiso', fldptr2=Fioi_meltw_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_evap_wiso', fldptr2=Faii_evap_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    do n = 1,size(Fioi_melth)
       wgt_o = is_local%wrap%mesh_info(compice)%areas(n) * ofrac(n)
       wgt_i = is_local%wrap%mesh_info(compice)%areas(n) * ifrac(n)

       if (is_local%wrap%mesh_info(compice)%lats(n) > 0.0_r8) then
          ic = c_inh_irecv
       else
          ic = c_ish_irecv
       endif

       budget_local(f_area  ,ic,ip) = budget_local(f_area  ,ic,ip) + wgt_i
       budget_local(f_heat_melt ,ic,ip) = budget_local(f_heat_melt ,ic,ip) - wgt_i*Fioi_melth(n)
       budget_local(f_watr_melt ,ic,ip) = budget_local(f_watr_melt ,ic,ip) - wgt_i*Fioi_meltw(n)
       budget_local(f_watr_salt ,ic,ip) = budget_local(f_watr_salt ,ic,ip) - wgt_i*Fioi_salt(n) * SFLXtoWFLX
       budget_local(f_heat_swnet,ic,ip) = budget_local(f_heat_swnet,ic,ip) + wgt_i*Faii_swnet(n) - wgt_i*Fioi_swpen(n)
       budget_local(f_heat_lwup ,ic,ip) = budget_local(f_heat_lwup ,ic,ip) + wgt_i*Faii_lwup(n)
       budget_local(f_heat_latv ,ic,ip) = budget_local(f_heat_latv ,ic,ip) + wgt_i*Faii_lat(n)
       budget_local(f_heat_sen  ,ic,ip) = budget_local(f_heat_sen  ,ic,ip) + wgt_i*Faii_sen(n)
       budget_local(f_watr_evap ,ic,ip) = budget_local(f_watr_evap ,ic,ip) + wgt_i*Faii_evap(n)

       if ( ice_wiso )then
          budget_local(f_watr_melt_16O,ic,ip) = budget_local(f_watr_melt_16O,ic,ip) - wgt_i*Fioi_meltw_wiso(1,n)
          budget_local(f_watr_melt_18O,ic,ip) = budget_local(f_watr_melt_18O,ic,ip) - wgt_i*Fioi_meltw_wiso(2,n)
          budget_local(f_watr_melt_HDO,ic,ip) = budget_local(f_watr_melt_HDO,ic,ip) - wgt_i*Fioi_meltw_wiso(3,n)

          budget_local(f_watr_evap_16O,ic,ip) = budget_local(f_watr_evap_16O,ic,ip) + wgt_i*Faii_evap_wiso(1,n)
          budget_local(f_watr_evap_18O,ic,ip) = budget_local(f_watr_evap_18O,ic,ip) + wgt_i*Faii_evap_wiso(2,n)
          budget_local(f_watr_evap_HDO,ic,ip) = budget_local(f_watr_evap_HDO,ic,ip) + wgt_i*Faii_evap_wiso(3,n)
       end if
    end do

  end subroutine med_phases_diag_ice_ice2med

  !===============================================================================

  subroutine med_phases_diag_ice_med2ice( gcomp, rc)

    ! ------------------------------------------------------------------
    ! Compute global ice input/output flux diagnostics
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(in)           :: gcomp
    integer             , intent(out)          :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n,ic,ip
    real(r8)            :: wgt_i, wgt_o
    real(r8), pointer   :: ofrac(:)
    real(r8), pointer   :: ifrac(:)
    real(r8), pointer   :: Faxa_lwdn(:)
    real(r8), pointer   :: Faxa_rain(:)
    real(r8), pointer   :: Faxa_snow(:)
    real(r8), pointer   :: Fixx_rofi(:)
    real(r8), pointer   :: Fioo_q(:)
    real(r8), pointer   :: Faxa_rain_wiso(:,:)
    real(r8), pointer   :: Faxa_snow_wiso(:,:)
    logical             :: ice_wiso = .false.
    character(*), parameter :: subName = '(med_phases_diag_ice_med2ice) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ifrac', fldptr1=ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ofrac', fldptr1=ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBExp(compice), 'Faxa_rain_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBExp(compice), 'Faxa_snow_wiso', rc=rc)) then
       ice_wiso = .true.
       call FB_getFldPtr(is_local%wrap%FBExp(compice), 'Faxa_rain_wiso', fldptr2=Faxa_rain_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBExp(compice), 'Faxa_snow_wiso', fldptr2=Faxa_snow_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    do n = 1,size(Faxa_lwdn)
       wgt_o = is_local%wrap%mesh_info(compice)%areas(n) * ofrac(n)
       wgt_i = is_local%wrap%mesh_info(compice)%areas(n) * ifrac(n)

       if (is_local%wrap%mesh_info(compice)%lats(n) > 0.0_r8) then
          ic = c_inh_isend
       else
          ic = c_ish_isend
       endif

       budget_local(f_area ,ic,ip) = budget_local(f_area ,ic,ip) + wgt_i
       budget_local(f_heat_lwdn,ic,ip) = budget_local(f_heat_lwdn,ic,ip) + wgt_i*Faxa_lwdn(n)
       budget_local(f_watr_rain,ic,ip) = budget_local(f_watr_rain,ic,ip) + wgt_i*Faxa_rain(n)
       budget_local(f_watr_snow,ic,ip) = budget_local(f_watr_snow,ic,ip) + wgt_i*Faxa_snow(n)
       budget_local(f_watr_ioff,ic,ip) = budget_local(f_watr_ioff,ic,ip) + wgt_i*Fixx_rofi(n)
       budget_local(f_heat_frz ,ic,ip) = budget_local(f_heat_frz ,ic,ip) - (wgt_o + wgt_i)*max(0.0_r8,Fioo_q(n))

       if ( ice_wiso ) then
          budget_local(f_watr_rain_16O,ic,ip) = budget_local(f_watr_rain_16O,ic,ip) + wgt_i*Faxa_rain_wiso(1,n)
          budget_local(f_watr_rain_18O,ic,ip) = budget_local(f_watr_rain_18O,ic,ip) + wgt_i*Faxa_rain_wiso(2,n)
          budget_local(f_watr_rain_HDO,ic,ip) = budget_local(f_watr_rain_HDO,ic,ip) + wgt_i*Faxa_rain_wiso(3,n)

          budget_local(f_watr_snow_16O,ic,ip) = budget_local(f_watr_snow_16O,ic,ip) + wgt_i*Faxa_snow_wiso(1,n)
          budget_local(f_watr_snow_18O,ic,ip) = budget_local(f_watr_snow_18O,ic,ip) + wgt_i*Faxa_snow_wiso(2,n)
          budget_local(f_watr_snow_HDO,ic,ip) = budget_local(f_watr_snow_HDO,ic,ip) + wgt_i*Faxa_snow_wiso(3,n)
       end if
    end do

    ic = c_inh_isend
    budget_local(f_heat_latf,ic,ip) = -budget_local(f_watr_snow,ic,ip)*shr_const_latice
    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice
    budget_local(f_watr_frz ,ic,ip) =  budget_local(f_heat_frz ,ic,ip)*HFLXtoWFLX

    ic = c_ish_isend
    budget_local(f_heat_latf,ic,ip) = -budget_local(f_watr_snow,ic,ip)*shr_const_latice
    budget_local(f_heat_ioff,ic,ip) = -budget_local(f_watr_ioff,ic,ip)*shr_const_latice
    budget_local(f_watr_frz ,ic,ip) =  budget_local(f_heat_frz ,ic,ip)*HFLXtoWFLX

  end subroutine med_phases_diag_ice_med2ice

  !===============================================================================

  subroutine med_phases_diag_print(gcomp, rc)

    ! ------------------------------------------------------------------
    ! Print global budget diagnostics.
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(in) :: gcomp
    integer             , intent(out):: rc

    ! local variables
    type(ESMF_Clock)        :: clock
    type(ESMF_Alarm)        :: stop_alarm
    integer                 :: cdate,sec    ! coded date, seconds
    integer                 :: yr,mon,day   ! date
    integer                 :: output_level ! print level
    logical                 :: sumdone      ! has a sum been computed yet
    character(CS)           :: cvalue
    logical                 :: first_call = .true.
    character(*), parameter :: subName = '(med_phases_diag_print) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !-------------------------------------------------------------------------------
    ! Get config variables on first call
    !-------------------------------------------------------------------------------

    if (first_call) then
       call NUOPC_CompAttributeGet(gcomp, name='budget_inst', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) budget_print_inst

       call NUOPC_CompAttributeGet(gcomp, name='budget_daily', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) budget_print_daily

       call NUOPC_CompAttributeGet(gcomp, name='budget_month', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) budget_print_month

       call NUOPC_CompAttributeGet(gcomp, name='budget_ann', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) budget_print_ann

       call NUOPC_CompAttributeGet(gcomp, name='budget_ltann', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) budget_print_ltann

       call NUOPC_CompAttributeGet(gcomp, name='budget_ltend', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) budget_print_ltend

       first_call = .false.
    end if

    !-------------------------------------------------------------------------------
    ! Get clock and alarm info
    !-------------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=curr_year, mm=curr_mon, dd=curr_day, s=curr_tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    cdate = curr_yr*10000 + curr_mon*100 + curr_day
    ! TODO: add stop alarm

    !-------------------------------------------------------------------------------
    ! Accumulate budget data
    !-------------------------------------------------------------------------------

    call med_diag_accum()

    !-------------------------------------------------------------------------------
    ! Print budget data
    !-------------------------------------------------------------------------------

    sumdone = .false.

    do ip = 1,period_ntypes

       ! Determine output level for this period type
       output_level = 0
       if (ip == period_inst) then
          output_level = max(output_level, budget_print_inst)
       endif
       if (ip==period_curr_day .and. curr_tod==0) then
          output_level = max(output_level, budget_print_daily)
       endif
       if (ip==period_mon .and. curr_day==1 .and. curr_tod==0) then
          output_level = max(output_level, budget_print_month)
       endif
       if (ip==period_ann .and. mon==1 .and. curr_day==1 .and. curr_tod==0) then
          output_level = max(output_level, budget_print_ann)
       endif
       if (ip==period_inf .and. mon==1 .and. curr_day==1 .and. curr_tod==0) then
          output_level = max(output_level, budget_print_ltann)
       endif
       if (ip==period_inf .and. stop_alarm) then
          output_level = max(output_level, budget_print_ltend)
       endif

       ! Currently output_level is limited to levels of 0,1,2, 3
       ! (see comment for print obtains at top)

       if (output_level > 0 .and. .not. sumdone) then
          ! Some budgets will be printed for this period type

          ! Determine sums if not already done
          call med_diag_sum_master()
          sumdone = .true.
          dataGpr(:,:,:) = budget_global(:,:,:)

          ! budget normalizations (global area and 1e6 for water)
          dataGpr = dataGpr/(4.0_r8*shr_const_pi)
          dataGpr(f_watr_beg:f_watr_end,:,:) = dataGpr(f_watr_beg:f_watr_end,:,:) * 1.0e6_r8
          if ( flds_wiso ) then
             dataGpr(iso0(1):isof(nisotopes),:,:) = dataGpr(iso0(1):isof(nisotopes),:,:) * 1.0e6_r8
          end if
          dataGpr(:,:,:) = dataGpr/budget_counter(:,:,:)
       endif

       if (output_level > 0 .and. masterproc) then
          ! Write diagnostic tables to logunit (masterproc only)
          if (output_level >= 3) then
             ! detail atm budgets and breakdown into components ---
             call med_diag_print_atm(dataGrp, ip)
          end if
          if (output_level >= 2) then
             ! detail lnd/ocn/ice component budgets ----
             call med_diag_print_lnd_ice_ocn(dataGrp, ip)
          end if
          if (output_level >= 1) then
             ! net summary budgets
             call med_diag_print_summary(dataGrp, ip)
          endif
          write(logunit,*) ' '
       endif ! output_level > 0 and masterproc

    enddo  ! ip = 1, period_types

    !-------------------------------------------------------------------------------
    ! Zero budget data
    !-------------------------------------------------------------------------------

    call med_diag_zero(gcomp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine med_phases_diag_print

  !===============================================================================

  subroutine med_diag_print_atm(data)

    ! ---------------------------------------------------------
    ! detail atm budgets and breakdown into components
    ! ---------------------------------------------------------

    ! intput/output variables
    real(r8), intent(in) :: data(:,:,:) ! values to print, scaled and such

    ! local variables
    integer           :: ic,nf,ip,is ! data array indicies
    integer           :: ica,icl
    integer           :: icn,ics,ico
    character(len=40) :: str         ! string
    character(*), parameter:: subName = '(med_phases_diag_print_level3) '
    ! ------------------------------------------------------------------

    do ic = 1,2
       if (ic == 1) then    ! from atm to mediator
          ica = c_atm_arecv ! total from atm
          icl = c_lnd_arecv ! from land   to med on atm grid
          icn = c_inh_arecv ! from ice-nh to med on atm grid
          ics = c_ish_arecv ! from ice-sh to med on atm grid
          ico = c_ocn_arecv ! from ocn to to med on atm grid
          str = "ATM_to_CPL"
       elseif (ic == 2) then ! from mediator to atm
          ica = c_atm_asend  ! merged to atm
          icl = c_lnd_asend  ! from land   to atm
          icn = c_inh_asend  ! from ice-nh to atm
          ics = c_ish_asend  ! from ice-sh to atm
          ico = c_ocn_asend  ! from ocn    to atm
          str = "CPL_TO_ATM"
       else
          call shr_sys_abort(subname//' ERROR in ic index code 411')
       endif

       write(logunit,*) ' '
       write(logunit,FAH) subname,trim(str)//' AREA BUDGET (m2/m2): period = ',trim(period_name(ip)),&
            ': date = ', cdate, curr_tod
       write(logunit,FA0) cdata(ica)%name,cdata(ici)%name,cdata%name(icn),cdata(ics)%name,cdata(ico)%name,' *SUM*  '
       write(logunit,FA1) budget_fields(f_area)%fname,&
            data(f_area,ica,ip), &
            data(f_area,icl,ip), &
            data(f_area,icn,ip), &
            data(f_area,ics,ip), &
            data(f_area,ico,ip), &
            data(f_area,ica,ip) + data(f_area,icl,ip) + &
            data(f_area,icn,ip) + data(f_area,ics,ip) + data(f_area,ico,ip)

       write(logunit,*) ' '
       write(logunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',trim(period_name(ip)),': date = ',cdate,curr_tod
       write(logunit,FA0) cdata(ica)%name,cdata(icl)%name,cdata(icn)%name,cdata(ics)%name,cdata(ico)%name,' *SUM*  '
       do nf = f_heat_beg, f_heat_end
          write(logunit,FA1) fname(nf),&
               data(nf,ica,ip), &
               data(nf,icl,ip), &
               data(nf,icn,ip), &
               data(nf,ics,ip), &
               data(nf,ico,ip), &
               data(nf,ica,ip) + data(nf,icl,ip) + data(nf,icn,ip) + data(nf,ics,ip) + data(nf,ico,ip)
       enddo
       write(logunit,FA1)    '   *SUM*'   ,&
            sum(data(f_heat_beg:f_heat_end,ica,ip)), &
            sum(data(f_heat_beg:f_heat_end,icl,ip)), &
            sum(data(f_heat_beg:f_heat_end,icn,ip)), &
            sum(data(f_heat_beg:f_heat_end,ics,ip)), &
            sum(data(f_heat_beg:f_heat_end,ico,ip)), &
            sum(data(f_heat_beg:f_heat_end,ica,ip)) + sum(data(f_heat_beg:f_heat_end,icl,ip)) + &
            sum(data(f_heat_beg:f_heat_end,icn,ip)) + sum(data(f_heat_beg:f_heat_end,ics,ip)) + &
            sum(data(f_heat_beg:f_heat_end,ico,ip))

       write(logunit,*) ' '
       write(logunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(period_name(ip)),': date = ',cdate,curr_tod
       write(logunit,FA0) cdata(ica)%name,cdata(icl)%name,cdata(icn)%name,cdata(ics)%name,cdata(ico)%name,' *SUM*  '
       do nf = f_watr_beg, f_watr_end
          write(logunit,FA1) fname(nf),&
               data(nf,ica,ip), &
               data(nf,icl,ip), &
               data(nf,icn,ip), &
               data(nf,ics,ip), &
               data(nf,ico,ip), &
               data(nf,ica,ip) + data(nf,icl,ip) + data(nf,icn,ip) + data(nf,ics,ip) + data(nf,ico,ip)
       enddo
       write(logunit,FA1)    '   *SUM*'   ,&
            sum(data(f_watr_beg:f_watr_end,ica,ip)), &
            sum(data(f_watr_beg:f_watr_end,icl,ip)), &
            sum(data(f_watr_beg:f_watr_end,icn,ip)), &
            sum(data(f_watr_beg:f_watr_end,ics,ip)), &
            sum(data(f_watr_beg:f_watr_end,ico,ip)), &
            sum(data(f_watr_beg:f_watr_end,ica,ip)) + sum(data(f_watr_beg:f_watr_end,icl,ip)) + &
            sum(data(f_watr_beg:f_watr_end,icn,ip)) + sum(data(f_watr_beg:f_watr_end,ics,ip)) + &
            sum(data(f_watr_beg:f_watr_end,ico,ip))

       if ( flds_wiso ) then
          do is = 1, nisotopes
             write(logunit,*) ' '
             write(logunit,FAH) subname,trim(str)//' '//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ', &
                  trim(period_name(ip)),': date = ',cdate,curr_tod
             write(logunit,FA0) cdata(ica)%name,cdata(icl)%name,cdata(icn)%name,cdata(ics)%name,cdata(ico)%name,' *SUM*  '
             do nf = iso0(is), isof(is)
                write(logunit,FA1)    fname(nf),&
                     data(nf,ica,ip), &
                     data(nf,icl,ip), &
                     data(nf,icn,ip), &
                     data(nf,ics,ip), &
                     data(nf,ico,ip), &
                     data(nf,ica,ip) + data(nf,icl,ip) + data(nf,icn,ip) + data(nf,ics,ip) + data(nf,ico,ip)
             enddo
             write(logunit,FA1)    '   *SUM*', &
                  sum(data(iso0(is):isof(is),ica,ip)), &
                  sum(data(iso0(is):isof(is),icl,ip)), &
                  sum(data(iso0(is):isof(is),icn,ip)), &
                  sum(data(iso0(is):isof(is),ics,ip)), &
                  sum(data(iso0(is):isof(is),ico,ip)), &
                  sum(data(iso0(is):isof(is),ica,ip)) + sum(data(iso0(is):isof(is),icl,ip)) + &
                  sum(data(iso0(is):isof(is),icn,ip)) + sum(data(iso0(is):isof(is),ics,ip)) + &
                  sum(data(iso0(is):isof(is),ico,ip))
          end do
       end if

    enddo

  end subroutine med_diag_print_atm

  !===============================================================================

  subroutine med_diag_print_lnd_ocn_ice(data, ip)

    ! ---------------------------------------------------------
    ! detail lnd/ocn/ice component budgets
    ! ---------------------------------------------------------

    ! intput/output variables
    real(r8), intent(in) :: data(:,:,:) ! values to print, scaled and such
    integer , intent(in) :: ip

    ! local variables
    integer           :: ic,nf,is ! data array indicies
    integer           :: icar,icas
    integer           :: icxs,icxr
    character(len=40) :: str      ! string
    character(*), parameter :: subName = '(med_diag_print_lnd_ocn_ice) '
    ! ------------------------------------------------------------------

    do ic = 1,4
       if (ic == 1) then
          icar = c_lnd_arecv
          icxs = c_lnd_lsend
          icxr = c_lnd_lrecv
          icas = c_lnd_asend
          str = "LND"
       elseif (ic == 2) then
          icar = c_ocn_arecv
          icxs = c_ocn_osend
          icxr = c_ocn_orecv
          icas = c_ocn_asend
          str = "OCN"
       elseif (ic == 3) then
          icar = c_inh_arecv
          icxs = c_inh_isend
          icxr = c_inh_irecv
          icas = c_inh_asend
          str = "ICE_NH"
       elseif (ic == 4) then
          icar = c_ish_arecv
          icxs = c_ish_isend
          icxr = c_ish_irecv
          icas = c_ish_asend
          str = "ICE_SH"
       else
          call shr_sys_abort(subname//' ERROR in ic index code 412')
       endif

       ! heat budgets atm<->lnd, atm<->ocn, atm<->ice_nh, atm<->ice_sh,

       write(logunit,*) ' '
       write(logunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',trim(period_name(ip)),': date = ',cdate,curr_tod
       write(logunit,FA0) cdata(icar)%name,cdata(icxs)%name,cdata(icxr)%name,cdata(icas)%name,' *SUM*  '
       do nf = f_heat_beg, f_heat_end
          write(logunit,FA1) fname(nf),&
               -data(nf,icar,ip), &
                data(nf,icxs,ip), &
                data(nf,icxr,ip), &
               -data(nf,icas,ip), &
               -data(nf,icar,ip) + data(nf,icxs,ip) + data(nf,icxr,ip) - data(nf,icas,ip)
       enddo
       write(logunit,FA1)'   *SUM*',&
            -sum(data(f_heat_beg:f_heat_end,icar,ip)), &
             sum(data(f_heat_beg:f_heat_end,icxs,ip)), &
             sum(data(f_heat_beg:f_heat_end,icxr,ip)), &
            -sum(data(f_heat_beg:f_heat_end,icas,ip)), &
            -sum(data(f_heat_beg:f_heat_end,icar,ip)) + sum(data(f_heat_beg:f_heat_end,icxs,ip)) + &
             sum(data(f_heat_beg:f_heat_end,icxr,ip)) - sum(data(f_heat_beg:f_heat_end,icas,ip))

       ! water budgets atm<->lnd, atm<->ocn, atm<->ice_nh, atm<->ice_sh,

       write(logunit,*) ' '
       write(logunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(period_name(ip)),': date = ',cdate,curr_tod
       write(logunit,FA0) cdata(icar)%name,cdata(icxs)%name,cdata(icxr)%name,cdata(icas)%name,' *SUM*  '
       do nf = f_watr_beg, f_watr_end
          write(logunit,FA1)  fname(nf),&
               -data(nf,icar,ip),&
                data(nf,icxs,ip), &
                data(nf,icxr,ip),&
               -data(nf,icas,ip), &
               -data(nf,icar,ip) + data(nf,icxs,ip) + data(nf,icxr,ip) - data(nf,icas,ip)
       enddo
       write(logunit,FA1)    '   *SUM*',&
            -sum(data(f_watr_beg:f_watr_end,icar,ip)), &
             sum(data(f_watr_beg:f_watr_end,icxs,ip)), &
             sum(data(f_watr_beg:f_watr_end,icxr,ip)), &
            -sum(data(f_watr_beg:f_watr_end,icas,ip)), &
            -sum(data(f_watr_beg:f_watr_end,icar,ip)) + sum(data(f_watr_beg:f_watr_end,icxs,ip)) + &
             sum(data(f_watr_beg:f_watr_end,icxr,ip)) - sum(data(f_watr_beg:f_watr_end,icas,ip))

       if ( flds_wiso ) then
          do is = 1, nisotopes

             ! heat budgets atm<->lnd, atm<->ocn, atm<->ice_nh, atm<->ice_sh for water isotopes

             write(logunit,*) ' '
             write(logunit,FAH) subname,trim(str)//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(period_name(ip)), &
                  ': date = ',cdate,curr_tod
             write(logunit,FA0) cdata(icar)%name,cdata(icxs)%name,cdata(icxr)%name,cdata(icas)%name,' *SUM*  '
             do nf = iso0(is), isof(is)
                write(logunit,FA1) fname(nf),&
                     -data(nf,icar,ip), &
                      data(nf,icxs,ip), &
                      data(nf,icxr,ip), &
                     -data(nf,icas,ip), &
                     -data(nf,icar,ip) + data(nf,icxs,ip) + data(nf,icxr,ip) - data(nf,icas,ip)
             enddo
             write(logunit,FA1)    '   *SUM*',&
                  -sum(data(iso0(is):isof(is),icar,ip)),&
                   sum(data(iso0(is):isof(is),icxs,ip)), &
                   sum(data(iso0(is):isof(is),icxr,ip)), &
                  -sum(data(iso0(is):isof(is),icas,ip)), &
                  -sum(data(iso0(is):isof(is),icar,ip)) + sum(data(iso0(is):isof(is),icxs,ip)) + &
                   sum(data(iso0(is):isof(is),icxr,ip)) - sum(data(iso0(is):isof(is),icas,ip))

             ! water budgets atm<->lnd, atm<->ocn, atm<->ice_nh, atm<->ice_sh for water isotopes

             write(logunit,*) ' '
             write(logunit,FAH) subname,trim(str)//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(period_name(ip)),&
                  ': date = ',cdate,curr_tod
             write(logunit,FA0) cdata(icar)%name,cdata(icxs)%name,cdata(icxr)%name,cdata(icas)%name,' *SUM*  '
             do nf = iso0(is), isof(is)
                write(logunit,FA1) fname(nf),&
                     -data(nf,icar,ip), &
                      data(nf,icxs,ip), &
                      data(nf,icxr,ip), &
                     -data(nf,icas,ip), &
                     -data(nf,icar,ip) + data(nf,icxs,ip) + data(nf,icxr,ip) - data(nf,icas,ip)
             enddo
             write(logunit,FA1)    '   *SUM*',                &
                  -sum(data(iso0(is):isof(is), icar, ip)), &
                   sum(data(iso0(is):isof(is), icxs, ip)), &
                   sum(data(iso0(is):isof(is), icxr, ip)), &
                  -sum(data(iso0(is):isof(is), icas, ip)), &
                  -sum(data(iso0(is):isof(is), icar, ip)) + sum(data(iso0(is):isof(is), icxs, ip)) + &
                   sum(data(iso0(is):isof(is), icxr, ip)) - sum(data(iso0(is):isof(is), icas, ip))
          end do
       end if
    enddo

  end subroutine med_diag_print_lnd_ocn_ice

  !===============================================================================

  subroutine med_diag_print_net_summary(data, ip)

    ! ---------------------------------------------------------
    ! net summary budgets
    ! ---------------------------------------------------------

    ! intput/output variables
    real(r8), intent(in) :: data(:,:,:) ! values to print, scaled and such
    integer , intent(in) :: ip

    ! local variables
    integer  :: ic,nf,is ! data array indicies
    real(r8) :: atm_area, lnd_area, ocn_area
    real(r8) :: ice_area_nh, ice_area_sh
    real(r8) :: net_water_atm, net_water_lnd, net_water_rof
    real(r8) :: net_water_ocn, net_water_ice_nh, net_water_ice_sh
    real(r8) :: net_water_glc
    real(r8) :: sum_net_water_atm, sum_net_water_lnd, sum_net_water_rof
    real(r8) :: sum_net_water_ocn, sum_net_water_ice_nh, sum_net_water_ice_sh
    real(r8) :: sum_net_water_glc
    real(r8) :: net_heat_atm, net_heat_lnd, net_heat_rof
    real(r8) :: net_heat_ocn,  net_heat_ice_nh, net_heat_ice_sh
    real(r8) :: net_heat_glc
    real(r8) :: sum_net_heat_atm, sum_net_heat_lnd, sum_net_heat_rof
    real(r8) :: sum_net_heat_ocn, sum_net_heat_ice_nh, sum_net_heat_ice_sh
    real(r8) :: sum_net_heat_glc
    character(len=40) :: str
    character(*), parameter:: subName = '(med_diag_print_level3) '
    ! ------------------------------------------------------------------

    ! write out areas

    write(logunit,*) ' '
    write(logunit,FAH) subname,'NET AREA BUDGET (m2/m2): period = ',trim(period_name(ip)),': date = ',cdate,curr_tod
    write(logunit,FA0) '     atm','     lnd','     ocn','  ice nh','  ice sh',' *SUM*  '
    atm_area    = data(f_area,c_atm_arecv,ip)
    lnd_area    = data(f_area,c_lnd_lrecv,ip)
    ocn_area    = data(f_area,c_ocn_orecv,ip)
    ice_area_nh = data(f_area,c_inh_irecv,ip)
    ice_area_sh = data(f_area,c_ish_irecv,ip)
    write(logunit,FA1) fname(f_area), atm_area, lnd_area, ocn_area, ice_area_nh, ice_area_sh

    ! write out net heat budgets

    write(logunit,*) ' '
    write(logunit,FAH) subname,'NET HEAT BUDGET (W/m2): period = ',trim(period_name(ip)),': date = ',cdate,curr_tod
    write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
    do nf = f_heat_beg, f_heat_end
       net_heat_atm    = data(nf, c_atm_arecv, ip) + data(nf, c_atm_asend, ip)
       net_heat_lnd    = data(nf, c_lnd_lrecv, ip) + data(nf, c_lnd_lsend, ip)
       net_heat_rof    = data(nf, c_rof_rrecv, ip) + data(nf, c_rof_rsend, ip)
       net_heat_ocn    = data(nf, c_ocn_orecv, ip) + data(nf, c_ocn_osend, ip)
       net_heat_ice_nh = data(nf, c_inh_irecv, ip) + data(nf, c_inh_isend, ip)
       net_heat_ice_sh = data(nf, c_ish_irecv, ip) + data(nf, c_ish_isend, ip)
       net_heat_glc    = data(nf, c_glc_grecv, ip) + data(nf, c_glc_gsend, ip)
       net_heat_tot    = net_heat_atm + net_heat_lnd + net_heat_rof + net_heat_ocn + &
                         net_heat_ice_nh + net_heat_ice_sh + net_heat_glc

       write(logunit,FA1r)  fname(nf), &
            net_heat_atm, net_heat_lnd, net_heat_rof, net_heat_ocn, &
            net_heat_ice_nh, net_heat_ice_sh, net_heat_glc, net_heat_tot
    end do

    ! Write out sum over all net heat budgets (sum over f_heat_beg -> f_heat_end)

    sum_net_heat_atm    = sum(data(f_heat_beg:f_heat_end, c_atm_arecv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_atm_asend, ip))
    sum_net_heat_lnd    = sum(data(f_heat_beg:f_heat_end, c_lnd_lrecv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_lnd_lsend, ip))
    sum_net_heat_rof    = sum(data(f_heat_beg:f_heat_end, c_rof_rrecv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_rof_rsend, ip))
    sum_net_heat_ocn    = sum(data(f_heat_beg:f_heat_end, c_ocn_orecv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_ocn_osend, ip))
    sum_net_heat_ice_nh = sum(data(f_heat_beg:f_heat_end, c_inh_irecv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_inh_isend, ip))
    sum_net_heat_ice_sh = sum(data(f_heat_beg:f_heat_end, c_ish_irecv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_ish_isend, ip))
    sum_net_heat_glc    = sum(data(f_heat_beg:f_heat_end, c_glc_grecv, ip)) + &
                          sum(data(f_heat_beg:f_heat_end, c_glc_gsend, ip))
    sum_net_heat_tot    = sum_net_heat_atm + sum_net_heat_lnd + sum_net_heat_rof + sum_net_heat_ocn + &
                          sum_net_heat_ice_nh + sum_net_heat_ice_sh + sum_net_heat_glc

    write(logunit,FA1r)'   *SUM*',&
         sum_net_heat_atm, sum_net_heat_lnd, sum_net_heat_rof, sum_net_heat_ocn, &
         sum_net_heat_ice_nh, sum_net_heat_ice_sh, sum_net_heat_glc, sum_net_heat_tot

    ! write out net water budgets

    write(logunit,*) ' '
    write(logunit,FAH) subname,'NET WATER BUDGET (kg/m2s*1e6): period = ',trim(period_name(ip)),': date = ',cdate,curr_tod
    write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
    do nf = f_watr_beg, f_watr_end
       net_water_atm    = data(nf, c_atm_arecv, ip) + data(nf, c_atm_asend, ip)
       net_water_lnd    = data(nf, c_lnd_lrecv, ip) + data(nf, c_lnd_lsend, ip)
       net_water_rof    = data(nf, c_rof_rrecv, ip) + data(nf, c_rof_rsend, ip)
       net_water_ocn    = data(nf, c_ocn_orecv, ip) + data(nf, c_ocn_osend, ip)
       net_water_ice_nh = data(nf, c_inh_irecv, ip) + data(nf, c_inh_isend, ip)
       net_water_ice_sh = data(nf, c_ish_irecv, ip) + data(nf, c_ish_isend, ip)
       net_water_glc    = data(nf, c_glc_grecv, ip) + data(nf, c_glc_gsend, ip)
       net_water_tot    = net_water_atm + net_water_lnd + net_water_rof + net_water_ocn + &
                          net_water_ice_nh + net_water_ice_sh + net_water_glc

       write(logunit,FA1r)   fname(nf),&
            net_water_atm, net_water_lnd, net_water_rof, net_water_ocn, &
            net_water_ice_nh, net_water_ice_sh, net_water_glc, net_water_tot
    enddo

    ! Write out sum over all net heat budgets (sum over f_watr_beg -> f_watr_end)

    sum_net_water_atm    = sum(data(f_watr_beg:f_watr_end, c_atm_arecv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_atm_asend, ip))
    sum_net_water_lnd    = sum(data(f_watr_beg:f_watr_end, c_lnd_lrecv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_lnd_lsend, ip))
    sum_net_water_rof    = sum(data(f_watr_beg:f_watr_end, c_rof_rrecv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_rof_rsend, ip))
    sum_net_water_ocn    = sum(data(f_watr_beg:f_watr_end, c_ocn_orecv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_ocn_osend, ip))
    sum_net_water_ice_nh = sum(data(f_watr_beg:f_watr_end, c_inh_irecv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_inh_isend, ip))
    sum_net_water_ice_sh = sum(data(f_watr_beg:f_watr_end, c_ish_irecv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_ish_isend, ip))
    sum_net_water_glc    = sum(data(f_watr_beg:f_watr_end, c_glc_grecv, ip)) + &
                           sum(data(f_watr_beg:f_watr_end, c_glc_gsend, ip))
    sum_net_water_tot    = sum_net_water_atm + sum_net_water_lnd + sum_net_water_rof + sum_net_water_ocn + &
                           sum_net_water_ice_nh + sum_net_water_ice_sh + sum_net_water_glc

    write(logunit,FA1r)'   *SUM*',&
         sum_net_water_atm, sum_net_water_lnd, sum_net_water_rof, sum_net_water_ocn, &
         sum_net_water_ice_nh, sum_net_water_ice_sh, sum_net_water_glc, sum_net_water_tot

    ! write out net water water-isoptope budgets

    if ( flds_wiso ) then

       do is = 1, nisotopes
          write(logunit,*) ' '
          write(logunit,FAH) subname,'NET '//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ', &
               trim(period_name(ip)),': date = ',cdate,curr_tod
          write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
          do nf = iso0(is), isof(is)
             net_water_atm    = data(nf, c_atm_arecv, ip) + data(nf, c_atm_asend, ip)
             net_water_lnd    = data(nf, c_lnd_lrecv, ip) + data(nf, c_lnd_lsend, ip)
             net_water_rof    = data(nf, c_rof_rrecv, ip) + data(nf, c_rof_rsend, ip)
             net_water_ocn    = data(nf, c_ocn_orecv, ip) + data(nf, c_ocn_osend, ip)
             net_water_ice_nh = data(nf, c_inh_irecv, ip) + data(nf, c_inh_isend, ip)
             net_water_ice_sh = data(nf, c_ish_irecv, ip) + data(nf, c_ish_isend, ip)
             net_water_glc    = data(nf, c_glc_grecv, ip) + data(nf, c_glc_gsend, ip)
             net_water_tot    = net_water_atm + net_water_lnd + net_water_rof + net_water_ocn + &
                                net_water_ice_nh + net_water_ice_sh + net_water_glc

             write(logunit,FA1r)   fname(nf),&
                  net_water_atm, net_water_lnd, net_water_rof, net_water_ocn, &
                  net_water_ice_nh, net_water_ice_sh, net_water_glc, net_water_tot
          enddo

          sum_net_water_atm    = sum(data(iso0(is):isof(is), c_atm_arecv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_atm_asend, ip))
          sum_net_water_lnd    = sum(data(iso0(is):isof(is), c_lnd_lrecv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_lnd_lsend, ip))
          sum_net_water_rof    = sum(data(iso0(is):isof(is), c_rof_rrecv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_rof_rsend, ip))
          sum_net_water_ocn    = sum(data(iso0(is):isof(is), c_ocn_orecv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_ocn_osend, ip))
          sum_net_water_ice_nh = sum(data(iso0(is):isof(is), c_inh_irecv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_inh_isend, ip))
          sum_net_water_ice_sh = sum(data(iso0(is):isof(is), c_ish_irecv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_ish_isend, ip))
          sum_net_water_glc    = sum(data(iso0(is):isof(is), c_glc_grecv, ip)) + &
                                 sum(data(iso0(is):isof(is), c_glc_gsend, ip))
          sum_net_water_tot    = sum_nnet_water_atm + sum_nnet_water_lnd + sum_nnet_water_rof + &
                                 sum_nnet_water_ocn + sum_nnet_water_ice_nh + sum_nnet_water_ice_sh + &
                                 sum_nnet_water_glc

          write(logunit,FA1r)'   *SUM*',&
               sum_net_water_atm, sum_net_water_lnd, sum_net_water_rof, sum_net_water_ocn, &
               sum_net_water_ice_nh, sum_net_water_ice_sh, sum_net_water_glc, sum_net_water_tot
       end do
    end if

  end subroutine med_diag_print_net_summary


  !===============================================================================

  subroutine add_to_budget_type(entries ,num, index, name, rc)

    ! input/output variables
    type(budget_type) , intent(inout) :: entries(:)
    integer           , intent(inout) :: num
    integer           , intent(out)   :: index
    character(len=*)  , intent(in)    :: name
    integer :: rc

    ! local variables
    character(len=*), parameter :: subname='(add_to_budget_type)'
    !----------------------------------------------------------------------

    num = num + 1
    if (num > budget_nmax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > budget_nmax "//&
            trim(name), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
      return
    endif
    entries(num)%name = trim(name)
    entries(num)%index = num

  end subroutine add_to_budget_type

end module med_diag_mod
