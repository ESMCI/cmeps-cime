module med_diag

  !----------------------------------------------------------------------------
  ! Compute spatial \& time averages of fluxed quatities
  !
  ! CESM sign convention for fluxes is positive downward with hierarchy being
  !    atm/glc/lnd/rof/ice/ocn
  ! Sign convention:
  !    positive value <=> the model is gaining water, heat, momentum, etc.
  ! Unit convention:
  !    heat flux     ~ W/m^2
  !    momentum flux ~ N/m^2
  !    water flux    ~ (kg/s)/m^2
  !    salt  flux    ~ (kg/s)/m^2
  !----------------------------------------------------------------------------

  use ESMF                  , only : ESMF_Clock
  use shr_kind_mod          , only : r8 => shr_kind_r8, in=>shr_kind_in, i8 => shr_kind_i8,  cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod           , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod         , only : shr_const_rearth, shr_const_pi, shr_const_latice
  use shr_const_mod         , only : shr_const_ice_ref_sal, shr_const_ocn_ref_sal, shr_const_isspval
  use shr_reprosum_mod      , only : shr_reprosum_calc
  use med_internal_state    , only : InternalState, logunit, mastertask
  use shr_nuopc_methods_mod , only : fldchk => shr_nuopc_methods_FB_FldChk

  implicit none
  private

  !PUBLIC MEMBER FUNCTIONS:

  public :: med_diag_zero
  public :: med_diag_atm
  public :: med_diag_lnd
  public :: med_diag_rof
  public :: med_diag_glc
  public :: med_diag_ocn
  public :: med_diag_ice
  public :: med_diag_accum
  public :: med_diag_sum0
  public :: med_diag_print

  !----------------------------------------------------------------------------
  ! Local data
  !----------------------------------------------------------------------------

  !----- local constants -----
  real(r8),parameter :: HFLXtoWFLX = & ! water flux implied by latent heat of fusion
       &  - (shr_const_ocn_ref_sal-shr_const_ice_ref_sal) / &
       &    (shr_const_ocn_ref_sal*shr_const_latice)

  real(r8),parameter :: SFLXtoWFLX = & ! water flux implied by salt flux (kg/m^2s)
       -1._r8/(shr_const_ocn_ref_sal*1.e-3_r8)

  ! ---------------------------------
  ! C for component
  ! ---------------------------------

  ! "r" is receive in the mediator
  ! "s" is send from the mediator

  integer    ,parameter :: c_atm_asend   = 1  ! model index: atm
  integer    ,parameter :: c_atm_arecv   = 2  ! model index: atm
  integer    ,parameter :: c_inh_isend   = 3  ! model index: ice, northern
  integer    ,parameter :: c_inh_irecv   = 4  ! model index: ice, northern
  integer    ,parameter :: c_ish_isend   = 5  ! model index: ice, southern
  integer    ,parameter :: c_ish_irecv   = 6  ! model index: ice, southern
  integer    ,parameter :: c_lnd_lsend   = 7  ! model index: lnd
  integer    ,parameter :: c_lnd_lrecv   = 8  ! model index: lnd
  integer    ,parameter :: c_ocn_osend   = 9  ! model index: ocn
  integer    ,parameter :: c_ocn_orecv   = 10 ! model index: ocn
  integer    ,parameter :: c_rof_rsend   = 11 ! model index: rof
  integer    ,parameter :: c_rof_rrecv   = 12 ! model index: rof
  integer    ,parameter :: c_glc_gsend   = 13 ! model index: glc
  integer    ,parameter :: c_glc_grecv   = 14 ! model index: glc

  ! The folowing is needed for detailing the atm budgets and breakdown into components
  integer    ,parameter :: c_inh_asend   = 15 ! model index: ice, northern, on atm grid
  integer    ,parameter :: c_inh_arecv   = 16 ! model index: ice, northern, on atm grid
  integer    ,parameter :: c_ish_asend   = 17 ! model index: ice, southern, on atm grid
  integer    ,parameter :: c_ish_arecv   = 18 ! model index: ice, southern, on atm grid
  integer    ,parameter :: c_lnd_asend   = 19 ! model index: lnd, on atm grid
  integer    ,parameter :: c_lnd_arecv   = 20 ! model index: lnd, on atm grid
  integer    ,parameter :: c_ocn_asend   = 21 ! model index: ocn, on atm grid
  integer    ,parameter :: c_ocn_arecv   = 22 ! model index: ocn, on atm grid

  ! Total size
  integer    ,parameter :: c_size     = 22

  ! Note - the character strings here MUST match the parameters above
  character(len=8),parameter :: cname(c_size) = &
       (/' c2a_atm', & ! 1
         ' a2c_atm', & ! 2
         ' c2i_inh', & ! 3
         ' i2c_inh', & ! 4
         ' c2i_ish', & ! 5
         ' i2c_ish', & ! 6
         ' c2l_lnd', & ! 7
         ' l2c_lnd', & ! 8
         ' c2o_ocn', & ! 9
         ' o2c_ocn', & ! 10
         ' c2r_rof', & ! 11
         ' r2c_rof', & ! 12
         ' c2g_glc', & ! 13
         ' g2c_glc', & ! 14
         ' c2a_inh', & ! 15
         ' a2c_inh', & ! 16
         ' c2a_ish', & ! 17
         ' a2c_ish', & ! 18
         ' c2a_lnd', & ! 19
         ' a2c_lnd', & ! 20
         ' c2a_ocn', & ! 21
         ' a2c_ocn'  & ! 22
         /)

  ! ---------------------------------
  ! F for field
  ! ---------------------------------

  integer , parameter :: f_area      =  1         ! area (wrt to unit sphere)
  integer , parameter :: f_hfrz      =  2         ! heat : latent, freezing
  integer , parameter :: f_hmelt     =  3         ! heat : latent, melting
  integer , parameter :: f_hswnet    =  4         ! heat : short wave, net
  integer , parameter :: f_hlwdn     =  5         ! heat : longwave down
  integer , parameter :: f_hlwup     =  6         ! heat : longwave up
  integer , parameter :: f_hlatv     =  7         ! heat : latent, vaporization
  integer , parameter :: f_hlatf     =  8         ! heat : latent, fusion, snow
  integer , parameter :: f_hioff     =  9         ! heat : latent, fusion, frozen runoff
  integer , parameter :: f_hsen      = 10         ! heat : sensible
  integer , parameter :: f_wfrz      = 11         ! water: freezing
  integer , parameter :: f_wmelt     = 12         ! water: melting
  integer , parameter :: f_wrain     = 13         ! water: precip, liquid
  integer , parameter :: f_wsnow     = 14         ! water: precip, frozen
  integer , parameter :: f_wevap     = 15         ! water: evaporation
  integer , parameter :: f_wsalt     = 16         ! water: water equivalent of salt flux
  integer , parameter :: f_wroff     = 17         ! water: runoff/flood
  integer , parameter :: f_wioff     = 18         ! water: frozen runoff
  integer , parameter :: f_wfrz_16O  = 19         ! water: freezing
  integer , parameter :: f_wmelt_16O = 20         ! water: melting
  integer , parameter :: f_wrain_16O = 21         ! water: precip, liquid
  integer , parameter :: f_wsnow_16O = 22         ! water: precip, frozen
  integer , parameter :: f_wevap_16O = 23         ! water: evaporation
  integer , parameter :: f_wroff_16O = 24         ! water: runoff/flood
  integer , parameter :: f_wioff_16O = 25         ! water: frozen runoff
  integer , parameter :: f_wfrz_18O  = 26         ! water: freezing
  integer , parameter :: f_wmelt_18O = 27         ! water: melting
  integer , parameter :: f_wrain_18O = 28         ! water: precip, liquid
  integer , parameter :: f_wsnow_18O = 29         ! water: precip, frozen
  integer , parameter :: f_wevap_18O = 30         ! water: evaporation
  integer , parameter :: f_wroff_18O = 31         ! water: runoff/flood
  integer , parameter :: f_wioff_18O = 32         ! water: frozen runoff
  integer , parameter :: f_wfrz_HDO  = 33         ! water: freezing
  integer , parameter :: f_wmelt_HDO = 34         ! water: melting
  integer , parameter :: f_wrain_HDO = 35         ! water: precip, liquid
  integer , parameter :: f_wsnow_HDO = 36         ! water: precip, frozen
  integer , parameter :: f_wevap_HDO = 37         ! water: evaporation
  integer , parameter :: f_wroff_HDO = 38         ! water: runoff/flood
  integer , parameter :: f_wioff_HDO = 39         ! water: frozen runoff
  integer , parameter :: f_size     = f_wioff_HDO ! Total array size of all elements

  integer , parameter :: f_a        = f_area      ! 1st index for area
  integer , parameter :: f_a_end    = f_area      ! last index for area
  integer , parameter :: f_h        = f_hfrz      ! 1st index for heat
  integer , parameter :: f_h_end    = f_hsen      ! Last index for heat
  integer , parameter :: f_w        = f_wfrz      ! 1st index for water
  integer , parameter :: f_w_end    = f_wioff     ! Last index for water
  integer , parameter :: f_16O      = f_wfrz_16O  ! 1st index for 16O water isotope
  integer , parameter :: f_18O      = f_wfrz_18O  ! 1st index for 18O water isotope
  integer , parameter :: f_HDO      = f_wfrz_HDO  ! 1st index for HDO water isotope
  integer , parameter :: f_16O_end  = f_wioff_16O ! Last index for 16O water isotope
  integer , parameter :: f_18O_end  = f_wioff_18O ! Last index for 18O water isotope
  integer , parameter :: f_HDO_end  = f_wioff_HDO ! Last index for HDO water isotope

  ! Note - the character strings here MUST match the parameters values above
  character(len=12), parameter :: fname(f_size) = &
       (/'      area', & ! 1
       '     hfreeze', & ! 2
       '       hmelt', & ! 3
       '      hnetsw', & ! 4
       '       hlwdn', & ! 5
       '       hlwup', & ! 6
       '     hlatvap', & ! 7
       '     hlatfus', & ! 8
       '      hiroff', & ! 9
       '        hsen', & ! 10
       '     wfreeze', & ! 11
       '       wmelt', & ! 12
       '       wrain', & ! 13
       '       wsnow', & ! 14
       '       wevap', & ! 15
       '    weqsaltf', & ! 16
       '     wrunoff', & ! 17
       '     wfrzrof', & ! 18
       ' wfreeze_16O', & ! 19
       '   wmelt_16O', & ! 20
       '   wrain_16O', & ! 21
       '   wsnow_16O', & ! 22
       '   wevap_16O', & ! 23
       ' wrunoff_16O', & ! 24
       ' wfrzrof_16O', & ! 25
       ' wfreeze_18O', & ! 26
       '   wmelt_18O', & ! 27
       '   wrain_18O', & ! 28
       '   wsnow_18O', & ! 29
       '   wevap_18O', & ! 30
       ' wrunoff_18O', & ! 31
       ' wfrzrof_18O', & ! 32
       ' wfreeze_HDO', & ! 33
       '   wmelt_HDO', & ! 34
       '   wrain_HDO', & ! 35
       '   wsnow_HDO', & ! 36
       '   wevap_HDO', & ! 37
       ' wrunoff_HDO', & ! 38
       ' wfrzrof_HDO'/)  ! 39

  ! ---------------------------------
  ! P for period
  ! ---------------------------------

  integer    ,parameter :: period_ntypes = 5

  integer    ,parameter :: period_inst = 1
  integer    ,parameter :: period_day  = 2
  integer    ,parameter :: period_mon  = 3
  integer    ,parameter :: period_ann  = 4
  integer    ,parameter :: period_inf  = 5

  character(len=8),parameter :: pname(period_ntypes) = &
       (/'    inst','   daily',' monthly','  annual','all_time' /)

  logical :: flds_wiso ! If water isotope fields are active

  ! ---------------------------------
  ! public data members
  ! ---------------------------------

  !--- time-averaged (annual?) global budge diagnostics ---
  !--- note: call sum0 then save budget_global and budget_ns on restart from/to root pe ---

  real(r8),public :: budget_local (f_size,c_size,period_ntypes) ! local sum, valid on all pes
  real(r8),public :: budget_global(f_size,c_size,period_ntypes) ! global sum, valid only on root pe
  real(r8),public :: budget_ns           (f_size,c_size,period_ntypes) ! counter, valid only on root pe

  character(len=*),parameter :: afldname  = 'aream'
  character(len=*),parameter :: latname   = 'lat'
  character(len=*),parameter :: afracname = 'afrac'
  character(len=*),parameter :: lfracname = 'lfrac'
  character(len=*),parameter :: ofracname = 'oofrac'
  character(len=*),parameter :: ifracname = 'ifrac'
  character(len=*),parameter :: modName   = "(med_diag) "
  integer         ,parameter :: debug     = 0 ! internal debug level

!===============================================================================
contains
!===============================================================================

  subroutine med_diag_zero( EClock, mode)

    ! !DESCRIPTION: Zero out global budget diagnostic data.

    ! input/output variables
    type(ESMF_Clock), intent(in),optional :: EClock
    character(len=*), intent(in),optional :: mode

    ! local variables
    integer     :: ip,yr,mon,day,sec
    character(*),parameter :: subName = '(med_diag_zero) '
    !-------------------------------------------------------------------------------

    if (.not. present(EClock) .and. .not. present(mode)) then
       call shr_sys_abort(subName//' ERROR EClock or mode should be present')
    endif

    if (present(EClock)) then
       call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet( currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       do ip = 1,period_ntypes
          if (ip == period_inst) then
             budget_local(:,:,ip) = 0.0_r8
             budget_global(:,:,ip) = 0.0_r8
             budget_ns(:,:,ip) = 0.0_r8
          endif
          if (ip==period_day .and. current_tod==0) then
             budget_local(:,:,ip) = 0.0_r8
             budget_global(:,:,ip) = 0.0_r8
             budget_ns(:,:,ip) = 0.0_r8
          endif
          if (ip==period_mon .and. current_day==1 .and. current_tod==0) then
             budget_local(:,:,ip) = 0.0_r8
             budget_global(:,:,ip) = 0.0_r8
             budget_ns(:,:,ip) = 0.0_r8
          endif
          if (ip==period_ann .and. current_mon==1 .and. current_day==1 .and. current_tod==0) then
             budget_local(:,:,ip) = 0.0_r8
             budget_global(:,:,ip) = 0.0_r8
             budget_ns(:,:,ip) = 0.0_r8
          endif
       enddo

    else if(present(mode)) then

       if (trim(mode) == 'inst') then
          budget_local(:,:,period_inst) = 0.0_r8
          budget_global(:,:,period_inst) = 0.0_r8
          budget_ns(:,:,period_inst) = 0.0_r8
       elseif (trim(mode) == 'day') then
          budget_local(:,:,period_day) = 0.0_r8
          budget_global(:,:,period_day) = 0.0_r8
          budget_ns(:,:,period_day) = 0.0_r8
       elseif (trim(mode) == 'mon') then
          budget_local(:,:,period_mon) = 0.0_r8
          budget_global(:,:,period_mon) = 0.0_r8
          budget_ns(:,:,period_mon) = 0.0_r8
       elseif (trim(mode) == 'ann') then
          budget_local(:,:,period_ann) = 0.0_r8
          budget_global(:,:,period_ann) = 0.0_r8
          budget_ns(:,:,period_ann) = 0.0_r8
       elseif (trim(mode) == 'inf') then
          budget_local(:,:,period_inf) = 0.0_r8
          budget_global(:,:,period_inf) = 0.0_r8
          budget_ns(:,:,period_inf) = 0.0_r8
       elseif (trim(mode) == 'all') then
          budget_local(:,:,:) = 0.0_r8
          budget_global(:,:,:) = 0.0_r8
          budget_ns(:,:,:) = 0.0_r8
       else
          call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
       endif
    endif

  end subroutine med_diag_zero

  !===============================================================================

  subroutine med_diag_accum(rc)

    ! ------------------------------------------------------------------
    ! Accumulate out global budget diagnostic data.
    ! ------------------------------------------------------------------

    ! input/output variables
    integer, intent(out) :: rc

    ! local variables
    integer     :: ip
    character(*),parameter :: subName = '(med_diag_accum) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    do ip = period_inst+1,period_ntypes
       budget_local(:,:,ip) = budget_local(:,:,ip) + budget_local(:,:,period_inst)
    enddo
    budget_ns(:,:,:) = budget_ns(:,:,:) + 1.0_r8

  end subroutine med_diag_accum

  !===============================================================================

  subroutine med_diag_sum0(gcomp, rc)

    ! ------------------------------------------------------------------
    ! Sum local values to global on root
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)       :: vm
    integer             :: count
    real(r8)            :: budget_global_tmp(f_size,c_size,period_ntypes) ! temporary sum
    character(*),parameter :: subName = '(med_diag_sum0) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp, vm, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    count = size(budget_global_tmp)
    budget_global_tmp(:,:,:) = 0.0_r8
    call ESMF_VMReduce(vm, budget_local, budget_global_tmp, count, &
         ESMF_REDUCE_SUM, 0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    budget_global(:,:,:) = budget_global(:,:,;) + budget_global_tmp(:,:,:)
    budget_local(:,:,:) = 0.0_r8

  end subroutine med_diag_sum0

  !===============================================================================

  subroutine med_diag_atm(gcomp,  rc)

    ! ------------------------------------------------------------------
    ! Compute global atm input/output flux diagnostics
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(in)           :: gcomp
    integer             , intent(out)          :: rc

    ! local variables
    type(InternalState)    :: is_local
    integer                :: n,ic,ip
    real(r8)               :: wgt, area
    logical                :: atm_wiso_recv = .false.
    logical                :: atm_wiso_send = .false.
    real(r8), pointer      :: afrac(:)
    real(r8), pointer      :: lfrac(:)
    real(r8), pointer      :: ifrac(:)
    real(r8), pointer      :: ofrac(:)
    real(r8), pointer      :: Faxa_swnet(:)        ! from atm to mediator
    real(r8), pointer      :: Faxa_lwdn(:)         ! from atm to mediator
    real(r8), pointer      :: Faxa_rainc(:)        ! from atm to mediator
    real(r8), pointer      :: Faxa_rainl(:)        ! from atm to mediator
    real(r8), pointer      :: Faxa_snowc(:)        ! from atm to mediator
    real(r8), pointer      :: Faxa_snowl(:)        ! from atm to mediator
    real(r8), pointer      :: Faxa_rainl_wiso(:,:) ! from atm to mediator
    real(r8), pointer      :: Faxa_snowl_wiso(:,:) ! from atm to mediator
    real(r8), pointer      :: Faxx_lwup(:)         ! from  mediator to atm
    real(r8), pointer      :: Faxx_lat(:)          ! from  mediator to atm
    real(r8), pointer      :: Faxx_sen(:)          ! from  mediator to atm
    real(r8), pointer      :: Faxx_evap(:)         ! from  mediator to atm
    real(r8), pointer      :: Faxx_evap_wiso(:,:)  ! from  mediator to atm
    character(*),parameter :: subName = '(med_diag_atm) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get fractions on atm mesh
    call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'afrac', afrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', lfrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ifrac', ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------
    ! from atm to mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_swnet', Faxa_swnet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lwdn' , Faxa_lwdn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainc', Faxa_rainc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl', Faxa_rainl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowc', Faxa_snowc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl', Faxa_snowl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(FBlnd, 'Faxa_rainl_wiso', rc=rc) .and. &
         fldchk(FBlnd, 'Faxa_snowl_wiso', rc=rc)) then
       atm_wiso_recv = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_rainl_wiso', Faxa_rainl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faxa_snowl_wiso', Faxa_snowl_wiso, rc=rc)
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

          budget_local(f_area  ,ic,ip) = budget_local(f_area  ,ic,ip) + wgt
          budget_local(f_hswnet,ic,ip) = budget_local(f_hswnet,ic,ip) + wgt*Faxa_swnet(n)
          budget_local(f_hlwdn ,ic,ip) = budget_local(f_hlwdn ,ic,ip) + wgt*Faxa_lwdn(n)
          budget_local(f_wrain ,ic,ip) = budget_local(f_wrain ,ic,ip) + wgt*(Faxa_rainc(n) Faxa_rainl(n))
          budget_local(f_wsnow ,ic,ip) = budget_local(f_wsnow ,ic,ip) + wgt*(Faxa_snowc(n) Faxa_snowl(n))

          if ( atm_wiso_recv )then
             budget_local(f_wrain_16O,ic,ip) = budget_local(f_wrain_16O,ic,ip) + wgt*(Faxa_rainc_wiso(1,n) + Faxa_rainl_wiso(1,n))
             budget_local(f_wrain_18O,ic,ip) = budget_local(f_wrain_18O,ic,ip) + wgt*(Faxa_rainc_wiso(2,n) + Faxa_rainl_wiso(2,n))
             budget_local(f_wrain_HDO,ic,ip) = budget_local(f_wrain_HDO,ic,ip) + wgt*(Faxa_rainc_wiso(3,n) + Faxa_rainl_wiso(3,n))
          end if
       enddo
    enddo

    ! heat implied by snow flux
    budget_local(f_hlatf,c_atm_arecv,ip) = -budget_local(f_wsnow,c_atm_arecv,ip)*shr_const_latice
    budget_local(f_hlatf,c_lnd_arecv,ip) = -budget_local(f_wsnow,c_lnd_arecv,ip)*shr_const_latice
    budget_local(f_hlatf,c_ocn_arecv,ip) = -budget_local(f_wsnow,c_ocn_arecv,ip)*shr_const_latice
    budget_local(f_hlatf,c_inh_arecv,ip) = -budget_local(f_wsnow,c_inh_arecv,ip)*shr_const_latice
    budget_local(f_hlatf,c_ish_arecv,ip) = -budget_local(f_wsnow,c_ish_arecv,ip)*shr_const_latice

    !-------------------------------
    ! from mediator to atm
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_lwup', Faxx_lwup, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_lat' , Faxx_lat , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_sen' , Faxx_sen , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_evap', Faxx_evap, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBExp(compatm), 'Faxx_evap_wiso', Faxa_evap_wiso, rc=rc)) then
       atm_wiso_send = .true.
       call FB_GetFldPtr(is_local%wrap%FBExp(compatm), 'Faxx_evap_wiso', Faxx_evap_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    do n = 1,size(Faxx_lwup)

       ! TODO: this loop must be over the number of components that are being merged to obtain the
       ! atm input plus 1
       do k=1,4  ! sum over merged, lfrac, ifrac, ofrac (3 inputs to the atm, 1 merged input)

          ! Determine weight
          if (k == 1) then
             wgt = -areas(n) * afrac(n)
          else if (k == 2) then
             wgt =  areas(n) * lfrac(n)
          elseif (k == 3) then
             wgt =  areas(n) * ifrac(n)
          elseif (k == 4) then
             wgt =  areas(n) * ofrac(n)
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

          budget_local(f_area ,ic,ip) = budget_local(f_area ,ic,ip) + wgt
          budget_local(f_hlwup,ic,ip) = budget_local(f_hlwup,ic,ip) + wgt*Faxx_lwup(n)
          budget_local(f_hlatv,ic,ip) = budget_local(f_hlatv,ic,ip) + wgt*Faxx_lat(n)
          budget_local(f_hsen ,ic,ip) = budget_local(f_hsen ,ic,ip) + wgt*Faxx_sen(n)
          budget_local(f_wevap,ic,ip) = budget_local(f_wevap,ic,ip) + wgt*Faxx_evap(n)

          if ( atm_wiso_send )then
             budget_local(f_wevap_16O,ic,ip) = budget_local(f_wevap_16O,ic,ip) + wgt*Faxx_evap_wiso(1,n)
             budget_local(f_wevap_18O,ic,ip) = budget_local(f_wevap_18O,ic,ip) + wgt*Faxx_evap_wiso(2,n)
             budget_local(f_wevap_HDO,ic,ip) = budget_local(f_wevap_HDO,ic,ip) + wgt*Faxx_evap_wiso(3,n)
          end if

       enddo
    enddo

  end subroutine med_diag_atm

  !===============================================================================

  subroutine med_diag_lnd( gcomp, rc)

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
    real(r8), pointer      :: Flrl_roftdo(:)       ! to mediator from land
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
    real(r8), pointer      :: Flrr_flood_wiso(:,:) ! from mediator to land
    character(*),parameter :: subName = '(med_diag_lnd) '
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

    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_swnet', Fall_swnet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_lwup', Fall_lwup, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_lat', Fall_lat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_sen', Fall_sen, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap', Fall_evap, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofsur', Flrl_rofsur, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofgwl', Flrl_rofgwl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofsub', Flrl_rofsub, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_roftdo', Flrl_roftdo, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofi', Flrl_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fldchk(is_local%wrap%FBImp(complnd,complnd), 'Flrl_irrig', rc=rc)) then
       lnd_irrig = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_irrig', Flrl_irrig, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       lnd_irrig = .false.
    end if

    if ( fldchk(FBlnd, 'Fall_evap_wiso', rc=rc) .and. &
         fldchk(FBlnd, 'Flrl_rofl_wiso', rc=rc) .and. &
         fldchk(FBlnd, 'Flrl_rofi_wiso', rc=rc)) then
       lnd_wiso_recv = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Fall_evap_wiso', Fall_evap_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofl_wiso', Flrl_rofl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(complnd,complnd), 'Flrl_rofi_wiso', Flrl_rofi_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    ic = c_lnd_lrecv
    do n = 1, size(Fall_swnet)
       wgt = is_local%wrap%mesh_info(complnd)%areas(n) * lfrac(n)

       budget_local(f_area  ,ic,ip) = budget_local(f_area  ,ic,ip) + wgt
       budget_local(f_hswnet,ic,ip) = budget_local(f_hswnet,ic,ip) + wgt*Fall_swnet(n)
       budget_local(f_hlwup ,ic,ip) = budget_local(f_hlwup ,ic,ip) + wgt*Fall_lwup(n)
       budget_local(f_hlatv ,ic,ip) = budget_local(f_hlatv ,ic,ip) + wgt*Fall_lat(n)
       budget_local(f_hsen  ,ic,ip) = budget_local(f_hsen  ,ic,ip) + wgt*Fall_sen(n)
       budget_local(f_wevap ,ic,ip) = budget_local(f_wevap ,ic,ip) + wgt*Fall_evap(n)

       budget_local(f_wroff ,ic,ip) = budget_local(f_wroff ,ic,ip) &
            - wgt*Flrl_rofsur(n) - wgt*Flrl_rofgwl(n) - wgt*Flrl_rofsub(n) - wgt*Flrl_rofdto(n)
       if (lnd_irrig) then
          budget_local(f_wroff,ic,ip) = budget_local(f_wroff,ic,ip) - wgt*Flrl_irrig(n)
       end if
       budget_local(f_wioff ,ic,ip) = budget_local(f_wioff,ic,ip)  - wgt*Flrl_rofi(n)

       if ( lnd_wiso_recv ) then
          budget_local(f_wevap_16O,ic,ip) = budget_local(f_wevap_16O,ic,ip) + wgt*Fall_evap_wiso(1,n)
          budget_local(f_wevap_18O,ic,ip) = budget_local(f_wevap_18O,ic,ip) + wgt*Fall_evap_wiso(2,n)
          budget_local(f_wevap_HDO,ic,ip) = budget_local(f_wevap_HDO,ic,ip) + wgt*Fall_evap_wiso(3,n)

          budget_local(f_wroff_16O,ic,ip) = budget_local(f_wroff_16O,ic,ip) - wgt*Flrl_rofl_wiso(1,n)
          budget_local(f_wroff_18O,ic,ip) = budget_local(f_wroff_18O,ic,ip) - wgt*Flrl_rofl_wiso(2,n)
          budget_local(f_wroff_HDO,ic,ip) = budget_local(f_wroff_HDO,ic,ip) - wgt*Flrl_rofl_wiso(3,n)

          budget_local(f_wioff_16O,ic,ip) = budget_local(f_wioff_16O,ic,ip) - wgt*Flrl_rofi_wiso(1,n)
          budget_local(f_wioff_18O,ic,ip) = budget_local(f_wioff_18O,ic,ip) - wgt*Flrl_rofi_wiso(2,n)
          budget_local(f_wioff_HDO,ic,ip) = budget_local(f_wioff_HDO,ic,ip) - wgt*Flrl_rofi_wiso(3,n)
       end if
    end do
    budget_local(f_hioff,ic,ip) = -budget_local(f_wioff,ic,ip)*shr_const_latice

    !-------------------------------
    ! to land from mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_lwdn', Faxa_lwdn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainc', Faxa_rainc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainl', Faxa_rainl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainc', Faxa_snowc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainl', Faxa_snowl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Flrl_flood', Flrl_flood, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(FBlnd, 'Faxa_rainl_wiso', rc=rc) .and. &
         fldchk(FBlnd, 'Faxa_snowl_wiso', rc=rc) .and. &
         fldchk(FBlnd, 'Flrl_flood_wiso', rc=rc)) then
       lnd_wiso_send = .true.
       call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainl_wiso', Faxa_rainl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_rainc_wiso', Faxa_rainc_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_snowl_wiso', Faxa_snowl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Faxa_snowc_wiso', Faxa_snowc_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(complnd), 'Flrl_flood_wiso', Flrl_flood_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    ic = c_lnd_lsend
    do n = 1,size(Faxa_lwdn)
       wgt = is_local%wrap%mesh_info(complnd)%areas(n) * lfrac(n)

       budget_local(f_area ,ic,ip) = budget_local(f_area ,ic,ip) + wgt
       budget_local(f_hlwdn,ic,ip) = budget_local(f_hlwdn,ic,ip) + wgt*Faxa_lwdn(n)
       budget_local(f_wrain,ic,ip) = budget_local(f_wrain,ic,ip) + wgt*Faxa_rainc(n) + wgt*Faxa_rainl(n)
       budget_local(f_wsnow,ic,ip) = budget_local(f_wsnow,ic,ip) + wgt*Faxa_snowc(n) + wgt*Faxa_snowl(n)
       budget_local(f_wroff,ic,ip) = budget_local(f_wroff,ic,ip) - wgt*Flrr_flood(n)

       if ( lnd_wiso_send )then
          budget_local(f_wrain_16O,ic,ip) = budget_local(f_wrain_16O,ic,ip) + wgt*Faxa_rainc_wiso(1,n) + wgt*Faxa_rainl_wiso(1,n)
          budget_local(f_wrain_18O,ic,ip) = budget_local(f_wrain_18O,ic,ip) + wgt*Faxa_rainc_wiso(2,n) + wgt*Faxa_rainl_wiso(2,n)
          budget_local(f_wrain_HDO,ic,ip) = budget_local(f_wrain_HDO,ic,ip) + wgt*Faxa_rainc_wiso(3,n) + wgt*Faxa_rainl_wiso(3,n)

          budget_local(f_wsnow_16O,ic,ip) = budget_local(f_wsnow_16O,ic,ip) + wgt*Faxa_snowc_wiso(1,n) + wgt*Faxa_snowl_wiso(1,n)
          budget_local(f_wsnow_18O,ic,ip) = budget_local(f_wsnow_18O,ic,ip) + wgt*Faxa_snowc_wiso(2,n) + wgt*Faxa_snowl_wiso(2,n)
          budget_local(f_wsnow_HDO,ic,ip) = budget_local(f_wsnow_HDO,ic,ip) + wgt*Faxa_snowc_wiso(3,n) + wgt*Faxa_snowl_wiso(3,n)

          budget_local(f_wroff_16O,ic,ip) = budget_local(f_wroff_16O,ic,ip) - wgt*Flrr_flood_wiso(1,n)
          budget_local(f_wroff_18O,ic,ip) = budget_local(f_wroff_18O,ic,ip) - wgt*Flrr_flood_wiso(2,n)
          budget_local(f_wroff_HDO,ic,ip) = budget_local(f_wroff_HDO,ic,ip) - wgt*Flrr_flood_wiso(3,n)
       end if

    end do
    budget_local(f_hlatf,ic,ip) = -budget_local(f_wsnow,ic,ip)*shr_const_latice

  end subroutine med_diag_lnd

  !===============================================================================

  subroutine med_diag_rof( gcomp, rc)

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
    integer                :: n,ic,ip
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
    character(*),parameter :: subName = '(med_diag_rof) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------
    ! to river from mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofsur', Flrl_rofsur, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofgwl', Flrl_rofgwl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofsub', Flrl_rofsub, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofdto', Flrl_rofdto, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofi', Flrl_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (fldchk(is_local%wrap%FBImp(comprof,comprof), 'Flrl_irrig', rc=rc)) then
       rof_irrig = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_irrig', Flrl_irrig, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if ( fldchk(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofl_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofi_wiso', rc=rc)) then
       rof_wiso = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofl_wiso', Flrl_rofl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(comprof,comprof), 'Flrl_rofi_wiso', Flrl_rofi_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    ic = c_rof_rrecv
    do n = 1,size(Flrl_rofdto(n)
       wgt = is_local%wrap%mesh_info(comprof)%areas(n)

       budget_local(f_wroff,ic,ip) = budget_local(f_wroff,ic,ip) &
            + wgt*Flrl_rofsur(n) + wgt*Flrl_rofgwl(n) + wgt*Flrl_rofsub(n) + wgt*Flrl_rofdto(n)
       if (rof_irrig /= 0) then
          budget_local(f_wroff,ic,ip) = budget_local(f_wroff,ic,ip) + wgt*Flrl_irrig(n)
       end if
       budget_local(f_wioff,ic,ip) = budget_local(f_wioff,ic,ip) + wgt*Flrl_rofi(n)

       if ( rof_wiso )then
          budget_local(f_wroff_16O,ic,ip) = budget_local(f_wroff_16O,ic,ip) + wgt*Flrl_rofl_wiso(1,n)
          budget_local(f_wroff_18O,ic,ip) = budget_local(f_wroff_18O,ic,ip) + wgt*Flrl_rofl_wiso(2,n)
          budget_local(f_wroff_HDO,ic,ip) = budget_local(f_wroff_HDO,ic,ip) + wgt*Flrl_rofl_wiso(3,n)

          budget_local(f_wioff_16O,ic,ip) = budget_local(f_wioff_16O,ic,ip) + wgt*Flrl_rofi_wiso(1,n)
          budget_local(f_wioff_18O,ic,ip) = budget_local(f_wioff_18O,ic,ip) + wgt*Flrl_rofi_wiso(2,n)
          budget_local(f_wioff_HDO,ic,ip) = budget_local(f_wioff_HDO,ic,ip) + wgt*Flrl_rofi_wiso(3,n)
       end if
    end do
    budget_local(f_hioff,ic,ip) = -budget_local(f_wioff,ic,ip)*shr_const_latice

    !-------------------------------
    ! from river to mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Forr_rofl', Forr_rofl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Flrl_flood', Flrl_flood, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Forr_rofi', Forr_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Firr_rofi', Firr_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBExp(comprof), 'Forr_rofl_wiso' , rc=rc) .and. &
         fldchk(is_local%wrap%FBExp(comprof), 'Forr_rofi_wiso' , rc=rc) .and. &
         fldchk(is_local%wrap%FBExp(comprof), 'Flrl_flood_wiso', rc=rc)) then
       rof_wiso = .true.
       call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Forr_rofl_wiso', Forr_rofl_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Forr_rofi_wiso', Forr_rofi_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBExp(comprof), 'Flrl_flood_wiso', Flrl_flood_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ip = period_inst
    ic = c_rof_rsend
    do n = 1,size(Flrr_flood(n)
       wgt = is_local%wrap%mesh_info(comprof)%areas(n)

       budget_local(f_wroff,ic,ip) = budget_local(f_wroff,ic,ip) - wgt*Forr_rofl(n) + wgt*Flrr_flood(n)
       budget_local(f_wioff,ic,ip) = budget_local(f_wioff,ic,ip) - wgt*Forr_rofi(n) - wgt*Firr_rofi(n)

       if ( rof_wiso )then
          budget_local(f_wroff_16O,ic,ip) = budget_local(f_wroff_16O,ic,ip) - wgt*Forr_rofl_wiso(1,n)
          budget_local(f_wroff_18O,ic,ip) = budget_local(f_wroff_18O,ic,ip) - wgt*Forr_rofl_wiso(2,n)
          budget_local(f_wroff_HDO,ic,ip) = budget_local(f_wroff_HDO,ic,ip) - wgt*Forr_rofl_wiso(3,n)

          budget_local(f_wioff_16O,ic,ip) = budget_local(f_wioff_16O,ic,ip) - wgt*Forr_rofi_wiso(1,n)
          budget_local(f_wioff_18O,ic,ip) = budget_local(f_wioff_18O,ic,ip) - wgt*Forr_rofi_wiso(2,n)
          budget_local(f_wioff_HDO,ic,ip) = budget_local(f_wioff_HDO,ic,ip) - wgt*Forr_rofi_wiso(3,n)

          budget_local(f_wroff_16O,ic,ip) = budget_local(f_wroff_16O,ic,ip) + wgt*Flrr_flood_wiso(1,n)
          budget_local(f_wroff_18O,ic,ip) = budget_local(f_wroff_18O,ic,ip) + wgt*Flrr_flood_wiso(2,n)
          budget_local(f_wroff_HDO,ic,ip) = budget_local(f_wroff_HDO,ic,ip) + wgt*Flrr_flood_wiso(3,n)
       end if
    end do
    budget_local(f_hioff,ic,ip) = -budget_local(f_wioff,ic,ip)*shr_const_latice

  end subroutine med_diag_rof

  !===============================================================================

  subroutine med_diag_glc( gcomp, rc)

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
    character(*),parameter :: subName = '(med_diag_glc) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------
    ! from glc to mediator
    !-------------------------------

    call FB_GetFldPtr(is_local%wrap%FBImp(compglc,compglc), 'Fogg_rofl', Fogg_rofl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compglc,compglc), 'Fogg_rofi', Fogg_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compglc,compglc), 'Figg_rofi', Figg_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ip = period_inst
    ic = c_glc_gsend
    do n = 1,size(Fogg_rofl)
       wgt = is_local%wrap%mesh_info(compglc)%areas(n)
       budget_local(f_wroff,ic,ip) = budget_local(f_wroff,ic,ip) - wgt*Fogg_rofl(n)
       budget_local(f_wioff,ic,ip) = budget_local(f_wioff,ic,ip) - wgt*Fogg_rofi(n) - wgt*Figg_rofi(n)
    end do
    budget_local(f_hioff,ic,ip) = -budget_local(f_wioff,ic,ip)*shr_const_latice

  end subroutine med_diag_glc

  !===============================================================================

  subroutine med_diag_ocn_ocn2med( gcomp, rc)

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
    character(*),parameter :: subName = '(med_diag_ocn) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac', ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac', ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ip = period_inst
    ic = c_ocn_orecv
    do n = 1,size(Fioo_q)
       wgt_o = is_local%wrap%mesh_info(compocn)%areas(n) * ofrac(n)
       wgt_i = is_local%wrap%mesh_info(compocn)%areas(n) * ifrac(n)

       budget_local(f_area,ic,ip) = budget_local(f_area,ic,ip) + wgt_o
       budget_local(f_hfrz,ic,ip) = budget_local(f_hfrz,ic,ip) + (wgt_o + wgt_i)*max(0.0_r8,Fioo_q(n))
    end do
    budget_local(f_wfrz,ic,ip) = budget_local(f_hfrz,ic,ip) * HFLXtoWFLX

    if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lwup', rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_lat' , rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_sen' , rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Faox_evap', rc=rc)) then
       ip = period_inst
       ic = c_ocn_orecv
       do n = 1,size(Faox_lwup)
          wgt_o = is_local%wrap%mesh_info(compocn)%areas(n) * ofrac(n)

          budget_local(f_hlwup,ic,ip) = budget_local(f_hlwup,ic,ip) + wgt_o*Faox_lwup(n)
          budget_local(f_hlatv,ic,ip) = budget_local(f_hlatv,ic,ip) + wgt_o*Faox_lat(n)
          budget_local(f_hsen ,ic,ip) = budget_local(f_hsen ,ic,ip) + wgt_o*Faox_sen(n)
          budget_local(f_wevap,ic,ip) = budget_local(f_wevap,ic,ip) + wgt_o*Faox_evap(n)

          if ( ocn_wiso) then
             budget_local(f_wevap_16O,ic,ip) = budget_local(f_wevap_16O,ic,ip) + wgt_o*Faox_evap_wiso(1,n)
             budget_local(f_wevap_18O,ic,ip) = budget_local(f_wevap_18O,ic,ip) + wgt_o*Faox_evap_wiso(2,n)
             budget_local(f_wevap_HDO,ic,ip) = budget_local(f_wevap_HDO,ic,ip) + wgt_o*Faox_evap_wiso(3,n)
          end if
       end do
    end if

  end subroutine med_diag_ocn_ocn2med

  !===============================================================================

  subroutine med_diag_ocn_med2ocn(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(in)          :: gcomp
    integer             , intent(out)         :: rc

    ! local variables
    type(InternalState)    :: is_local
    character(*),parameter :: subName = '(med_diag_ocn_ocn2med) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac', ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac', ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBMed_aoflux_o, 'Foxx_lwup', rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Foxx_lat' , rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Foxx_sen' , rc=rc) .and. &
         fldchk(is_local%wrap%FBMed_aoflux_o, 'Foxx_evap', rc=rc)) then

       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_lwup', Foxx_lwup, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_lat', Foxx_lat, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_sen', Foxx_sen, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_evap', Foxx_evap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ip = period_inst
       ic = c_ocn_orecv
       do n = 1,size(Foxx_lwup)
          wgt_o = is_local%wrap%mesh_info(compocn)%areas(n) * ofrac(n)
          wgt_i = is_local%wrap%mesh_info(compocn)%areas(n) * ifrac(n)

          nf = f_hlwup; budget_local(nf,ic,ip) = budget_local(nf,ic,ip) + (wgt_o+wgt_i)*Foxx_lwup(n)
          nf = f_hlatv; budget_local(nf,ic,ip) = budget_local(nf,ic,ip) + (wgt_o+wgt_i)*Foxx_lat(n)
          nf = f_hsen ; budget_local(nf,ic,ip) = budget_local(nf,ic,ip) + (wgt_o+wgt_i)*Foxx_sen(n)
          nf = f_wevap; budget_local(nf,ic,ip) = budget_local(nf,ic,ip) + (wgt_o+wgt_i)*Foxx_evap(n)
       end do
    endif

    call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Fioi_meltw', Fioi_meltw, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Fioi_melth', Fioi_melth, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBExp(compocn), 'Fioi_bergw', rc=rc)) then
       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Fioi_bergw', Fioi_bergw, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ice_bergw = .true.
    end if
    if ( fldchk(is_local%wrap%FBExp(compocn), 'Fioi_bergh', rc=rc)) then
       call FB_getFldPtr(is_local%wrap%FBExp(compocn), 'Fioi_bergh', Fioi_bergh, rc=rc)
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
          budget_local(f_wmelt,ic,ip) = budget_local(nf,ic,ip) + (wgt_o+wgt_i)*Fioi_meltw(n)
       else
          budget_local(f_wmelt,ic,ip) = budget_local(nf,ic,ip) + (wgt_o+wgt_i)*(Fioi_meltw(n)+Fioi_bergw(n))
       endif

       if (.not. ice_bergh) then
          budget_local(f_hmelt,ic,ip) = budget_local(f_hmelt,ic,ip) + (wgt_o+wgt_i)*Fioi_melth(n)
       else
          budget_local(f_hmelt,ic,ip) = budget_local(f_hmelt,ic,ip) + (wgt_o+wgt_i)*(Fioi_melth(n) + Fioi_bergh(n))
       endif

       budget_local(f_wsalt ,ic,ip) = budget_local(f_wsalt ,ic,ip) + (wgt_o+wgt_i)*Fioi_salt(n) * SFLXtoWFLX
       budget_local(f_hswnet,ic,ip) = budget_local(f_hswnet,ic,ip) + (wgt_o+wgt_i)*Foxx_swnet(n)
       budget_local(f_hlwdn ,ic,ip) = budget_local(f_hlwdn ,ic,ip) + (wgt_o+wgt_i)*Faxa_lwdn(n)
       budget_local(f_wrain ,ic,ip) = budget_local(f_wrain ,ic,ip) + (wgt_o+wgt_i)*Faxa_rain(n)
       budget_local(f_wsnow ,ic,ip) = budget_local(f_wsnow ,ic,ip) + (wgt_o+wgt_i)*Faxa_snow(n)
       budget_local(f_wroff ,ic,ip) = budget_local(f_wroff ,ic,ip) + (wgt_o+wgt_i)*Foxx_rofl(n)
       budget_local(f_wioff ,ic,ip) = budget_local(f_wioff ,ic,ip) + (wgt_o+wgt_i)*Foxx_rofi(n)

       if ( ocn_wiso )then
          budget_local(f_wmelt_16O,ic,ip) = budget_local(f_wmelt_16O,ic,ip) + (wgt_o+wgt_i)*Fioi_meltw_wiso(1,n)
          budget_local(f_wmelt_18O,ic,ip) = budget_local(f_wmelt_18O,ic,ip) + (wgt_o+wgt_i)*Fioi_meltw_wiso(2,n)
          budget_local(f_wmelt_HDO,ic,ip) = budget_local(f_wmelt_HDO,ic,ip) + (wgt_o+wgt_i)*Fioi_meltw_wiso(3,n)

          budget_local(f_wrain_16O,ic,ip) = budget_local(f_wrain_16O,ic,ip) + (wgt_o+wgt_i)*Faxa_rain_wiso(1,n)
          budget_local(f_wrain_18O,ic,ip) = budget_local(f_wrain_18O,ic,ip) + (wgt_o+wgt_i)*Faxa_rain_wiso(2,n)
          budget_local(f_wrain_HDO,ic,ip) = budget_local(f_wrain_HDO,ic,ip) + (wgt_o+wgt_i)*Faxa_rain_wiso(3,n)

          budget_local(f_wsnow_16O,ic,ip) = budget_local(f_wsnow_16O,ic,ip) + (wgt_o+wgt_i)*Faxa_snow_wiso(1,n)
          budget_local(f_wsnow_18O,ic,ip) = budget_local(f_wsnow_18O,ic,ip) + (wgt_o+wgt_i)*Faxa_snow_wiso(2,n)
          budget_local(f_wsnow_HDO,ic,ip) = budget_local(f_wsnow_HDO,ic,ip) + (wgt_o+wgt_i)*Faxa_snow_wiso(3,n)

          budget_local(f_wroff_16O,ic,ip) = budget_local(f_wroff_16O,ic,ip) + (wgt_o+wgt_i)*Foxx_rofl_wiso(1,n)
          budget_local(f_wroff_18O,ic,ip) = budget_local(f_wroff_18O,ic,ip) + (wgt_o+wgt_i)*Foxx_rofl_wiso(2,n)
          budget_local(f_wroff_HDO,ic,ip) = budget_local(f_wroff_HDO,ic,ip) + (wgt_o+wgt_i)*Foxx_rofl_wiso(3,n)

          budget_local(f_wioff_16O,ic,ip) = budget_local(f_wioff_16O,ic,ip) + (wgt_o+wgt_i)*Foxx_rofi_wiso(1,n)
          budget_local(f_wioff_18O,ic,ip) = budget_local(f_wioff_18O,ic,ip) + (wgt_o+wgt_i)*Foxx_rofi_wiso(2,n)
          budget_local(f_wioff_HDO,ic,ip) = budget_local(f_wioff_HDO,ic,ip) + (wgt_o+wgt_i)*Foxx_rofi_wiso(3,n)
       end if
    end do
    budget_local(f_hlatf,ic,ip) = -budget_local(f_wsnow,ic,ip)*shr_const_latice
    budget_local(f_hioff,ic,ip) = -budget_local(f_wioff,ic,ip)*shr_const_latice

    ! EBK -- isotope r2x_Forr_rofl/i?

  end subroutine med_diag_ocn_med2ocn

  !===============================================================================

  subroutine med_diag_ice_ice2med( gcomp, rc)

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
    character(*),parameter :: subName = '(med_diag_ice_ice2med) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ifrac', ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ofrac', ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Fioi_melth', Fioi_melth, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Fioi_meltw', Fioi_meltw, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Fioi_salt', Fioi_salt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Fioi_swpen', Fioi_swpen, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_swnet', Faii_swnet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_lwup', Faii_lwup, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_sen', Faii_sen, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_lat', Faii_lat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_evap', Faii_evap, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_meltw_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBImp(compice,compice), 'Faii_evap_wiso' , rc=rc)) then
       ice_wiso = .true.
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Fioi_meltw_wiso', Fioi_meltw_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compatm), 'Faii_evap_wiso', Faii_evap_wiso, rc=rc)
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
       budget_local(f_hmelt ,ic,ip) = budget_local(f_hmelt ,ic,ip) - wgt_i*Fioi_melth(n)
       budget_local(f_wmelt ,ic,ip) = budget_local(f_wmelt ,ic,ip) - wgt_i*Fioi_meltw(n)
       budget_local(f_wsalt ,ic,ip) = budget_local(f_wsalt ,ic,ip) - wgt_i*Fioi_salt(n) * SFLXtoWFLX
       budget_local(f_hswnet,ic,ip) = budget_local(f_hswnet,ic,ip) + wgt_i*Faii_swnet(n) - wgt_i*Fioi_swpen(n)
       budget_local(f_hlwup ,ic,ip) = budget_local(f_hlwup ,ic,ip) + wgt_i*Faii_lwup(n)
       budget_local(f_hlatv ,ic,ip) = budget_local(f_hlatv ,ic,ip) + wgt_i*Faii_lat(n)
       budget_local(f_hsen  ,ic,ip) = budget_local(f_hsen  ,ic,ip) + wgt_i*Faii_sen(n)
       budget_local(f_wevap ,ic,ip) = budget_local(f_wevap ,ic,ip) + wgt_i*Faii_evap(n)

       if ( ice_wiso )then
          budget_local(f_wmelt_16O,ic,ip) = budget_local(f_wmelt_16O,ic,ip) - wgt_i*Fioi_meltw_wiso(1,n)
          budget_local(f_wmelt_18O,ic,ip) = budget_local(f_wmelt_18O,ic,ip) - wgt_i*Fioi_meltw_wiso(2,n)
          budget_local(f_wmelt_HDO,ic,ip) = budget_local(f_wmelt_HDO,ic,ip) - wgt_i*Fioi_meltw_wiso(3,n)

          budget_local(f_wevap_16O,ic,ip) = budget_local(f_wevap_16O,ic,ip) + wgt_i*Faii_evap_wiso(1,n)
          budget_local(f_wevap_18O,ic,ip) = budget_local(f_wevap_18O,ic,ip) + wgt_i*Faii_evap_wiso(2,n)
          budget_local(f_wevap_HDO,ic,ip) = budget_local(f_wevap_HDO,ic,ip) + wgt_i*Faii_evap_wiso(3,n)
       end if
    end do

  end subroutine med_diag_ice_ice2med

  !===============================================================================

  subroutine med_diag_ice_med2ice( gcomp, rc)

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
    character(*),parameter :: subName = '(med_diag_ice_med2ice) '
    ! ------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ifrac', ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ofrac', ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( fldchk(is_local%wrap%FBExp(compice), 'Faxa_rain_wiso', rc=rc) .and. &
         fldchk(is_local%wrap%FBExp(compice), 'Faxa_snow_wiso', rc=rc)) then
       ice_wiso = .true.
       call FB_getFldPtr(is_local%wrap%FBExp(compice), 'Faxa_rain_wiso', Faxa_rain_wiso, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBExp(compice), 'Faxa_snow_wiso', Faxa_snow_wiso, rc=rc)
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
       budget_local(f_hlwdn,ic,ip) = budget_local(f_hlwdn,ic,ip) + wgt_i*Faxa_lwdn(n)
       budget_local(f_wrain,ic,ip) = budget_local(f_wrain,ic,ip) + wgt_i*Faxa_rain(n)
       budget_local(f_wsnow,ic,ip) = budget_local(f_wsnow,ic,ip) + wgt_i*Faxa_snow(n)
       budget_local(f_wioff,ic,ip) = budget_local(f_wioff,ic,ip) + wgt_i*Fixx_rofi(n)
       budget_local(f_hfrz ,ic,ip) = budget_local(f_hfrz ,ic,ip) - (wgt_o + wgt_i)*max(0.0_r8,Fioo_q(n))

       if ( ice_wiso ) then
          budget_local(f_wrain_16O,ic,ip) = budget_local(f_wrain_16O,ic,ip) + wgt_i*Faxa_rain_wiso(1,n)
          budget_local(f_wrain_18O,ic,ip) = budget_local(f_wrain_18O,ic,ip) + wgt_i*Faxa_rain_wiso(2,n)
          budget_local(f_wrain_HDO,ic,ip) = budget_local(f_wrain_HDO,ic,ip) + wgt_i*Faxa_rain_wiso(3,n)

          budget_local(f_wsnow_16O,ic,ip) = budget_local(f_wsnow_16O,ic,ip) + wgt_i*Faxa_snow_wiso(1,n)
          budget_local(f_wsnow_18O,ic,ip) = budget_local(f_wsnow_18O,ic,ip) + wgt_i*Faxa_snow_wiso(2,n)
          budget_local(f_wsnow_HDO,ic,ip) = budget_local(f_wsnow_HDO,ic,ip) + wgt_i*Faxa_snow_wiso(3,n)
       end if
    end do

    ic = c_inh_isend
    budget_local(f_hlatf,ic,ip) = -budget_local(f_wsnow,ic,ip)*shr_const_latice
    budget_local(f_hioff,ic,ip) = -budget_local(f_wioff,ic,ip)*shr_const_latice
    budget_local(f_wfrz ,ic,ip) =  budget_local(f_hfrz ,ic,ip)*HFLXtoWFLX

    ic = c_ish_isend
    budget_local(f_hlatf,ic,ip) = -budget_local(f_wsnow,ic,ip)*shr_const_latice
    budget_local(f_hioff,ic,ip) = -budget_local(f_wioff,ic,ip)*shr_const_latice
    budget_local(f_wfrz ,ic,ip) =  budget_local(f_hfrz ,ic,ip)*HFLXtoWFLX

  end subroutine med_diag_ice_med2ice

  !===============================================================================

  subroutine med_diag_print(gcomp, rc=rc)

    ! ------------------------------------------------------------------
    ! Print global budget diagnostics.
    ! ------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp), intent(in) :: gcomp
    logical          , intent(in) :: stop_alarm
    integer          , intent(out):: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Alarm)            :: stop_alarm
    integer                     :: budg_print_inst
    integer                     :: budg_print_daily
    integer                     :: budg_print_month
    integer                     :: budg_print_ann
    integer                     :: budg_print_ltann
    integer                     :: budg_print_ltend
    integer                     :: ic,nf,ip,is                           ! data array indicies
    integer                     :: ica,icl,icn,ics,ico
    integer                     :: icar,icxs,icxr,icas
    integer                     :: cdate,sec                             ! coded date, seconds
    integer                     :: yr,mon,day                            ! date
    integer                     :: iam                                   ! pe number
    integer                     :: period_type                           ! print level
    logical                     :: sumdone                               ! has a sum been computed yet
    character(len=40)           :: str                                   ! string
    real(r8)                    :: dataGpr (f_size,c_size,period_ntypes) ! values to print, scaled and such
    integer, parameter          :: nisotopes = 3
    character(len=5), parameter :: isoname(nisotopes) = (/ 'H216O',   'H218O',   '  HDO'   /)
    integer, parameter          :: iso0(nisotopes)    = (/ f_16O,     f_18O,     f_hdO     /)
    integer, parameter          :: isof(nisotopes)    = (/ f_16O_end, f_18O_end, f_hdO_end /)
    character(*),parameter      :: F00   = "('(med_diag_print) ',4a)"
    character(*),parameter      :: FAH="(4a,i9,i6)"
    character(*),parameter      :: FA0= "('    ',12x,6(6x,a8,1x))"
    character(*),parameter      :: FA1= "('    ',a12,6f15.8)"
    character(*),parameter      :: FA0r="('    ',12x,8(6x,a8,1x))"
    character(*),parameter      :: FA1r="('    ',a12,8f15.8)"
    character(*),parameter      :: subName = '(med_diag_print) '
    ! ------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    ! print instantaneous budget data
    !-------------------------------------------------------------------------------

    sumdone = .false.

    call NUOPC_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=curr_year, mm=curr_mon, dd=curr_day, s=curr_tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    cdate = curr_yr*10000 + curr_mon*100 + curr_day

    do ip = 1,period_ntypes
       period_type = 0
       if (ip == period_inst) then
          period_type = max(period_type,budg_print_inst)
       endif
       if (ip==period_curr_day .and. curr_tod==0) then
          period_type = max(period_type,budg_print_daily)
       endif
       if (ip==period_mon .and. curr_day==1 .and. curr_tod==0) then
          period_type = max(period_type,budg_print_month)
       endif
       if (ip==period_ann .and. mon==1 .and. curr_day==1 .and. curr_tod==0) then
          period_type = max(period_type,budg_print_ann)
       endif
       if (ip==period_inf .and. mon==1 .and. curr_day==1 .and. curr_tod==0) then
          period_type = max(period_type,budg_print_ltann)
       endif
       if (ip==period_inf .and. stop_alarm) then
          period_type = max(period_type,budg_print_ltend)
       endif

       if (period_type > 0) then
          ! ---- doprint ---- doprint ---- doprint ----
          if (.not.sumdone) then

             call med_diag_sum0()
             dataGpr = budget_global
             sumdone = .true.

             ! old budget normalizations (global area and 1e6 for water)
             dataGpr = dataGpr/(4.0_r8*shr_const_pi)
             dataGpr(f_w:f_w_end,:,:) = dataGpr(f_w:f_w_end,:,:) * 1.0e6_r8
             if ( flds_wiso ) then
                dataGpr(iso0(1):isof(nisotopes),:,:) = dataGpr(iso0(1):isof(nisotopes),:,:) * 1.0e6_r8
             end if
             dataGpr = dataGpr/budget_ns

             if (.not. masterproc) return
          endif

          ! ---------------------------------------------------------
          ! ---- detail atm budgets and breakdown into components ---
          ! ---------------------------------------------------------

          if (period_type >= 3) then
             do ic = 1,2
                if (ic == 1) then    ! from atm to mediator
                   ica = c_atm_arecv ! total from atm
                   icl = c_lnd_arecv ! sent to land from med on atm grid
                   icn = c_inh_arecv ! sent to ice nh from med on atm grid
                   ics = c_ish_arecv ! sent to ice sh from med on atm grid
                   ico = c_ocn_arecv ! sent to ocean from med on atm grid
                   str = "ATM_to_CPL"
                elseif (ic == 2) then ! from mediator to atm
                   ica = c_atm_asend  ! merged on atm grid
                   icl = c_lnd_asend  ! from land on atm grid
                   icn = c_inh_asend  ! from ice nh on atm grid
                   ics = c_ish_asend  ! from ice sh on atm grid
                   ico = c_ocn_asend  ! from ocn on atm grid
                   str = "CPL_TO_ATM"
                else
                   call shr_sys_abort(subname//' ERROR in ic index code 411')
                endif

                write(logunit,*) ' '
                write(logunit,FAH) subname,trim(str)//' AREA BUDGET (m2/m2): period = ',trim(pname(ip)),&
                     ': date = ', cdate, curr_tod
                write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
                do nf = f_a, f_a_end
                   write(logunit,FA1) fname(nf),&
                        dataGpr(nf,ica,ip), &
                        dataGpr(nf,icl,ip), &
                        dataGpr(nf,icn,ip), &
                        dataGpr(nf,ics,ip), &
                        dataGpr(nf,ico,ip), &
                        dataGpr(nf,ica,ip) + dataGpr(nf,icl,ip) + &
                        dataGpr(nf,icn,ip) + dataGpr(nf,ics,ip) + dataGpr(nf,ico,ip)
                enddo

                write(logunit,*) ' '
                write(logunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,curr_tod
                write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
                do nf = f_h, f_h_end
                   write(logunit,FA1) fname(nf),&
                        dataGpr(nf,ica,ip), &
                        dataGpr(nf,icl,ip), &
                        dataGpr(nf,icn,ip), &
                        dataGpr(nf,ics,ip), &
                        dataGpr(nf,ico,ip), &
                        dataGpr(nf,ica,ip) + dataGpr(nf,icl,ip) + &
                        dataGpr(nf,icn,ip) + dataGpr(nf,ics,ip) + dataGpr(nf,ico,ip)
                enddo
                write(logunit,FA1)    '   *SUM*'   ,&
                     sum(dataGpr(f_h:f_h_end,ica,ip)), &
                     sum(dataGpr(f_h:f_h_end,icl,ip)), &
                     sum(dataGpr(f_h:f_h_end,icn,ip)), &
                     sum(dataGpr(f_h:f_h_end,ics,ip)), &
                     sum(dataGpr(f_h:f_h_end,ico,ip)), &
                     sum(dataGpr(f_h:f_h_end,ica,ip)) + sum(dataGpr(f_h:f_h_end,icl,ip)) + &
                     sum(dataGpr(f_h:f_h_end,icn,ip)) + sum(dataGpr(f_h:f_h_end,ics,ip)) + sum(dataGpr(f_h:f_h_end,ico,ip))

                write(logunit,*) ' '
                write(logunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,curr_tod
                write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
                do nf = f_w, f_w_end
                   write(logunit,FA1) fname(nf),&
                        dataGpr(nf,ica,ip), &
                        dataGpr(nf,icl,ip), &
                        dataGpr(nf,icn,ip), &
                        dataGpr(nf,ics,ip), &
                        dataGpr(nf,ico,ip), &
                        dataGpr(nf,ica,ip) + dataGpr(nf,icl,ip) + &
                        dataGpr(nf,icn,ip) + dataGpr(nf,ics,ip) + dataGpr(nf,ico,ip)
                enddo
                write(logunit,FA1)    '   *SUM*'   ,&
                     sum(dataGpr(f_w:f_w_end,ica,ip)), &
                     sum(dataGpr(f_w:f_w_end,icl,ip)), &
                     sum(dataGpr(f_w:f_w_end,icn,ip)), &
                     sum(dataGpr(f_w:f_w_end,ics,ip)), &
                     sum(dataGpr(f_w:f_w_end,ico,ip)), &
                     sum(dataGpr(f_w:f_w_end,ica,ip)) + sum(dataGpr(f_w:f_w_end,icl,ip)) + &
                     sum(dataGpr(f_w:f_w_end,icn,ip)) + sum(dataGpr(f_w:f_w_end,ics,ip)) + sum(dataGpr(f_w:f_w_end,ico,ip))

                if ( flds_wiso ) then
                   do is = 1, nisotopes
                      write(logunit,*) ' '
                      write(logunit,FAH) subname,trim(str)//' '//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ', &
                           trim(pname(ip)),': date = ',cdate,curr_tod
                      write(logunit,FA0) cname(ica),cname(icl),cname(icn),cname(ics),cname(ico),' *SUM*  '
                      do nf = iso0(is), isof(is)
                         write(logunit,FA1)    fname(nf),&
                              dataGpr(nf,ica,ip), &
                              dataGpr(nf,icl,ip), &
                              dataGpr(nf,icn,ip), &
                              dataGpr(nf,ics,ip), &
                              dataGpr(nf,ico,ip), &
                              dataGpr(nf,ica,ip) + dataGpr(nf,icl,ip) + &
                              dataGpr(nf,icn,ip) + dataGpr(nf,ics,ip) + dataGpr(nf,ico,ip)
                      enddo
                      write(logunit,FA1)    '   *SUM*', &
                           sum(dataGpr(iso0(is):isof(is),ica,ip)), &
                           sum(dataGpr(iso0(is):isof(is),icl,ip)), &
                           sum(dataGpr(iso0(is):isof(is),icn,ip)), &
                           sum(dataGpr(iso0(is):isof(is),ics,ip)), &
                           sum(dataGpr(iso0(is):isof(is),ico,ip)), &
                           sum(dataGpr(iso0(is):isof(is),ica,ip)) + sum(dataGpr(iso0(is):isof(is),icl,ip)) + &
                           sum(dataGpr(iso0(is):isof(is),icn,ip)) + sum(dataGpr(iso0(is):isof(is),ics,ip)) + &
                           sum(dataGpr(iso0(is):isof(is),ico,ip))
                   end do
                end if

             enddo
          endif   ! period_type

          ! ---------------------------------------------------------
          ! ---- detail lnd/ocn/ice component budgets ----
          ! ---------------------------------------------------------

          if (period_type /= period_inst) then
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

                write(logunit,*) ' '
                write(logunit,FAH) subname,trim(str)//' HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,curr_tod
                write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
                do nf = f_h, f_h_end
                   write(logunit,FA1) fname(nf),&
                        -dataGpr(nf,icar,ip), &
                         dataGpr(nf,icxs,ip), &
                         dataGpr(nf,icxr,ip), &
                        -dataGpr(nf,icas,ip), &
                        -dataGpr(nf,icar,ip) + dataGpr(nf,icxs,ip) + dataGpr(nf,icxr,ip) - dataGpr(nf,icas,ip)
                enddo
                write(logunit,FA1)'   *SUM*',&
                     -sum(dataGpr(f_h:f_h_end,icar,ip)), &
                      sum(dataGpr(f_h:f_h_end,icxs,ip)), &
                      sum(dataGpr(f_h:f_h_end,icxr,ip)), &
                     -sum(dataGpr(f_h:f_h_end,icas,ip)), &
                     -sum(dataGpr(f_h:f_h_end,icar,ip)) + sum(dataGpr(f_h:f_h_end,icxs,ip)) + &
                      sum(dataGpr(f_h:f_h_end,icxr,ip)) - sum(dataGpr(f_h:f_h_end,icas,ip))

                write(logunit,*) ' '
                write(logunit,FAH) subname,trim(str)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,curr_tod
                write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
                do nf = f_w, f_w_end
                   write(logunit,FA1)  fname(nf),&
                        -dataGpr(nf,icar,ip),&
                         dataGpr(nf,icxs,ip), &
                         dataGpr(nf,icxr,ip),&
                        -dataGpr(nf,icas,ip), &
                        -dataGpr(nf,icar,ip) + dataGpr(nf,icxs,ip) + dataGpr(nf,icxr,ip) - dataGpr(nf,icas,ip)
                enddo
                write(logunit,FA1)    '   *SUM*',&
                     -sum(dataGpr(f_w:f_w_end,icar,ip)), &
                      sum(dataGpr(f_w:f_w_end,icxs,ip)), &
                      sum(dataGpr(f_w:f_w_end,icxr,ip)), &
                     -sum(dataGpr(f_w:f_w_end,icas,ip)), &
                     -sum(dataGpr(f_w:f_w_end,icar,ip)) + sum(dataGpr(f_w:f_w_end,icxs,ip)) + &
                      sum(dataGpr(f_w:f_w_end,icxr,ip)) - sum(dataGpr(f_w:f_w_end,icas,ip))

                if ( flds_wiso ) then
                   do is = 1, nisotopes
                      write(logunit,*) ' '
                      write(logunit,FAH) subname,trim(str)//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)), &
                           ': date = ',cdate,curr_tod
                      write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
                      do nf = iso0(is), isof(is)
                         write(logunit,FA1) fname(nf),&
                              -dataGpr(nf,icar,ip), &
                               dataGpr(nf,icxs,ip), &
                               dataGpr(nf,icxr,ip), &
                              -dataGpr(nf,icas,ip), &
                              -dataGpr(nf,icar,ip) + dataGpr(nf,icxs,ip) + dataGpr(nf,icxr,ip) - dataGpr(nf,icas,ip)
                      enddo
                      write(logunit,FA1)    '   *SUM*',-sum(dataGpr(iso0(is):isof(is),icar,ip)),&
                            sum(dataGpr(iso0(is):isof(is),icxs,ip)), &
                            sum(dataGpr(iso0(is):isof(is),icxr,ip)), &
                           -sum(dataGpr(iso0(is):isof(is),icas,ip)), &
                           -sum(dataGpr(iso0(is):isof(is),icar,ip)) + sum(dataGpr(iso0(is):isof(is),icxs,ip)) + &
                            sum(dataGpr(iso0(is):isof(is),icxr,ip)) - sum(dataGpr(iso0(is):isof(is),icas,ip))
                      write(logunit,*) ' '
                      write(logunit,FAH) subname,trim(str)//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),&
                           ': date = ',cdate,curr_tod
                      write(logunit,FA0) cname(icar),cname(icxs),cname(icxr),cname(icas),' *SUM*  '
                      do nf = iso0(is), isof(is)
                         write(logunit,FA1) fname(nf),&
                              -dataGpr(nf,icar,ip), &
                               dataGpr(nf,icxs,ip), &
                               dataGpr(nf,icxr,ip), &
                              -dataGpr(nf,icas,ip), &
                              -dataGpr(nf,icar,ip) + dataGpr(nf,icxs,ip) + dataGpr(nf,icxr,ip) - dataGpr(nf,icas,ip)
                      enddo
                      write(logunit,FA1)    '   *SUM*',                &
                           -sum(dataGpr(iso0(is):isof(is), icar, ip)), &
                            sum(dataGpr(iso0(is):isof(is), icxs, ip)), &
                            sum(dataGpr(iso0(is):isof(is), icxr, ip)), &
                           -sum(dataGpr(iso0(is):isof(is), icas, ip)), &
                           -sum(dataGpr(iso0(is):isof(is), icar, ip)) + sum(dataGpr(iso0(is):isof(is), icxs, ip)) + &
                            sum(dataGpr(iso0(is):isof(is), icxr, ip)) - sum(dataGpr(iso0(is):isof(is), icas, ip))
                   end do
                end if
             enddo
          endif   ! period_type

          ! ---------------------------------------------------------
          ! ---- net summary budgets ----
          ! ---------------------------------------------------------

          if (period_type >= 1) then

             write(logunit,*) ' '
             write(logunit,FAH) subname,'NET AREA BUDGET (m2/m2): period = ',trim(pname(ip)),': date = ',cdate,curr_tod
             write(logunit,FA0) '     atm','     lnd','     ocn','  ice nh','  ice sh',' *SUM*  '
             do nf = f_a,f_a_end
                write(logunit,FA1) fname(nf),&
                     dataGpr(nf,c_atm_arecv,ip), &
                     dataGpr(nf,c_lnd_lrecv,ip), &
                     dataGpr(nf,c_ocn_orecv,ip), &
                     dataGpr(nf,c_inh_irecv,ip), &
                     dataGpr(nf,c_ish_irecv,ip), &
                     dataGpr(nf,c_atm_arecv,ip)+ &
                     dataGpr(nf,c_lnd_lrecv,ip)+ &
                     dataGpr(nf,c_ocn_orecv,ip)+ &
                     dataGpr(nf,c_inh_irecv,ip)+ &
                     dataGpr(nf,c_ish_irecv,ip)
             enddo

             write(logunit,*) ' '
             write(logunit,FAH) subname,'NET HEAT BUDGET (W/m2): period = ',trim(pname(ip)),': date = ',cdate,curr_tod
             write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
             do nf = f_h, f_h_end
                write(logunit,FA1r)   fname(nf),&
                     dataGpr(nf,c_atm_arecv,ip)+dataGpr(nf,c_atm_asend,ip), &
                     dataGpr(nf,c_lnd_lrecv,ip)+dataGpr(nf,c_lnd_lsend,ip), &
                     dataGpr(nf,c_rof_rrecv,ip)+dataGpr(nf,c_rof_rsend,ip), &
                     dataGpr(nf,c_ocn_orecv,ip)+dataGpr(nf,c_ocn_osend,ip), &
                     dataGpr(nf,c_inh_irecv,ip)+dataGpr(nf,c_inh_isend,ip), &
                     dataGpr(nf,c_ish_irecv,ip)+dataGpr(nf,c_ish_isend,ip), &
                     dataGpr(nf,c_glc_grecv,ip)+dataGpr(nf,c_glc_gsend,ip), &
                     dataGpr(nf,c_atm_arecv,ip)+dataGpr(nf,c_atm_asend,ip)+ &
                     dataGpr(nf,c_lnd_lrecv,ip)+dataGpr(nf,c_lnd_lsend,ip)+ &
                     dataGpr(nf,c_rof_rrecv,ip)+dataGpr(nf,c_rof_rsend,ip)+ &
                     dataGpr(nf,c_ocn_orecv,ip)+dataGpr(nf,c_ocn_osend,ip)+ &
                     dataGpr(nf,c_inh_irecv,ip)+dataGpr(nf,c_inh_isend,ip)+ &
                     dataGpr(nf,c_ish_irecv,ip)+dataGpr(nf,c_ish_isend,ip)+ &
                     dataGpr(nf,c_glc_grecv,ip)+dataGpr(nf,c_glc_gsend,ip)
             enddo
             write(logunit,FA1r)'   *SUM*',&
                  sum(dataGpr(f_h:f_h_end,c_atm_arecv,ip))+sum(dataGpr(f_h:f_h_end,c_atm_asend,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_lnd_lrecv,ip))+sum(dataGpr(f_h:f_h_end,c_lnd_lsend,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_rof_rrecv,ip))+sum(dataGpr(f_h:f_h_end,c_rof_rsend,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_ocn_orecv,ip))+sum(dataGpr(f_h:f_h_end,c_ocn_osend,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_inh_irecv,ip))+sum(dataGpr(f_h:f_h_end,c_inh_isend,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_ish_irecv,ip))+sum(dataGpr(f_h:f_h_end,c_ish_isend,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_glc_grecv,ip))+sum(dataGpr(f_h:f_h_end,c_glc_gsend,ip)), &
                  sum(dataGpr(f_h:f_h_end,c_atm_arecv,ip))+sum(dataGpr(f_h:f_h_end,c_atm_asend,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_lnd_lrecv,ip))+sum(dataGpr(f_h:f_h_end,c_lnd_lsend,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_rof_rrecv,ip))+sum(dataGpr(f_h:f_h_end,c_rof_rsend,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_ocn_orecv,ip))+sum(dataGpr(f_h:f_h_end,c_ocn_osend,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_inh_irecv,ip))+sum(dataGpr(f_h:f_h_end,c_inh_isend,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_ish_irecv,ip))+sum(dataGpr(f_h:f_h_end,c_ish_isend,ip))+ &
                  sum(dataGpr(f_h:f_h_end,c_glc_grecv,ip))+sum(dataGpr(f_h:f_h_end,c_glc_gsend,ip))

             write(logunit,*) ' '
             write(logunit,FAH) subname,'NET WATER BUDGET (kg/m2s*1e6): period = ',trim(pname(ip)),': date = ',cdate,curr_tod
             write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
             do nf = f_w, f_w_end
                write(logunit,FA1r)   fname(nf),&
                     dataGpr(nf,c_atm_arecv,ip)+dataGpr(nf,c_atm_asend,ip), &
                     dataGpr(nf,c_lnd_lrecv,ip)+dataGpr(nf,c_lnd_lsend,ip), &
                     dataGpr(nf,c_rof_rrecv,ip)+dataGpr(nf,c_rof_rsend,ip), &
                     dataGpr(nf,c_ocn_orecv,ip)+dataGpr(nf,c_ocn_osend,ip), &
                     dataGpr(nf,c_inh_irecv,ip)+dataGpr(nf,c_inh_isend,ip), &
                     dataGpr(nf,c_ish_irecv,ip)+dataGpr(nf,c_ish_isend,ip), &
                     dataGpr(nf,c_glc_grecv,ip)+dataGpr(nf,c_glc_gsend,ip), &
                     dataGpr(nf,c_atm_arecv,ip)+dataGpr(nf,c_atm_asend,ip)+ &
                     dataGpr(nf,c_lnd_lrecv,ip)+dataGpr(nf,c_lnd_lsend,ip)+ &
                     dataGpr(nf,c_rof_rrecv,ip)+dataGpr(nf,c_rof_rsend,ip)+ &
                     dataGpr(nf,c_ocn_orecv,ip)+dataGpr(nf,c_ocn_osend,ip)+ &
                     dataGpr(nf,c_inh_irecv,ip)+dataGpr(nf,c_inh_isend,ip)+ &
                     dataGpr(nf,c_ish_irecv,ip)+dataGpr(nf,c_ish_isend,ip)+ &
                     dataGpr(nf,c_glc_grecv,ip)+dataGpr(nf,c_glc_gsend,ip)
             enddo
             write(logunit,FA1r)'   *SUM*',&
                  sum(dataGpr(f_w:f_w_end,c_atm_arecv,ip))+sum(dataGpr(f_w:f_w_end,c_atm_asend,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_lnd_lrecv,ip))+sum(dataGpr(f_w:f_w_end,c_lnd_lsend,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_rof_rrecv,ip))+sum(dataGpr(f_w:f_w_end,c_rof_rsend,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_ocn_orecv,ip))+sum(dataGpr(f_w:f_w_end,c_ocn_osend,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_inh_irecv,ip))+sum(dataGpr(f_w:f_w_end,c_inh_isend,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_ish_irecv,ip))+sum(dataGpr(f_w:f_w_end,c_ish_isend,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_glc_grecv,ip))+sum(dataGpr(f_w:f_w_end,c_glc_gsend,ip)), &
                  sum(dataGpr(f_w:f_w_end,c_atm_arecv,ip))+sum(dataGpr(f_w:f_w_end,c_atm_asend,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_lnd_lrecv,ip))+sum(dataGpr(f_w:f_w_end,c_lnd_lsend,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_rof_rrecv,ip))+sum(dataGpr(f_w:f_w_end,c_rof_rsend,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_ocn_orecv,ip))+sum(dataGpr(f_w:f_w_end,c_ocn_osend,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_inh_irecv,ip))+sum(dataGpr(f_w:f_w_end,c_inh_isend,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_ish_irecv,ip))+sum(dataGpr(f_w:f_w_end,c_ish_isend,ip))+ &
                  sum(dataGpr(f_w:f_w_end,c_glc_grecv,ip))+sum(dataGpr(f_w:f_w_end,c_glc_gsend,ip))

             if ( flds_wiso ) then

                do is = 1, nisotopes
                   write(logunit,*) ' '
                   write(logunit,FAH) subname,'NET '//isoname(is)//' WATER BUDGET (kg/m2s*1e6): period = ', &
                        trim(pname(ip)),': date = ',cdate,curr_tod
                   write(logunit,FA0r) '     atm','     lnd','     rof','     ocn','  ice nh','  ice sh','     glc',' *SUM*  '
                   do nf = iso0(is), isof(is)
                      write(logunit,FA1r)   fname(nf),&
                           dataGpr(nf,c_atm_arecv,ip)+dataGpr(nf,c_atm_asend,ip), &
                           dataGpr(nf,c_lnd_lrecv,ip)+dataGpr(nf,c_lnd_lsend,ip), &
                           dataGpr(nf,c_rof_rrecv,ip)+dataGpr(nf,c_rof_rsend,ip), &
                           dataGpr(nf,c_ocn_orecv,ip)+dataGpr(nf,c_ocn_osend,ip), &
                           dataGpr(nf,c_inh_irecv,ip)+dataGpr(nf,c_inh_isend,ip), &
                           dataGpr(nf,c_ish_irecv,ip)+dataGpr(nf,c_ish_isend,ip), &
                           dataGpr(nf,c_glc_grecv,ip)+dataGpr(nf,c_glc_gsend,ip), &
                           dataGpr(nf,c_atm_arecv,ip)+dataGpr(nf,c_atm_asend,ip)+ &
                           dataGpr(nf,c_lnd_lrecv,ip)+dataGpr(nf,c_lnd_lsend,ip)+ &
                           dataGpr(nf,c_rof_rrecv,ip)+dataGpr(nf,c_rof_rsend,ip)+ &
                           dataGpr(nf,c_ocn_orecv,ip)+dataGpr(nf,c_ocn_osend,ip)+ &
                           dataGpr(nf,c_inh_irecv,ip)+dataGpr(nf,c_inh_isend,ip)+ &
                           dataGpr(nf,c_ish_irecv,ip)+dataGpr(nf,c_ish_isend,ip)+ &
                           dataGpr(nf,c_glc_grecv,ip)+dataGpr(nf,c_glc_gsend,ip)
                   enddo
                   write(logunit,FA1r)'   *SUM*',&
                        sum(dataGpr(iso0(is):isof(is), c_atm_arecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_atm_asend, ip)), &
                        sum(dataGpr(iso0(is):isof(is), c_lnd_lrecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_lnd_lsend, ip)), &
                        sum(dataGpr(iso0(is):isof(is), c_rof_rrecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_rof_rsend, ip)), &
                        sum(dataGpr(iso0(is):isof(is), c_ocn_orecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_ocn_osend, ip)), &
                        sum(dataGpr(iso0(is):isof(is), c_inh_irecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_inh_isend, ip)), &
                        sum(dataGpr(iso0(is):isof(is), c_ish_irecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_ish_isend, ip)), &
                        sum(dataGpr(iso0(is):isof(is), c_glc_grecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_glc_gsend, ip)), &
                        sum(dataGpr(iso0(is):isof(is), c_atm_arecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_atm_asend, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_lnd_lrecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_lnd_lsend, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_rof_rrecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_rof_rsend, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_ocn_orecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_ocn_osend, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_inh_irecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_inh_isend, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_ish_irecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_ish_isend, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_glc_grecv, ip))+ &
                        sum(dataGpr(iso0(is):isof(is), c_glc_gsend, ip))
                end do
             end if

          endif

          write(logunit,*) ' '
          ! ---- doprint ---- doprint ---- doprint ----
       endif  ! period_type > 0
    enddo  ! ip = 1,period_ntypes

  end subroutine med_diag_print

end module med_diag
