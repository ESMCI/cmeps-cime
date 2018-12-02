module med_phases_prep_glc_mod

  !-----------------------------------------------------------------------------
  ! Mediator Phase - prepare glc input
  !-----------------------------------------------------------------------------

#include "shr_assert.h"

  use ESMF                  , only : ESMF_FieldBundle
  use med_constants_mod     , only : R8, CS, czero => med_constants_czero
  use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
  use esmFlds               , only : compglc, complnd, compname
  use shr_sys_mod           , only : shr_sys_abort
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use glc_elevclass_mod     , only : glc_all_elevclass_strings, GLC_ELEVCLASS_STRLEN
  use glc_elevclass_mod     , only : glc_elevclass_as_string
  use shr_nuopc_fldList_mod , only : mapconsf, mapbilnr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_accum
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFldPtr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFieldN
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_average
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_clean
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use med_internalstate_mod , only : InternalState, mastertask, logunit
  use med_infodata_mod      , only : med_infodata_set_valid_glc_input, med_infodata
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_phases_prep_glc_avg
  public  :: med_phases_prep_glc_accum_fast

  private :: med_phases_prep_glc_accum_init

  ! glc fields with multiple elevation classes: lnd->glc
  ! - fields sent from lnd->med to glc    ARE     IN multiple elevation classes
  ! - fields sent from med->glc from land ARE NOT IN multiple elevation classes
  ! Need to keep track of the lnd->med fields destined for glc in the
  ! FBLndAccum field(:) bundles. 
  ! The size of the FBLndAccum(:) field bundle array will be the same size as 
  ! the number of fields in is_local%wrap%FBExp(glc)
  ! For each field bundle in FBLandAccum(:), the number of fields in that field bundle
  ! will equal the number of elevation classes coming from the land

  type(ESMF_FieldBundle), allocatable :: FBLndAccum(:)     ! Accumulator for various components import
  integer               , allocatable :: FBLndAccumCnt(:)  ! Accumulator counter for each FBLndAccum
  character(len=CS)     , allocatable :: fieldNameList(:)  ! Field names without elevation class suffix for each field bundle
  integer                             :: fieldCount        ! Size of fieldNameList
  integer                             :: nEC               ! Number of elevation classes 

  ! Whether to renormalize the SMB for conservation.
  ! Should be set to true for 2-way coupled runs with evolving ice sheets.
  ! Does not need to be true for 1-way coupling.
  logical :: smb_renormalize

  ! Feild names without elevation class suffix
  character(len=*), parameter :: qice_fieldname   = 'Flgl_qice' ! Name of flux field giving surface mass balance
  character(len=*), parameter :: Sg_frac_field    = 'Sg_ice_covered'
  character(len=*), parameter :: Sg_topo_field    = 'Sg_topo'
  character(len=*), parameter :: Sg_icemask_field = 'Sg_icemask'

  character(*), parameter :: u_FILE_u = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
  
  subroutine med_phases_prep_glc_accum_init(gcomp, rc)

    ! Initialize module field bundles needed for accumulation

    ! uses
    use ESMF              , only : ESMF_GridComp
    use ESMF              , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF              , only : ESMF_FieldBundleGet
    use glc_elevclass_mod , only : glc_elevclass_as_string
    use glc_elevclass_mod , only : glc_get_num_elevation_classes

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)            :: is_local
    integer                        :: nfld
    integer                        :: ec
    logical                        :: glc_coupled_fluxes
    logical                        :: lnd_prognostic
    character(len=CS)              :: glc_renormalize_smb
    character(len=CS)              :: do_renormalize_smb
    character(len=CS), allocatable :: accum_fields(:,:)
    character(len=:) , allocatable :: elevclass_as_string
    integer                        :: dbrc
    character(len=*) , parameter   :: subname='(init_lnd2glc_accum)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! -- Determine FBlndAccum(:) and FBLndAccumCnt(:)
    !---------------------------------------

    if (mastertask) then
       write(logunit,*) subname,' initializing accumulated FBs for '//trim(compname(complnd))
    end if

    ! TODO: make sure cpl_scalar is not used below

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fieldCount=fieldCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldNameList(fieldCount))
    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fieldNameList=fieldNamelist, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine number of elevation classes
    nEC = glc_get_num_elevation_classes()

    allocate(accum_fields(nEC, fieldCount))
    do nfld = 1, fieldCount 
       do ec = 1, nEC
          elevclass_as_string = glc_elevclass_as_string(ec)
          accum_fields(ec,nfld) = trim(fieldnameList(nfld)) // trim(elevclass_as_string)
       end do
    end do

    allocate(FBLndAccum(fieldCount))
    allocate(FBLndAccumCnt(fieldCount))
    do nfld = 1,fieldCount
       call shr_nuopc_methods_FB_init(FBLndAccum(nfld), flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(complnd), fieldNameList=accum_fields(:,nfld), name='FBLndAccum', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_reset(FBLndAccum(nfld), value=czero, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       FBLndAccumCnt(nfld) = 0
    end do

    deallocate(accum_fields)

    !---------------------------------------
    ! Determine if will do smb renormalization
    !---------------------------------------

    ! For now hardwire these
    ! TODO (mvertens, 2018-11-25) : put in the correct logic for the following 
    ! for now simply set them to false

    glc_renormalize_smb = 'off'
    glc_coupled_fluxes  = .false.
    lnd_prognostic      = .false.
    smb_renormalize = prep_glc_do_renormalize_smb(glc_renormalize_smb, glc_coupled_fluxes, lnd_prognostic)

  end subroutine med_phases_prep_glc_accum_init

  !================================================================================================

  function prep_glc_do_renormalize_smb(glc_renormalize_smb, glc_coupled_fluxes, lnd_prognostic) &
       result(do_renormalize_smb)
    
    ! Returns a logical saying whether we should do the smb renormalization

    ! function return
    logical :: do_renormalize_smb   ! function return value

    ! input/output variables
    character(len=*), intent(in) :: glc_renormalize_smb  ! namelist option saying whether to do smb renormalization
    logical         , intent(in) :: glc_coupled_fluxes   ! does glc send fluxes to other components?
    logical         , intent(in) :: lnd_prognostic       ! is lnd a prognostic component?

    ! local variables
    character(len=*), parameter :: subname = '(prep_glc_do_renormalize_smb)'
    !---------------------------------------------------------------

    select case (glc_renormalize_smb)
    case ('on')
       do_renormalize_smb = .true.

    case ('off')
       do_renormalize_smb = .false.

    case ('on_if_glc_coupled_fluxes')
       if (.not. lnd_prognostic) then
          ! Do not renormalize if running glc with dlnd (T compsets): In this case
          ! there is no feedback from glc to lnd, and conservation is not important
          do_renormalize_smb = .false.
       else if (.not. glc_coupled_fluxes) then
          ! Do not renormalize if glc does not send fluxes to other components: In this
          ! case conservation is not important
          do_renormalize_smb = .false.
       else
          ! lnd_prognostic is true and glc_coupled_fluxes is true
          do_renormalize_smb = .true.
       end if

    case default
       write(logunit,*) subname,' ERROR: unknown value for glc_renormalize_smb: ', &
            trim(glc_renormalize_smb)
       call shr_sys_abort(subname//' ERROR: unknown value for glc_renormalize_smb')

    end select

  end function prep_glc_do_renormalize_smb

  !================================================================================================

  subroutine med_phases_prep_glc_accum_fast(gcomp, rc)

    ! Carry out fast accumulation for the land-ice (glc) component
    ! Accumulation and averaging is done on the land input to the river component on the land grid
    ! Mapping from the land to the glc grid is then done with the time averaged fields

    ! uses
    use ESMF             , only : ESMF_GridComp, ESMF_GridCompGet 
    use ESMF             , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF             , only : ESMF_Clock, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet, ESMF_ClockPrint
    use ESMF             , only : ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)    :: clock
    type(ESMF_Time)     :: time
    character(len=64)   :: timestr
    type(InternalState) :: is_local
    integer             :: i,j,n,n1,ncnt
    integer             :: dbrc
    character(len=*),parameter  :: subname='(med_phases_prep_glc_accum_fast)'
    !---------------------------------------

    call t_startf('MED:'//subname)
<<<<<<< HEAD
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
=======

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
>>>>>>> changes to get TG compsets at least running - validation must still be done
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Initialize module field bundles if necessary
    !---------------------------------------

    if (.not. allocated(FBLndAccum)) then
       call med_phases_prep_glc_accum_init(gcomp, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !---------------------------------------
    !--- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt == 0) then
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compglc), returning", &
               ESMF_LOGMSG_INFO, rc=dbrc)
       endif
<<<<<<< HEAD
    else

       !---------------------------------------
       !--- Get the current time from the clock
       !---------------------------------------

       call ESMF_GridCompGet(gcomp, clock=clock)
=======
       RETURN
    end if

    !---------------------------------------
    ! accumulator land input to glc component (on land grid)
    !---------------------------------------

    do n = 1,fieldCount
       call shr_nuopc_methods_FB_accum(FBLndAccum(n), is_local%wrap%FBImp(complnd,complnd), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_diagnose(FBLndAccum(n), string=trim(subname)//&
            ' FBLndAccum for field ' // trim(fieldNameList(n)), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       FBLndAccumCnt(n) = FBLndAccumCnt(n) + 1
    end do

    !---------------------------------------
    ! update local scalar data - set valid input flag to .false.
    !---------------------------------------

    call med_infodata_set_valid_glc_input(.false., med_infodata, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! write diagnostic output
    !---------------------------------------

    if (dbug_flag > 1) then
       call ESMF_GridCompGet(gcomp, clock=clock)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGet(clock,currtime=time,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(time,timestring=timestr)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)

       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_accum_fast

  !================================================================================================

  subroutine med_phases_prep_glc_avg(gcomp, rc)

    ! Prepares the GLC export Fields from the mediator

    use NUOPC                 , only : NUOPC_IsConnected
    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time, ESMF_Array, ESMF_ArrayGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_ClockPrint
    use ESMF                  , only : ESMF_FieldBundleGet, ESMF_FieldBundle, ESMF_FieldBundleIsCreated 
    use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_Mesh, ESMF_MeshGet
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use glc_elevclass_mod     , only : glc_elevclass_as_string

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_FieldBundle) :: FBlnd_topo_l
    type(ESMF_FieldBundle) :: FBlnd_topo_g
    type(ESMF_FieldBundle) :: FBlnd_l
    type(ESMF_FieldBundle) :: FBlnd_g
    real(r8), pointer      :: dataptr_l_in(:)    ! pointer to data coming from FBImp(complnd,complnd) for a given elevation class
    real(r8), pointer      :: dataptr_l(:)       ! temporary data pointer for one elevation class
    real(r8), pointer      :: dataptr_g(:)       ! temporary data pointer for one elevation class
    real(r8), pointer      :: topolnd_l(:)       ! temporary data pointer
    real    , pointer      :: topolnd_g_EC(:,:)  ! topo in elevation classes
    real(r8), pointer      :: topoglc_g(:)       ! ice topographic height on the glc grid extracted from glc import
    real(r8), pointer      :: data_g_bareland(:) 
    real    , pointer      :: data_g_EC(:,:)     ! remapped field in each glc cell, in each EC
    real(r8), pointer      :: data_g_ice_covered(:) ! data for ice-covered regions on the GLC grid
    type(ESMF_Field)       :: lfield
    type(ESMF_Mesh)        :: lmesh
    type(ESMF_Clock)       :: clock
    type(ESMF_Time)        :: time
    character(len=64)      :: timestr
    type(InternalState)    :: is_local
    integer                :: nfld 
    integer                :: ec
    integer                :: i,j,n,g,ncnt
    integer                :: lsize_g
    type(ESMF_Array)       :: elemAreaArray_g
    type(ESMF_Array)       :: elemAreaArray_l
    real(r8), pointer      :: aream_g(:)         ! cell areas on glc grid, for mapping
    real(r8), pointer      :: aream_l(:)         ! cell areas on glc grid, for mapping
    real(r8), pointer      :: area_g(:)          ! cell areas on glc grid, according to glc model
    real(r8), pointer      :: glc_ice_covered(:) ! if points on the glc grid is ice-covered (1) or ice-free (0)
    integer , pointer      :: glc_elevclass(:)   ! elevation classes glc grid
    real(r8), pointer      :: topoptr_g(:)       ! mapped topo from land for a given elevation class
    real(r8), pointer      :: dataexp_g(:)       ! pointer into 
    integer                :: dbrc
    real(r8)               :: elev_l, elev_u     ! lower and upper elevations in interpolation range
    real(r8)               :: d_elev             ! elev_u - elev_l
    character(len=CS), allocatable :: fieldname_EC(:)
    character(len=CS), allocatable :: toponame_EC(:)
    character(len=:) , allocatable :: elevclass_as_string
    character(len=*) , parameter   :: subname='(med_phases_prep_glc_avg)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Initialize module field bundles if necessary
    !---------------------------------------

    if (.not. allocated(FBLndAccum)) then
       call med_phases_prep_glc_accum_init(gcomp, rc=rc)
>>>>>>> changes to get TG compsets at least running - validation must still be done
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

<<<<<<< HEAD
       call ESMF_ClockGet(clock,currtime=time,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet(time,timestring=timestr)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
       endif

       if (mastertask) then
          call ESMF_ClockPrint(clock, options="currTime", preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- mapping
       !---------------------------------------

       do n1 = 1,ncomps
          if (is_local%wrap%med_coupling_active(n1,compglc)) then
             call med_map_FB_Regrid_Norm( &
                  fldListFr(n1)%flds, n1, compglc, &
                  is_local%wrap%FBImp(n1,n1), &
                  is_local%wrap%FBImp(n1,compglc), &
                  is_local%wrap%FBFrac(n1), &
                  is_local%wrap%FBNormOne(n1,compglc,:), &
                  is_local%wrap%RH(n1,compglc,:), &
                  string=trim(compname(n1))//'2'//trim(compname(compglc)), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       enddo

       !---------------------------------------
       !--- auto merges
       !---------------------------------------

       call med_merge_auto(trim(compname(compglc)), &
            is_local%wrap%FBExp(compglc), is_local%wrap%FBFrac(compglc), &
            is_local%wrap%FBImp(:,compglc), fldListTo(compglc), &
            document=first_call, string='(merge_to_lnd)', mastertask=mastertask, rc=rc)
=======
    !---------------------------------------
    ! Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ncnt == 0) then
       call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBExp(compglc), returning", &
            ESMF_LOGMSG_INFO, rc=dbrc)
       RETURN
    end if

    !---------------------------------------
    ! Average import from land accumuled FB
    !---------------------------------------

    do nfld = 1,fieldCount
       call shr_nuopc_methods_FB_average(FBlndAccum(nfld), FBlndAccumCnt(nfld), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_diagnose(FBlndAccum(nfld), string=trim(subname) // &
            ' FBlndAccum for after avg for field bundle '//trim(fieldNameList(nfld)), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! ------------------------------------------------------------------------
    ! Determine elevation class of each glc point on glc grid (output is topoglc_g)
    ! ------------------------------------------------------------------------

    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), Sg_frac_field, glc_ice_covered, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), Sg_topo_field, topoglc_g, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    lsize_g = size(glc_ice_covered)
    allocate(glc_elevclass(lsize_g))
    call get_glc_elevation_classes(glc_ice_covered, topoglc_g, glc_elevclass)

    ! ------------------------------------------------------------------------
    ! Map the topo field from the land grid (in multiple elevation classes)
    ! to the glc grid (in multiple elevation classes) (output is topolnd_g_EC)
    ! ------------------------------------------------------------------------
    
    allocate(toponame_ec(0:nEC))

    do ec = 0,nEC
       elevclass_as_string = glc_elevclass_as_string(ec)
       toponame_EC(ec)  = 'Sl_topo' // elevclass_as_string
    end do

    ! initialize FB on lnd grid for the field coming from land (with elevation classes)
    call shr_nuopc_methods_FB_init(FBlnd_topo_l, flds_scalar_name, &
         FBgeom=is_local%wrap%FBImp(complnd,complnd), fieldNameList=toponame_EC, rc=rc)
    if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

    ! initialize FB on glc grid for the field coming from land (with elevation classes)
    call shr_nuopc_methods_FB_init(FBlnd_topo_g, flds_scalar_name, &
         FBgeom=is_local%wrap%FBImp(complnd,compglc), fieldNameList=toponame_EC, rc=rc)
    if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

    ! Set the data in FBlnd_topo to be just the topo fields coming over from the land import
    do ec = 1,nEC
       call shr_nuopc_methods_FB_getFldPtr(FBlnd_topo_l, toponame_EC(ec), dataptr_l, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(complnd,complnd), toponame_EC(ec), dataptr_l_in, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr_l(:) = dataptr_l_in(:)
    end do
         
    ! map from the land grid to the glc grid to obtain FBlnd_topo_g
    call med_map_FB_Regrid_Norm(FBlnd_topo_l, FBlnd_topo_g, toponame_EC, &
         is_local%wrap%FBFrac(complnd), 'lfrac', &
         is_local%wrap%RH(complnd,compglc,mapbilnr), &
         string='mapping normalized elevation class 0 (bare land) from lnd to to glc', rc=rc)

    ! save the result in topolnd_g_EC(:)
    allocate(topolnd_g_EC (lsize_g, nEC))
    do ec = 1, nEC
       call shr_nuopc_methods_FB_getFldPtr(FBlnd_topo_g, toponame_EC(ec), dataptr_g, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,lsize_g
          topolnd_g_EC(n,ec) = real(dataptr_g(n))
       end do
    enddo

    ! Clean field bundles
    call shr_nuopc_methods_FB_clean(FBlnd_topo_l, rc=rc)
    if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
    call shr_nuopc_methods_FB_clean(FBlnd_topo_g, rc=rc)
    if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

    deallocate(toponame_EC)

    ! ------------------------------------------------------------------------
    ! loop over field in export field bundle to GLC
    ! ------------------------------------------------------------------------

    allocate(fieldname_EC(0:nEC))

    do nfld = 1, fieldCount 

       ! ------------------------------------------------------------------------
       ! Set a pointer to the data for the field that will be sent to glc (without elevation classes)
       ! ------------------------------------------------------------------------

       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBExp(compglc), fieldNameList(nfld), dataexp_g, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! ------------------------------------------------------------------------
       ! Initialize temporary field bundles - and set data in FBlnd_l
       ! ------------------------------------------------------------------------

       do ec = 0, nEC
          elevclass_as_string = glc_elevclass_as_string(ec)
          fieldname_EC(ec) = trim(fieldnameList(nfld)) // elevclass_as_string
       end do

       ! Initialize FB on lnd grid for the field coming from land (with elevation classes)
       call shr_nuopc_methods_FB_init(FBlnd_l, flds_scalar_name, &
            FBgeom=is_local%wrap%FBImp(complnd,complnd), fieldNameList=fieldname_EC, rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
       
       ! Initialize FB on glc grid for the field coming from land (with elevation classes)
       call shr_nuopc_methods_FB_init(FBlnd_g, flds_scalar_name, &
            FBgeom=is_local%wrap%FBImp(complnd,compglc), fieldNameList=fieldname_EC, rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

       ! Set the data in FBlnd_l to be the fields coming over from the land import
       do ec = 0,nEC
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(complnd,complnd), &
               fieldname_EC(ec), dataptr_l_in, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_getFldPtr(FBlnd_l, fieldname_EC(ec), dataptr_l, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          dataptr_l(:) = dataptr_l_in(:)
       end do

       ! ------------------------------------------------------------------------
       ! Map elevation class 0 (bare land) from the land to glc grid and 
       ! set the output data equal to the bare land value everywhere 
       ! This will later get overwritten in places where we have ice
       ! ------------------------------------------------------------------------
       
       ! Map bare land from land to glc 
       call med_map_FB_Regrid_Norm(FBlnd_l, FBlnd_g, (/fieldname_EC(0)/), &
            is_local%wrap%FBFrac(complnd), 'lfrac', &
            is_local%wrap%RH(complnd,compglc,mapbilnr), &
            string='mapping normalized elevation class 0 (bare land) from lnd to to glc', rc=rc)

       call shr_nuopc_methods_FB_getFldPtr(FBlnd_g, fieldname_EC(0), data_g_bareland, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       
       ! TODO(wjs, 2015-01-20) This implies that we pass data to CISM even in places that
       ! CISM says is ocean (so CISM will ignore the incoming value). This differs from the
       ! current glint implementation, which sets acab and artm to 0 over ocean (although
       ! notes that this could lead to a loss of conservation). Figure out how to handle
       ! this case.
       
       allocate(data_g_EC (lsize_g, 0:nEC))
       do ec = 0,nEC
          data_g_EC(:,ec) = real(data_g_bareland(:))
       end do

       ! ------------------------------------------------------------------------
       ! Map the SMB field from the land grid (in multiple elevation classes)
       ! to the glc grid (in multiple elevation classes)
       ! ------------------------------------------------------------------------
       
       ! Regrid the above fields (corresponding to the above fldnames) from the land to the glc grid
       ! using bilinear interpolation (are regridding the SMB field and the topo values for each EC
       ! from the land grid to the glc grid)
       ! The fields in FBlnd_l and in FBlnd_g are in multiple elevation classes 
       
       call med_map_FB_Regrid_Norm(FBlnd_l, FBlnd_g, fieldname_EC(1:nEC), &
            is_local%wrap%FBFrac(complnd), 'lfrac', &
            is_local%wrap%RH(complnd,compglc,mapbilnr), &
            string='mapping normalized elevation class 0 (bare land) from lnd to to glc', rc=rc)

       do ec = 1, nEC
          call shr_nuopc_methods_FB_getFldPtr(FBlnd_g, fieldname_EC(ec), dataptr_g, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          data_g_EC(:,ec) = real(dataptr_g)
       enddo

       ! Perform vertical interpolation of data onto ice sheet topography
       ! This maps all of the input elevation classes into an export to glc without elevation classes

       allocate(data_g_ice_covered(lsize_g))
       data_g_ice_covered(:) = 0._r8
       do n = 1, lsize_g

          ! For each ice sheet point, find bounding EC values...
          if (topoglc_g(n) < topolnd_g_EC(n,1)) then
             ! lower than lowest mean EC elevation value
             data_g_ice_covered(n) = data_g_EC(n,1)

          else if (topoglc_g(n) >= topolnd_g_EC(n,nEC)) then
             ! higher than highest mean EC elevation value
             data_g_ice_covered(n) = data_g_EC(n,nEC)

          else
             ! do linear interpolation of data in the vertical
             do ec = 2, nEC
                if (topoglc_g(n) < topolnd_g_EC(n, ec)) then
                   elev_l = topolnd_g_EC(n, ec-1)
                   elev_u = topolnd_g_EC(n, ec)
                   d_elev = elev_u - elev_l
                   if (d_elev <= 0) then
                      ! This shouldn't happen, but handle it in case it does. In this case,
                      ! let's arbitrarily use the mean of the two elevation classes, rather
                      ! than the weighted mean.
                      write(logunit,*) subname//' WARNING: topo diff between elevation classes <= 0'
                      write(logunit,*) 'n, ec, elev_l, elev_u = ', n, ec, elev_l, elev_u
                      write(logunit,*) 'Simply using mean of the two elevation classes,'
                      write(logunit,*) 'rather than the weighted mean.'
                      data_g_ice_covered(n) = data_g_EC(n,ec-1) * 0.5_r8 &
                                            + data_g_EC(n,ec)   * 0.5_r8
                   else

                      data_g_ice_covered(n) =  data_g_EC(n,ec-1) * (elev_u - topoglc_g(n)) / d_elev  &
                                             + data_g_EC(n,ec)   * (topoglc_g(n) - elev_l) / d_elev
                   end if

                   exit
                end if
             end do
          end if  ! topoglc_g(n)
       end do  ! lsize_g

       deallocate(data_g_EC)

       ! Clean field bundles
       call shr_nuopc_methods_FB_clean(FBlnd_l, rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
       call shr_nuopc_methods_FB_clean(FBlnd_g, rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

       ! ------------------------------------------------------------------------
       ! Overwrite the output data for the ice-covered cells and clean up
       ! NOTE: this updates FBexp data
       ! ------------------------------------------------------------------------

       where (glc_elevclass /= 0)
          dataexp_g = data_g_ice_covered
       end where

       deallocate(data_g_ice_covered)

       ! ------------------------------------------------------------------------
       ! If the field is the flux field from the land giving the surface mass balance to glc
       ! ------------------------------------------------------------------------

       if (trim(fieldNameList(nfld)) == trim(qice_fieldname)) then

          ! Make a preemptive adjustment to qice_g to account for area differences between CISM and the coupler.
          ! - When sending back fluxes to glc, need to multiple the fluxes by by aream_g/area_g for 
          !   conservation purposes.  Where CISM areas are larger (area_g > aream_g), the fluxes are reduced, 
          !   and where CISM areas are smaller, the fluxes are increased.
          ! - As a result, an SMB of 1 m/yr in CLM would be converted to an SMB ranging from ~0.9 to 1.05 m/yr 
          !   in CISM (with smaller values where CISM areas are larger, and larger values where CISM areas are smaller).
          ! - Here, to keep CISM values close to the CLM values in the corresponding locations,
          !   we anticipate the later correction and multiply qice_g by area_g/aream_g.
          !   In the CISM cap, the multiplication by aream_g/area_g will bring qice back to the original values 
          !   obtained from bilinear remapping.
          !
          ! Note that we are free to do this or any other adjustments we want to qice at this
          ! point in the remapping, because the conservation correction will ensure that we
          ! still conserve globally despite these adjustments (and smb_renormalize = .false.
          ! should only be used in cases where conservation doesn't matter anyway).
          
          ! Determine aream_g
          call ESMF_FieldBundleGet(is_local%wrap%FBExp(compglc), qice_fieldname, field=lfield, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          ! TODO (mvertens, 2018-11-30) Cannot query the element area here - unless you create the grid with the addUserArea
          ! The following needs to be commented as a result and needs to be fixed in ESMF

          ! call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
          ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          ! call ESMF_MeshGet(lmesh, elemAreaArray=elemAreaArray_g, rc=rc)
          ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          ! call ESMF_ArrayGet(elemAreaArray_g, farrayptr=aream_g, rc=rc)
          ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          ! Determine areag_g
          ! TODO (mvertens, 2018-11-21) : need to determine area_g - since this is not available no
          ! actually should do this scaling in the component cap when this gets updated - scale this
          ! fluxes in the component cap once it has been received
          ! lsize_g = size(aream_g)
          ! allocate(area_g(lsize_g))
          ! area_g(:) = aream_g(:)
          ! Multiply qice_g (i.e. dataexp_g) by the ratio of area_g/area_m
          ! IMPORTANT NOTE: dataexp_g here is simply a pointer to FBExp for the target fieldname and if
          ! the field name is qice_field name than this is identical to qice_g below
          ! do n = 1, lsize_g
          !    if (aream_g(n) > 0.0_r8) then
          !       dataexp_g(n) = dataexp_g(n) * area_g(n)/aream_g(n)
          !    else
          !       dataexp_g(n) = 0.0_r8
          !    endif
          ! enddo
          ! deallocate(area_g)
          
          ! Renormalizes surface mass balance (smb, here named dataexp_g) so that the global
          ! integral on the glc grid is equal to the global integral on the land grid.

          if (smb_renormalize) then
             call prep_glc_renormalize_smb(gcomp, dataexp_g, rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

       end if ! end of if the field is qice_fieldname

    end do  ! end of loop over fields

    deallocate(fieldname_EC)
    deallocate(topolnd_g_EC)

    !---------------------------------------
    ! NOTE: all merges have already been done - there are no auto merges
    !---------------------------------------

    if (dbug_flag > 1) then
       call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compglc), string=trim(subname)//' FBexp(compglc) ', rc=rc)
>>>>>>> changes to get TG compsets at least running - validation must still be done
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

<<<<<<< HEAD
       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compglc), string=trim(subname)//' FBexp(compglc) ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       !---------------------------------------
       !--- custom calculations
       !---------------------------------------

       !---------------------------------------
       !--- update local scalar data
       !---------------------------------------

       !is_local%wrap%scalar_data(1) =

       !---------------------------------------
       !--- clean up
       !---------------------------------------

       first_call = .false.
    endif
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    call t_stopf('MED:'//subname)
=======
    !---------------------------------------
    ! update local scalar data - set valid input flag to .true.
    !---------------------------------------

    call med_infodata_set_valid_glc_input(.true., med_infodata, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- diagnostic output
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_GridCompGet(gcomp, clock=clock)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGet(clock,currtime=time,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(time,timestring=timestr)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_glc_avg

  !================================================================================================

  subroutine prep_glc_renormalize_smb(gcomp, qice_g, rc)

    !------------------
    ! Renormalizes surface mass balance (smb, here named qice_g) so that the global
    ! integral on the glc grid is equal to the global integral on the land grid.
    !
    ! (1) Map Sg_icemask from the glc grid to the land grid.
    !     Because of coupler lags, the current Sg_icemask_l might not be up to date with Sg_icemask_g.
    ! (2) Map Sg_ice_covered from the glc grid to the land grid
    !     This gives the fields Sg_ice_covered00, Sg_ice_covered01, etc. on the land grid.
    !
    ! This is required for conservation - although conservation is only necessary if we
    ! are running with a fully-interactive, two-way-coupled glc.
    !
    ! For high-level design, see:
    ! https://docs.google.com/document/d/1H_SuK6SfCv1x6dK91q80dFInPbLYcOkUj_iAa6WRnqQ/edit
    !------------------

    ! uses
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet 
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field, ESMF_FieldGet
    use ESMF                  , only : ESMF_Mesh, ESMF_MeshGet, ESMF_Array, ESMF_ArrayGet 
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use med_map_glc2lnd_mod   , only : med_map_glc2lnd_elevclass
    use shr_mpi_mod           , only : shr_mpi_sum, shr_mpi_bcast

    ! input/output variables
    type(ESMF_GridComp)      :: gcomp
    real(r8) , pointer       :: qice_g(:)  ! qice data on glc grid
    integer  , intent(out)   :: rc         ! return error code

    ! local variables
    ! Note: Sg_icemask defines where the ice sheet model can receive a nonzero SMB from the land model.
    type(InternalState)    :: is_local
    type(ESMF_FieldBundle) :: FBlnd_frac_ec
    type(ESMF_FieldBundle) :: FBlnd_icemask
    type(ESMF_Mesh)        :: lmesh
    type(ESMF_Field)       :: lfield
    type(ESMF_Array)       :: elemAreaArray_l
    type(ESMF_Array)       :: elemAreaArray_g
    logical                :: isPresent 
    integer                :: nfld
    real(r8), pointer      :: aream_l(:)      ! cell areas on land grid, for mapping
    real(r8), pointer      :: aream_g(:)      ! cell areas on glc grid, for mapping
    integer                :: lsize_l         ! number of points on land grid
    integer                :: lsize_g         ! number of points on glc grid
    real(r8), pointer      :: qice_l(:,:)     ! SMB (Flgl_qice) on land grid
    real(r8), pointer      :: frac_l(:,:)     ! EC fractions (Sg_ice_covered) on land grid
    real(r8), pointer      :: tmp_field_l(:)  ! temporary field on land grid
    real(r8), pointer      :: Sg_icemask_g(:) ! icemask on glc grid
    real(r8), pointer      :: Sg_icemask_l(:) ! icemask on land grid
    real(r8), pointer      :: lfrac(:)        ! land fraction on land grid
    integer                :: ec              ! loop index over elevation classes 
    integer                :: n
    integer                :: dbrc

    ! various strings for building field names
    character(len=:) , allocatable :: elevclass_as_string
    character(len=CS), allocatable :: qice_fields_ec(:)
    character(len=CS), allocatable :: frac_fields_ec(:)

    ! local and global sums of accumulation and ablation; used to compute renormalization factors
    real(r8) :: local_accum_on_land_grid
    real(r8) :: global_accum_on_land_grid
    real(r8) :: local_accum_on_glc_grid
    real(r8) :: global_accum_on_glc_grid

    real(r8) :: local_ablat_on_land_grid
    real(r8) :: global_ablat_on_land_grid
    real(r8) :: local_ablat_on_glc_grid
    real(r8) :: global_ablat_on_glc_grid

    ! renormalization factors (should be close to 1, e.g. in range 0.95 to 1.05)
    real(r8) :: accum_renorm_factor ! ratio between global accumulation on the two grids
    real(r8) :: ablat_renorm_factor ! ratio between global ablation on the two grids

    real(r8) :: effective_area      ! grid cell area multiplied by min(lfrac,Sg_icemask_l).
                                    ! This is the area that can contribute SMB to the ice sheet model.
    character(len=*), parameter  :: subname='(prep_glc_renormalize_smb)'
    !---------------------------------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Get required pointers
    !---------------------------------------

    ! determine fraction on land grid
    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBFrac(complnd), 'lfrac', lfrac, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine area on land grid (in the mediator)
    call shr_nuopc_methods_FB_getFieldN(is_local%wrap%FBImp(complnd,complnd), fieldnum=1, field=lfield, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh, elemAreaArray=elemAreaArray_l, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ArrayGet(elemAreaArray_l, farrayptr=aream_l, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Map Sg_icemask from the glc grid to the land grid.
    !---------------------------------------

    ! This may not be necessary, if Sg_icemask_l has already been mapped from Sg_icemask_g.
    ! It is done here for two reasons:
    ! (1) The mapping will *not* have been done if we are running with dlnd (e.g., a TG case).
    ! (2) Because of coupler lags, the current Sg_icemask_l might not be up to date with
    !     Sg_icemask_g. This probably isn't a problem in practice, but doing the mapping
    !     here ensures the mask is up to date.
    !
    ! This mapping uses the same options as the standard glc -> lnd mapping done in
    ! prep_lnd_calc_g2x_lx. If that mapping ever changed (e.g., changing norm to
    ! .false.), then we should change this mapping, too.
    !
    ! BUG(wjs, 2017-05-11, #1516) I think we actually want norm = .false. here, but this
    ! requires some more thought

    call shr_nuopc_methods_FB_init(FBlnd_icemask, flds_scalar_name, &
         FBgeom=is_local%wrap%FBImp(compglc,complnd), fieldNameList=(/Sg_icemask_field/), rc=rc)
    if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

    call med_map_FB_Regrid_Norm( &
         is_local%wrap%FBImp(compglc,compglc), FBlnd_icemask, (/Sg_icemask_field/), &
         is_local%wrap%FBNormOne(compglc,compglc,mapconsf ), 'one', &
         is_local%wrap%RH(compglc, complnd, mapconsf), &
         string='mapping Sg_imask_g to Sg_imask_l (from glc to land)', rc=rc)

    ! Note that Sg_icemask_l will not have elevation classes
    call shr_nuopc_methods_FB_getFldPtr(FBlnd_icemask, Sg_icemask_field, Sg_icemask_l)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Map Sg_ice_covered from the glc grid (no elevation classes) to the land grid (with elevation classes) 
    !---------------------------------------

    ! This gives the fields Sg_ice_covered00, Sg_ice_covered01, etc. on the land grid.
    ! These fields are needed to integrate the total SMB on the land grid, for conservation purposes.
    ! As above, the mapping may not be necessary, because Sg_ice_covered might already have been mapped.
    ! However, the mapping will not have been done in a TG case with dlnd, and it might not
    ! be up to date because of coupler lags (though the latter probably isn't a problem
    ! in practice).
    !
    ! Note that, for a case with full two-way coupling, we will only conserve if the
    ! actual land cover used over the course of the year matches these currently-remapped
    ! values. This should generally be the case with the current coupling setup.
    !
    ! One could argue that it would be safer (for conservation purposes) if LND sent its
    ! grid cell average SMB values, or if it sent its own notion of the area in each
    ! elevation class for the purpose of creating grid cell average SMB values here. But
    ! these options cause problems if we're not doing full two-way coupling (e.g., in a TG
    ! case with dlnd, or in the common case where GLC is a diagnostic component that
    ! doesn't cause updates in the glacier areas in LND). In these cases without full
    ! two-way coupling, if we use the LND's notion of the area in each elevation class,
    ! then the conservation corrections would end up correcting for discrepancies in
    ! elevation class areas between LND and GLC, rather than just correcting for
    ! discrepancies arising from the remapping of SMB. (And before you get worried: It
    ! doesn't matter that we are not conserving in these cases without full two-way
    ! coupling, because GLC isn't connected with the rest of the system in terms of energy
    ! and mass in these cases. So in these cases, it's okay that the LND integral computed
    ! here differs from the integral that LND itself would compute.)

    ! Map frac_field on glc grid to frac_fieldNN on land grid (NN are the elevation classes)

    allocate(frac_fields_ec(nEC))
    allocate(qice_fields_ec(nEC))
    do nec = 0, nEC
       elevclass_as_string = glc_elevclass_as_string(ec)
       frac_fields_ec(ec) = Sg_frac_field  // elevclass_as_string   ! Sg_ice_covered01, etc.
       qice_fields_ec(ec) = qice_fieldname // elevclass_as_string   ! Flgl_qice01, etc.
    enddo

    call shr_nuopc_methods_FB_init(FBlnd_frac_ec, flds_scalar_name, &
         FBgeom=is_local%wrap%FBImp(complnd,complnd), fieldNameList=frac_fields_ec, rc=rc)
    if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

    call med_map_glc2lnd_elevclass( FBglc=is_local%wrap%FBImp(compglc, compglc), &
         frac_field=Sg_frac_field, topo_field=Sg_topo_field, icemask_field=Sg_icemask_field, &
         RouteHandle=is_local%wrap%RH(compglc,complnd,mapconsf), FBlnd_ec=FBlnd_frac_ec, rc=rc)
    if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

    do ec = 0, nEC
       call shr_nuopc_methods_FB_getFldPtr(FBlnd_frac_ec, frac_fields_ec(ec), tmp_field_l, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       frac_l(:,ec) = tmp_field_l(:)

       do nfld = 1,fieldCount
          call ESMF_FieldBundleGet(FBLndAccum(nfld), qice_fields_ec(ec), isPresent=isPresent, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          if (isPresent) then
             call shr_nuopc_methods_FB_getFldPtr(FBLndAccum(nfld), qice_fields_ec(ec), tmp_field_l, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             qice_l(:,ec) = tmp_field_l(:)
          end if
       end do
    enddo

    ! Create local arrays over elevation class for qice_l and frac_l
    ! Note: qice comes from l2gacc_lx; frac comes from g2x_lx.

    lsize_l = size(lfrac)
    allocate(qice_l(lsize_l,0:nEC))
    allocate(frac_l(lsize_l,0:nEC))

    ! Sum qice over all elevation classes for each local land grid cell

    local_accum_on_land_grid = 0.0_r8
    local_ablat_on_land_grid = 0.0_r8
    do n = 1, lsize_l

       ! Calculate effective area for sum -  need the mapped Sg_icemask_l
       effective_area = min(lfrac(n),Sg_icemask_l(n)) * aream_l(n)

       do ec = 0, nEC
          if (qice_l(n,ec) >= 0.0_r8) then
             local_accum_on_land_grid = local_accum_on_land_grid + effective_area * frac_l(n,ec) * qice_l(n,ec)
          else
             local_ablat_on_land_grid = local_ablat_on_land_grid + effective_area * frac_l(n,ec) * qice_l(n,ec)
          endif
       enddo  ! ec
    enddo  ! n

    call shr_mpi_sum(local_accum_on_land_grid, global_accum_on_land_grid, is_local%wrap%mpicom, 'accum_l')
    call shr_mpi_sum(local_ablat_on_land_grid, global_ablat_on_land_grid, is_local%wrap%mpicom, 'ablat_l')
    call shr_mpi_bcast(global_accum_on_land_grid, is_local%wrap%mpicom)
    call shr_mpi_bcast(global_ablat_on_land_grid, is_local%wrap%mpicom)

    ! Determine aream_g
    call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(lmesh, elemAreaArray=elemAreaArray_g, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ArrayGet(elemAreaArray_g, farrayptr=aream_g, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Sum qice_g over local glc grid cells.
    ! Note: This sum uses the coupler areas (aream_g), which differ from the native CISM areas.
    !       But since the original qice_g (from bilinear remapping) has been multiplied by
    !         area_g/aream_g above, this calculation is equivalent to multiplying the original qice_g
    !         by the native CISM areas (area_g).
    !       If Flgl_qice were changed to a state (and not included in seq_flds_x2g_fluxes),
    !         then it would be appropriate to use the native CISM areas in this sum.

    ! Determine Sg_icemask_g from the export FB to glc
    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBexp(compglc), Sg_icemask_field, Sg_icemask_g, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    local_accum_on_glc_grid = 0.0_r8
    local_ablat_on_glc_grid = 0.0_r8
    lsize_g = size(qice_g)
    do n = 1, lsize_g
       if (qice_g(n) >= 0.0_r8) then
          local_accum_on_glc_grid = local_accum_on_glc_grid + Sg_icemask_g(n) * aream_g(n) * qice_g(n)
       else
          local_ablat_on_glc_grid = local_ablat_on_glc_grid + Sg_icemask_g(n) * aream_g(n) * qice_g(n)
       endif
    enddo  ! n

    call shr_mpi_sum(local_accum_on_glc_grid, global_accum_on_glc_grid, is_local%wrap%mpicom, 'accum_g')
    call shr_mpi_sum(local_ablat_on_glc_grid, global_ablat_on_glc_grid, is_local%wrap%mpicom, 'ablat_g')
    call shr_mpi_bcast(global_accum_on_glc_grid, is_local%wrap%mpicom)
    call shr_mpi_bcast(global_ablat_on_glc_grid, is_local%wrap%mpicom)

    ! Renormalize
    if (global_accum_on_glc_grid > 0.0_r8) then
       accum_renorm_factor = global_accum_on_land_grid / global_accum_on_glc_grid
    else
       accum_renorm_factor = 0.0_r8
    endif
>>>>>>> changes to get TG compsets at least running - validation must still be done

    if (global_ablat_on_glc_grid < 0.0_r8) then  ! negative by definition
       ablat_renorm_factor = global_ablat_on_land_grid / global_ablat_on_glc_grid
    else
       ablat_renorm_factor = 0.0_r8
    endif

    if (mastertask) then
       write(logunit,*) 'accum_renorm_factor = ', accum_renorm_factor
       write(logunit,*) 'ablat_renorm_factor = ', ablat_renorm_factor
    endif

    do n = 1, lsize_g
       if (qice_g(n) >= 0.0_r8) then
          qice_g(n) = qice_g(n) * accum_renorm_factor
       else
          qice_g(n) = qice_g(n) * ablat_renorm_factor
       endif
    enddo

    deallocate(qice_l)
    deallocate(frac_l)

    call shr_nuopc_methods_FB_clean(FBlnd_icemask, rc=rc)
    if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

    call shr_nuopc_methods_FB_clean(FBlnd_frac_ec, rc=rc)
    if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

  end subroutine prep_glc_renormalize_smb

  !================================================================================================

  subroutine get_glc_elevation_classes(glc_ice_covered, glc_topo, glc_elevclass)

    !-----------------------------------------------------------------------
    ! Get the elevation class of each point on the glc grid.
    ! For grid cells that are ice-free, the elevation class is set to 0.
    ! All arguments (glc_ice_covered, glc_topo and glc_elevclass) must be the same size.
    !-----------------------------------------------------------------------

    use glc_elevclass_mod, only : glc_get_elevation_class
    use glc_elevclass_mod, only : GLC_ELEVCLASS_ERR_NONE, GLC_ELEVCLASS_ERR_TOO_LOW
    use glc_elevclass_mod, only : GLC_ELEVCLASS_ERR_TOO_HIGH, glc_errcode_to_string

    ! input/output variables
    real(r8), intent(in)  :: glc_ice_covered(:) ! ice-covered (1) vs. ice-free (0)
    real(r8), intent(in)  :: glc_topo(:)        ! ice topographic height
    integer , intent(out) :: glc_elevclass(:)   ! elevation class

    ! local variables
    integer :: npts
    integer :: glc_pt
    integer :: err_code

    ! Tolerance for checking whether ice_covered is 0 or 1
    real(r8), parameter :: ice_covered_tol = 1.e-13

    character(len=*), parameter :: subname = 'get_glc_elevation_classes'
    !-----------------------------------------------------------------------

    npts = size(glc_elevclass)
    SHR_ASSERT_FL((size(glc_ice_covered) == npts), __FILE__, __LINE__)
    SHR_ASSERT_FL((size(glc_topo) == npts), __FILE__, __LINE__)

    do glc_pt = 1, npts
       if (abs(glc_ice_covered(glc_pt) - 1._r8) < ice_covered_tol) then

          ! This is an ice-covered point
          call glc_get_elevation_class(glc_topo(glc_pt), glc_elevclass(glc_pt), err_code)

          if ( err_code == GLC_ELEVCLASS_ERR_NONE .or. &
               err_code == GLC_ELEVCLASS_ERR_TOO_LOW .or. &
               err_code == GLC_ELEVCLASS_ERR_TOO_HIGH) then

             ! These are all acceptable "errors" - it is even okay for these purposes if
             ! the elevation is lower than the lower bound of elevation class 1, or
             ! higher than the upper bound of the top elevation class.

             ! Do nothing
          else
             write(logunit,*) subname, ': ERROR getting elevation class for ', glc_pt
             write(logunit,*) glc_errcode_to_string(err_code)
             call shr_sys_abort(subname//': ERROR getting elevation class')
          end if

       else if (abs(glc_ice_covered(glc_pt) - 0._r8) < ice_covered_tol) then

          ! This is a bare land point (no ice)
          glc_elevclass(glc_pt) = 0

       else

          ! glc_ice_covered is some value other than 0 or 1
          ! The lnd -> glc downscaling code would need to be reworked if we wanted to
          ! handle a continuous fraction between 0 and 1.
          write(logunit,*) subname, ': ERROR: glc_ice_covered must be 0 or 1'
          write(logunit,*) 'glc_pt, glc_ice_covered = ', glc_pt, glc_ice_covered(glc_pt)
          call shr_sys_abort(subname//': ERROR: glc_ice_covered must be 0 or 1')
       end if

    end do

  end subroutine get_glc_elevation_classes

end module med_phases_prep_glc_mod
