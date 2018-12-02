module med_map_glc2lnd_mod

  !---------------------------------------------------------------------
  ! This module contains routines for mapping fields from the GLC grid onto the LND grid
  ! (separated by GLC elevation class)
  !
  ! For high-level design, see:
  ! https://docs.google.com/document/d/1sjsaiPYsPJ9A7dVGJIHGg4rVIY2qF5aRXbNzSXVAafU/edit?usp=sharing
  !---------------------------------------------------------------------

#include "shr_assert.h"

  use med_constants_mod, only : CS

  implicit none
  private

  ! interfaces
  public  :: med_map_glc2lnd_elevclass ! map all fields from GLC -> LND grid will be separated by elevation class

  private :: get_glc_elevation_classes ! get elevation class of each glc cell

  ! private module variables
  character(len=CS) , allocatable :: fields_to_map(:)      ! fields_to_map without elev class suffixes
  character(len=CS) , allocatable :: fields_to_map_ec(:,:) ! fields_to_map with elev class suffixes
  character(len=CS) , allocatable :: frac_field_ec(:)      ! name of frac_fields with elev class suffix

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================================
contains
!================================================================================================
  
  subroutine med_map_glc2lnd_elevclass(FBglc, frac_field, topo_field, &
       icemask_field, extra_fields, RouteHandle, FBlnd_ec, rc)

    !------------------
    ! Maps fields from the GLC grid to the LND grid. 
    ! On the GLC grid the fields will not have elevation classes. 
    ! On the LND grid they will have elevation classes.
    !
    ! Maps frac_field, topo_field, plus all fields defined in extra_fields. extra_fields
    ! should be a colon-delimited list of fields, giving the field name in the g2x_g
    ! attribute vector (i.e., without the elevation class suffixes).
    !
    ! Assumes that FBglc contains:
    ! - frac_field
    ! - topo_field
    ! - icemask_field (Note: this is NOT mapped here, but is needed as an input to the mapping)
    ! - each field in extra_fields
    !
    ! Assumes that Fbout (on land grid) contains:
    ! - <frac_field>00, <frac_field>01, <frac_field>02, ...
    ! - <topo_field>00, <topo_field>01, <topo_field>02, ...
    ! - And similarly for each field in extra_fields
    !
    ! Currently assumes that all fields are mapped using the same mapper, which should be
    ! a conservative mapper (i.e., a flux mapper).
    !------------------

    ! uses
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_RouteHandle
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_clean
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_scalars_mod , only : flds_scalar_name
    use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
    use med_constants_mod     , only : R8
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use glc_elevclass_mod     , only : glc_mean_elevation_virtual
    use glc_elevclass_mod     , only : glc_get_num_elevation_classes, glc_get_elevation_class
    use glc_elevclass_mod     , only : glc_elevclass_as_string
    use glc_elevclass_mod     , only : glc_get_num_elevation_classes
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)           :: FBglc
    character(len=*)       , intent(in)           :: frac_field    ! name of field in g2x_g containing glc ice fraction
    character(len=*)       , intent(in)           :: topo_field    ! name of field in g2x_g containing glc topo
    character(len=*)       , intent(in)           :: icemask_field ! name of field in g2x_g containing ice mask
    character(len=*)       , intent(in), optional :: extra_fields(:)
    type(ESMF_RouteHandle) , intent(inout)        :: RouteHandle
    type(ESMF_FieldBundle) , intent(inout)        :: FBlnd_ec
    integer                , intent(out)          :: rc

    ! local variables
    character(len=:), allocatable :: elevclass_as_string
    integer                :: nec, nfld
    integer                :: numEC                               ! number of elevation classes
    integer                :: lsize_g
    integer                :: num_elevation_classes
    integer                :: num_extra_fields
    integer                :: ntopo_field 
    real(r8)               :: topo_virtual
    type(ESMF_FieldBundle) :: FBlnd_other                         ! temporary field bundle
    type(ESMF_FieldBundle) :: FBglc_other                         ! temporary field bundle
    type(ESMF_FieldBundle) :: FBglc_norm                          ! temporary field bundle
    type(ESMF_FieldBundle) :: FBglc_frac_ec                       ! temporary field bundle
    integer  , allocatable :: glc_elevclass(:)                    ! elevation class of each glc cell (assuming cell is ice-covered)
    real(r8) , pointer     :: glc_frac_g(:)                       ! total ice fraction in each glc cell
    real(r8) , pointer     :: glc_topo_g(:)                       ! topographic height of each glc cell
    real(r8) , pointer     :: glc_icemask_g(:)                    ! glc ice mask field on glc grid
    real(r8) , pointer     :: glc_frac_this_ec_g(:)               ! glc fractions in one elev class, on the glc grid
    real(r8) , pointer     :: glc_frac_this_ec_times_icemask_g(:) ! (glc fraction in one elev class) x (icemask), on the glc grid
    real(r8) , pointer     :: glc_frac_this_ec_l(:)               ! glc fractions in one elev class mapped to the land grid
    real(r8) , pointer     :: glc_topo_this_ec_l(:)               ! glc topo field mapped to land grid
    real(r8) , pointer     :: dataptr2(:) 
    real(r8) , pointer     :: dataptr1(:) 
    integer                :: dbrc
    integer                :: fieldcount
    character(len=CS), allocatable :: fieldnameList(:)
    logical                        :: first_call = .true.
    character(len=*), parameter    :: subname = 'map_glc2lnd_ec'
    !-----------------------------------------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS

    ! ------------------------------------------------------------------------
    ! Allocate and set module variables for character arrays
    ! ------------------------------------------------------------------------

    if (first_call) then

       call ESMF_FieldBundleGet(FBglc, fieldcount=fieldcount, rc=rc)
       allocate(fieldnamelist(fieldcount))
       call ESMF_FieldBundleGet(FBglc, fieldNameList=fieldNameList, rc=rc)

       ! TODO (mvertens, 2018-11-25): implement the folowing, the call to 
       ! ESMF_FieldBundleGet is incorret
       ! reauire that frac_field, topo_field and icemask_field be in FBglc
       ! call ESMF_FieldBundleGet(FBglc, trim(frac_field), isPresent=isPresent, rc=rc)
       ! if (.not. isPresent) then
       !    ! TODO: error
       ! end if
       ! call ESMF_FieldBundleGet(FBglc, trim(topo_field), isPresent=isPresent, rc=rc)
       ! if (.not. isPresent) then
       !    ! TODO: error
       ! end if
       ! call ESMF_FieldBundleGet(FBglc, trim(icemask_field), isPresent=isPresent, rc=rc)
       ! if (.not. isPresent) then
       !    ! TODO: error
       ! end if

       numEC = glc_get_num_elevation_classes()

       allocate(frac_field_ec(numEC))
       do nec = 0, numEC
          elevclass_as_string = glc_elevclass_as_string(nec)
          frac_field_ec(nec) = trim(frac_field) // elevclass_as_string    ! Sg_ice_covered01, etc.
       end do

       if (present(extra_fields)) then
          num_extra_fields = size(extra_fields)

          allocate(fields_to_map(num_extra_fields))  
          allocate(fields_to_map_ec(num_extra_fields, num_elevation_classes))

          do nfld = 1,size(extra_fields)
             fields_to_map(nfld) = extra_fields(nfld)
             do nec = 0, num_elevation_classes
                fields_to_map_ec(nfld,nec) = trim(extra_fields(nfld)) // glc_elevclass_as_string(nec)
             end do
             if (trim(fields_to_map(nfld)) == trim(topo_field)) then
                ntopo_field = nfld
             end if
          end do

          ! TODO: check that each extra field is in FBglc 
          ! TODO: check that topo_field is contained in the set of extra fields
       end if

       first_call = .false.
    end if

    ! ------------------------------------------------------------------------
    ! Determine ice mask field on glc grid (glc_icemask)
    ! ------------------------------------------------------------------------

    call shr_nuopc_methods_FB_getFldPtr(FBglc, trim(icemask_field), glc_icemask_g)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! Determine elevation class of each glc point on the glc grid (glc_elevclass)
    ! ------------------------------------------------------------------------

    ! glc_frac is the total ice fraction in each glc gridcell
    call shr_nuopc_methods_FB_getFldPtr(FBglc, trim(frac_field), glc_frac_g, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! glc_topo is the topographic height of each glc gridcell
    call shr_nuopc_methods_FB_getFldPtr(FBglc, trim(topo_field), glc_topo_g,  rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! given glc_topo for each glc grid cell, set the glc elevation classes, glc_elevclass
    lsize_g = size(glc_frac_g)
    allocate(glc_elevclass(lsize_g))
    call get_glc_elevation_classes(glc_topo_g, glc_elevclass)

    ! ------------------------------------------------------------------------
    ! Map fraction in each elevation class to the land grid
    ! ------------------------------------------------------------------------

    ! Loop over all elevation classes
    do nec = 0, num_elevation_classes
      
       !---------------------------------
       ! Determine fractional ice coverage for each elevation class on the land grid
       !---------------------------------

       ! create temporary field bundle on glc grid containing field name frac_field_ec(nec)
       call shr_nuopc_methods_FB_init(FBglc_frac_ec, flds_scalar_name, FBgeom=FBglc, fieldNameList=(/frac_field_ec(nec)/), rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

       ! set pointer to array holding glc fraction in one elev class, on the glc grid
       call shr_nuopc_methods_FB_getFldPtr(FBglc_frac_ec, frac_field_ec(nec), glc_frac_this_ec_g, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! set fractional ice coverage for a given elevation class (THIS IS WHAT WILL GET MAPPED)
       ! this will set the array glc_frac_this_ec_g, ice fractions in this elevation class, for each glc gridcell.
       ! note that glc_elevclass gives the elevation class of each glc grid cell, assuming that 
       ! the grid cell is ice-covered.
       if (nec == 0) then
          glc_frac_this_ec_g(:) = 1._r8 - glc_frac_g(:)
       else
          where (glc_elevclass == nec)
             glc_frac_this_ec_g = glc_frac_g
          elsewhere
             glc_frac_this_ec_g = 0._r8
          end where
       end if

       ! map fraction in each elevation class from the glc grid to the land grid
       ! the result will be contained in glc_frac_this_ec_l(:)
       call med_map_FB_Regrid_Norm( FBSrc=FBglc_frac_ec, FBDst=FBlnd_ec, fldnames=(/frac_field_ec(nec)/), &
            FBNorm=FBglc, fldnorm=trim(icemask_field), RouteHandle=RouteHandle, &
            string='mapping elevation class '//trim(glc_elevclass_as_string(nec))//' fraction from glc to lnd', rc=rc)
       
       ! clean temporary field bundle
       call shr_nuopc_methods_FB_clean(FBglc_frac_ec, rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

       !---------------------------------
       ! Map the other fields if needed
       !---------------------------------

       if (num_extra_fields > 0) then

          ! Create a mask that is (fraction in this elevation class) x (icemask). So, only
          ! grid cells that are both (a) within the icemask and (b) in this elevation class
          ! will be included in the following mapping.
          
          call shr_nuopc_methods_FB_init(FBglc_norm, flds_scalar_name, &
               FBgeom=FBglc, fieldNameList=(/'Sg_frac_times_icemask'/), rc=rc)
          if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

          call shr_nuopc_methods_FB_getFldPtr(FBglc_norm, 'Sg_frac_times_icemask', glc_frac_this_ec_times_icemask_g, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          glc_frac_this_ec_times_icemask_g(:) = glc_frac_this_ec_g(:) * glc_icemask_g(:)

          ! Map other fields to the land grid
          ! Note that bare land values are mapped in the same way as ice-covered values

          call shr_nuopc_methods_FB_init(FBglc_other, flds_scalar_name, &
               FBgeom=FBglc, fieldNameList=fields_to_map, rc=rc)
          if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

          call shr_nuopc_methods_FB_init(FBlnd_other, flds_scalar_name, &
               FBgeom=FBglc, fieldNameList=fields_to_map, rc=rc)
          if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

          ! map all other fields normalizing by glc_frac_this_ec_times_icemask_g - this is what introduces
          ! elevation class information from the glc grid (without elevation classes) to the land grid (with 
          ! elevation classes)

          call med_map_FB_Regrid_Norm(FBSrc=FBglc_other, FBDst=FBlnd_other, fldnames=fields_to_map, &
               FBNorm=FBglc_norm, fldnorm='Sg_frac_times_icemask', RouteHandle=RouteHandle, &
               string='mapping other fields from glc to land elevation class'//trim(elevclass_as_string), rc=rc)

          ! Fill in FBlnd_ec for other fields

          do nfld = 1,size(fields_to_map)
             call shr_nuopc_methods_FB_getFldPtr(FBlnd_other, trim(fields_to_map(nfld)), dataptr1, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             call shr_nuopc_methods_FB_getFldPtr(FBlnd_ec, trim(fields_to_map_ec(nfld,nec)), dataptr2, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             dataptr2(:) = dataptr1(:)
          end do

          call shr_nuopc_methods_FB_clean(FBglc_other, rc=rc)
          if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
          
          call shr_nuopc_methods_FB_clean(FBlnd_other, rc=rc)
          if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

          ! ------------------------------------------------------------------------
          ! set the topo field for virtual columns, in a given elevation class.
          ! ------------------------------------------------------------------------

          ! This is needed because virtual columns (i.e., elevation classes that have no
          ! contributing glc grid cells) won't have any topographic information mapped onto
          ! them, so would otherwise end up with an elevation of 0.

          call shr_nuopc_methods_FB_getFldPtr(FBlnd_ec, trim(frac_field_ec(nec)), glc_frac_this_ec_l, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call shr_nuopc_methods_FB_getFldPtr(FBlnd_ec, trim(fields_to_map_ec(ntopo_field,nec)), glc_topo_this_ec_l, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          topo_virtual = glc_mean_elevation_virtual(nec)
          where (glc_frac_this_ec_l <= 0)
             glc_topo_this_ec_l = topo_virtual
          end where

       end if ! if num_extra_fields > 0

    end do  ! loop over elevation classes

    ! clean up
    deallocate(glc_elevclass)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
       call t_stopf('MED:'//subname)
    end if

  end subroutine med_map_glc2lnd_elevclass

!================================================================================================

  subroutine get_glc_elevation_classes(glc_topo, glc_elevclass)

    !------------------
    ! Get elevation class of each grid cell on the glc grid.
    !
    ! This does not consider glc_frac: it simply gives the elevation class that the grid
    ! cell would be in if it were ice-covered. So it never returns an elevation class of
    ! 0 (bare land). (This design would allow us, in the future, to have glc grid cells
    ! that are part ice-covered, part ice-free.)
    !------------------

    ! uses
    use med_constants_mod     , only : R8
    use med_internalstate_mod , only : logunit
    use glc_elevclass_mod     , only : GLC_ELEVCLASS_ERR_NONE, GLC_ELEVCLASS_ERR_TOO_LOW
    use glc_elevclass_mod     , only : GLC_ELEVCLASS_ERR_TOO_HIGH, glc_errcode_to_string
    use glc_elevclass_mod     , only : glc_get_elevation_class
    use shr_sys_mod           , only : shr_sys_abort

    ! input/output variables
    real(r8), intent(in)  :: glc_topo(:)      ! topographic height
    integer , intent(out) :: glc_elevclass(:) ! elevation class
    !
    ! local variables
    integer :: npts
    integer :: glc_pt
    integer :: err_code
    character(len=*), parameter :: subname = 'get_glc_elevation_classes'
    !-----------------------------------------------------------------------

    npts = size(glc_elevclass)
    SHR_ASSERT_FL((size(glc_topo) == npts), __FILE__, __LINE__)

    do glc_pt = 1, npts
       call glc_get_elevation_class(glc_topo(glc_pt), glc_elevclass(glc_pt), err_code)
       select case (err_code)
       case (GLC_ELEVCLASS_ERR_NONE)
          ! Do nothing
       case (GLC_ELEVCLASS_ERR_TOO_LOW, GLC_ELEVCLASS_ERR_TOO_HIGH)
          write(logunit,*) subname, ': WARNING, for glc_pt, topo = ', glc_pt, glc_topo(glc_pt)
          write(logunit,*) glc_errcode_to_string(err_code)
       case default
          write(logunit,*) subname, ': ERROR getting elevation class for glc_pt = ', glc_pt
          write(logunit,*) glc_errcode_to_string(err_code)
          call shr_sys_abort(subname//': ERROR getting elevation class')
       end select
    end do

  end subroutine get_glc_elevation_classes

end module med_map_glc2lnd_mod
