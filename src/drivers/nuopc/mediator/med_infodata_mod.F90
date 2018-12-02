module med_infodata_mod

  ! !DESCRIPTION: A module to get, put, and store some standard scalar data

  ! !USES:

  use med_constants_mod , only : CL, R8
  use esmFlds           , only : ncomps

  implicit none
  private  ! default private

  ! !PUBLIC TYPES:

  public :: med_infodata_type

  ! !PUBLIC MEMBER FUNCTIONS

  public :: med_infodata_GetData         ! Get values from infodata object
  public :: med_infodata_set_valid_glc_input
  public :: med_infodata_CopyStateToInfodata
  public :: med_infodata_CopyInfodataToState

  ! !PUBLIC DATA MEMBERS:
  public :: med_infodata                  ! instance of infodata datatype

  ! InputInfo derived type
  type med_infodata_type
     private

     ! Set via components at initialization
     integer :: nx(ncomps) = -1              ! global nx
     integer :: ny(ncomps) = -1              ! global ny

     ! Set by rof at initialization
     logical :: rofice_present = .false.     ! does rof have iceberg coupling on
     logical :: rof_prognostic = .false.     ! does rof component need input data
     logical :: flood_present = .false.      ! does rof have flooding on

     ! Set by ice at initialization
     logical :: iceberg_prognostic = .false. ! does the ice model support icebergs

     ! Set by atm during run time
     real(R8) :: nextsw_cday = -1.0_R8 ! calendar of next atm shortwave

     ! Set by ocn during runtime
     real(R8) :: precip_fact =  1.0_R8 ! precip factor

     ! Set by mediator during runtime
     logical  :: glc_valid_input = .true. ! is valid accumulated data being sent to prognostic glc

  end type med_infodata_type

  type (med_infodata_type), target :: med_infodata ! single instance for cpl and all comps

  ! used/reused in module

  character(*),parameter :: u_FILE_u = &
    __FILE__

!===============================================================================
CONTAINS
!===============================================================================

  subroutine med_infodata_set_valid_glc_input(glc_valid_input, infodata, rc)

    use ESMF, only : ESMF_SUCCESS

    logical                , intent(in)    :: glc_valid_input
    type(med_infodata_type), intent(inout) :: infodata
    integer   ,   optional , intent(out)   :: rc 
    
    rc = ESMF_SUCCESS

    infodata%glc_valid_input = glc_valid_input

  end subroutine med_infodata_set_valid_glc_input

  !================================================================================

  subroutine med_infodata_CopyStateToInfodata(State, infodata, type, mpicom, rc)

    use mpi  ! TODO (mvertens 2018-12-01) use EMSF communication routines rather than mpi
    use ESMF                  , only : ESMF_State, ESMF_Field, ESMF_StateItem_Flag
    use ESMF                  , only : ESMF_StateGet, ESMF_FieldGet, ESMF_LogWrite
    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGMSG_INFO
    use ESMF                  , only : ESMF_STATEITEM_NOTFOUND, operator(==)
    use esmFlds               , only : compname
    use shr_nuopc_scalars_mod , only : flds_scalar_num, flds_scalar_name
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nx, flds_scalar_index_ny
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nextsw_cday
    use shr_nuopc_scalars_mod , only : flds_scalar_index_flood_present
    use shr_nuopc_scalars_mod , only : flds_scalar_index_rofice_present
    use shr_nuopc_scalars_mod , only : flds_scalar_index_precip_fact
    use shr_nuopc_scalars_mod , only : flds_scalar_index_valid_glc_input
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkErr

    ! ----------------------------------------------
    ! Copy scalar data from State to local data on root then broadcast data
    ! to all PETs in component.
    ! ----------------------------------------------

    type(ESMF_State),        intent(in)     :: State
    type(med_infodata_type), intent(inout)  :: infodata
    character(len=*),        intent(in)     :: type
    integer,                 intent(in)     :: mpicom
    integer,                 intent(inout)  :: rc

    ! local variables
    integer                         :: n
    integer                         :: mytask, ierr, len
    character(MPI_MAX_ERROR_STRING) :: lstring
    type(ESMF_Field)                :: field
    type(ESMF_StateItem_Flag)       :: itemType
    real(R8), pointer               :: farrayptr(:,:)
    real(R8)                        :: data(flds_scalar_num)
    character(len=32)               :: ntype
    integer                         :: dbrc
    character(len=1024)             :: msgString
    character(len=*), parameter     :: subname='(med_infodata_CopyStateToInfodata)'
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    call MPI_COMM_RANK(mpicom, mytask, rc)
    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), itemType=itemType, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (itemType == ESMF_STATEITEM_NOTFOUND) then
       call ESMF_LogWrite(trim(subname)//": "//trim(flds_scalar_name)//" not found", ESMF_LOGMSG_INFO, &
            line=__LINE__, file=u_FILE_u, rc=dbrc)
    else
      call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (mytask == 0) then
        call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (size(data) < flds_scalar_num .or. size(farrayptr) < flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR on data size", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
        data(1:flds_scalar_num) = farrayptr(1,1:flds_scalar_num)
      endif

      call MPI_BCAST(data, flds_scalar_num, MPI_REAL8, 0, mpicom, rc)
      if (rc /= MPI_SUCCESS) then
        call MPI_ERROR_STRING(rc,lstring,len,ierr)
        call ESMF_LogWrite(trim(subname)//": ERROR "//trim(lstring), ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

      do n = 1,ncomps
         ntype = trim(compname(n))//'2cpli'
         if (trim(type) == trim(ntype)) then
            infodata%nx(n) = nint(data(flds_scalar_index_nx))
            infodata%ny(n) = nint(data(flds_scalar_index_ny))
            write(msgString,'(2i8,2l4)') nint(data(flds_scalar_index_nx)),nint(data(flds_scalar_index_ny))
            call ESMF_LogWrite(trim(subname)//":"//trim(type)//":"//trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
         endif
      enddo

      if (type == 'atm2cpli') then
        infodata%nextsw_cday = data(flds_scalar_index_nextsw_cday)

      elseif (type == 'ocn2cpli') then
        infodata%precip_fact=data(flds_scalar_index_precip_fact)

      elseif (type == 'ice2cpli') then
        ! nothing

      elseif (type == 'lnd2cpli') then
        ! nothing

      elseif (type == 'rof2cpli') then
        infodata%flood_present  = (nint(data(flds_scalar_index_flood_present))  /= 0)
        infodata%rofice_present = (nint(data(flds_scalar_index_rofice_present)) /= 0)

      elseif (type == 'wav2cpli') then
        ! nothing

      elseif (type == 'glc2cpli') then
        ! nothing

      elseif (type == 'atm2cpl') then
         infodata%nextsw_cday=data(flds_scalar_index_nextsw_cday)

      elseif (type == 'ocn2cpl') then
         infodata%precip_fact=data(flds_scalar_index_precip_fact)

      elseif (type == 'ice2cpl') then
        ! nothing

      elseif (type == 'lnd2cpl') then
        ! nothing

      elseif (type == 'rof2cpl') then
        ! nothing

      elseif (type == 'wav2cpl') then
        ! nothing

      elseif (type == 'glc2cpl') then
        ! nothing

      else
         call ESMF_LogWrite(trim(subname)//": ERROR in type = "//trim(type), &
              ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

    endif

  end subroutine med_infodata_CopyStateToInfodata

  !================================================================================

  subroutine med_infodata_CopyInfodataToState(infodata, State, type, mpicom, rc)

    ! ----------------------------------------------
    ! Copy local scalar data into State, root only,
    ! but called on all PETs in mediator
    ! ----------------------------------------------

    use ESMF                  , only : ESMF_State, ESMF_StateGet, ESMF_Field, ESMF_StateItem_Flag, ESMF_FieldGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_STATEITEM_NOTFOUND
    use ESMF                  , only : operator(==), ESMF_FAILURE
    use mpi                   , only : mpi_comm_rank
    use shr_nuopc_scalars_mod , only : flds_scalar_num, flds_scalar_name
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nx, flds_scalar_index_ny
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nextsw_cday
    use shr_nuopc_scalars_mod , only : flds_scalar_index_precip_fact
    use shr_nuopc_scalars_mod , only : flds_scalar_index_valid_glc_input
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkErr

    ! iput/output variables
    type(med_infodata_type), intent(in)    :: infodata
    type(ESMF_State),        intent(inout) :: State
    character(len=*),        intent(in)    :: type
    integer,                 intent(in)    :: mpicom
    integer,                 intent(inout) :: rc

    ! local variables
    integer                   :: mytask
    type(ESMF_Field)          :: field
    type(ESMF_StateItem_Flag) :: ItemType
    real(R8), pointer         :: farrayptr(:,:)
    integer                   :: dbrc
    character(len=*), parameter :: subname='(med_infodata_CopyInfodataToState)'
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    call MPI_COMM_RANK(mpicom, mytask, rc)

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), itemType=itemType, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (itemType == ESMF_STATEITEM_NOTFOUND) then

       call ESMF_LogWrite(trim(subname)//": "//trim(flds_scalar_name)//" not found", &
            ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)

    else

      call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (mytask == 0) then

        call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        if (size(farrayptr) < flds_scalar_num) then
           call ESMF_LogWrite(trim(subname)//": ERROR on data size", &
                ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          RETURN
        endif

        ! Set scalar data into farraptr array
        farrayptr(1,flds_scalar_index_nextsw_cday) = infodata%nextsw_cday

        farrayptr(1,flds_scalar_index_precip_fact) = infodata%precip_fact

        if (infodata%glc_valid_input) then
           farrayptr(1,flds_scalar_index_valid_glc_input) = 1._r8 
        else
           farrayptr(1,flds_scalar_index_valid_glc_input) = 0._r8 
        end if

      endif ! my_task == 0

    endif

  end subroutine med_infodata_CopyInfodataToState

  !===============================================================================

  subroutine med_infodata_GetData( infodata, ncomp, flux_epbal, flux_epbalfact, nx, ny, rc)

    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogWrite
    use med_internalstate_mod , only : logunit, loglevel

    ! Get values out of the infodata object.
  
    ! input/output variables
    type(med_infodata_type) , intent(in)  :: infodata       ! Input CCSM structure
    integer    ,   optional , intent(in)  :: ncomp          ! Component ID
    character(CL), optional , intent(in)  :: flux_epbal     ! selects E,P,R adjustment technique
    real(R8),      optional , intent(out) :: flux_epbalfact ! adjusted precip factor
    integer    ,   optional , intent(out) :: nx             ! nx
    integer    ,   optional , intent(out) :: ny             ! ny
    integer    ,   optional , intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname = '(med_infodata_GetData) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if ( present(flux_epbalfact)) then
       if (.not. present(flux_epbal)) then
          call ESMF_LogWrite(trim(subname)//&
               "Must provide flux_epbal as an input argument to determine infodata%precip_fact")
          rc = ESMF_FAILURE
          RETURN
       end if

       flux_epbalfact = 1.0_R8
       if (trim(flux_epbal) == 'ocn') then
          flux_epbalfact = infodata%precip_fact
          if (flux_epbalfact <= 0.0_R8) then
             if (loglevel > 0) then
                write(logunit,'(2a,e16.6)') trim(subname),' WARNING: factor from ocn = ',flux_epbalfact
                write(logunit,'(2a)') trim(subname),' WARNING: resetting flux_epbalfact to 1.0'
             end if
             flux_epbalfact = 1.0_R8
          end if
       end if
    endif

    if ( (present(nx) .and. .not.present(ncomp)) .or. &
         (present(ny) .and. .not.present(ncomp)) ) then
       call ESMF_LogWrite(trim(subname)//"Must provide ncomp")
       rc = ESMF_FAILURE
       RETURN
    endif
    ny = infodata%ny(ncomp)
    nx = infodata%nx(ncomp)

  end subroutine med_infodata_GetData

end module med_infodata_mod
