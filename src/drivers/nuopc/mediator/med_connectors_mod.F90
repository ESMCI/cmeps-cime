module med_connectors_mod

  !-----------------------------------------------------------------------------
  ! Connector phases
  !-----------------------------------------------------------------------------

  use ESMF                  , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_Failure
  use ESMF                  , only : ESMF_State, ESMF_GridComp, ESMF_GridCompGet
  use esmFlds               , only : compatm, compocn, compice, complnd, comprof, compwav, compglc
  use med_internalstate_mod , only : InternalState
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use med_constants_mod     , only : spval => med_constants_spval
  use med_constants_mod     , only : czero => med_constants_czero
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  character(*)      , parameter :: u_FILE_u = &
       __FILE__

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public med_connectors_prep_med2atm
  public med_connectors_prep_med2ocn
  public med_connectors_prep_med2ice
  public med_connectors_prep_med2lnd
  public med_connectors_prep_med2rof
  public med_connectors_prep_med2wav
  public med_connectors_prep_med2glc

  public med_connectors_post_atm2med
  public med_connectors_post_ocn2med
  public med_connectors_post_ice2med
  public med_connectors_post_lnd2med
  public med_connectors_post_rof2med
  public med_connectors_post_wav2med
  public med_connectors_post_glc2med

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private med_connectors_prep_generic
  private med_connectors_post_generic
  private med_connectors_diagnose

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_connectors_prep_generic(gcomp, comp_index, compname, rc)

    use med_infodata_mod      , only : med_infodata
    use med_infodata_mod      , only : med_infodata_CopyInfodataToState
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_copy
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_reset

    ! input/output variables
    type(ESMF_GridComp)          :: gcomp
    integer, intent(in)          :: comp_index
    character(len=*), intent(in) :: compname
    integer, intent(out)         :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_generic)'
    !---------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(compname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !-------------------------
    ! Get the internal state from Component.
    !-------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------
    ! diagnose export state update scalar data in Exp and Imp State
    !-------------------------
      
    is_local%wrap%conn_prep_cnt(comp_index) = is_local%wrap%conn_prep_cnt(comp_index) + 1

    call shr_nuopc_methods_State_reset(is_local%wrap%NStateExp(comp_index), value=spval, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_copy(is_local%wrap%NStateExp(comp_index), is_local%wrap%FBExp(comp_index), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call med_connectors_diagnose(is_local%wrap%NStateExp(comp_index), is_local%wrap%conn_prep_cnt(comp_index), &
         "med_to_"// trim(compname), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateExp(comp_index), &
         'cpl2' //trim(compname), is_local%wrap%mpicom,rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO (mvertens, 2018-12-01): why are we copying to the import state?
    call med_infodata_CopyInfodataToState(med_infodata,is_local%wrap%NStateImp(comp_index), &
         'cpl2'//trim(compname), is_local%wrap%mpicom,rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//trim(compname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)
    
  end subroutine med_connectors_prep_generic

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_generic(gcomp, comp_index, compname, rc)

    use med_infodata_mod      , only : med_infodata
    use med_infodata_mod      , only : med_infodata_CopyStateToInfodata
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_copy
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    integer          , intent(in)  :: comp_index 
    character(len=*) , intent(in)  :: compname
    integer          , intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_generic)'
    !---------------------------------------------

    rc = ESMF_SUCCESS

    ! Note: for information obtained by the mediator always write out the state
    ! if statewrite_flag is .true.
    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(compname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    !-------------------------
    ! Get the internal state from Component.
    !-------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------
    ! diagnose import state
    ! copy import state scalar data to local datatype
    !-------------------------

    is_local%wrap%conn_post_cnt(comp_index) = is_local%wrap%conn_post_cnt(comp_index) + 1

    ! The following will write out the fields in is_local%wrap%NStateImp(comp_index) to local netcdf files
    ! and set the time index in the files to is_local%wrap%conn_post_cnt(comp_index)
    call med_connectors_diagnose(is_local%wrap%NStateImp(comp_index), is_local%wrap%conn_post_cnt(comp_index), &
         " med_from_"//trim(compname), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! The following copies the scalar data from State into the med_infodata instance on the root pe 
    ! then broadcasts the data to all PETs in the mediator
    call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(comp_index), med_infodata, &
         trim(compname) // '2cpl', is_local%wrap%mpicom, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_reset(is_local%wrap%FBImp(comp_index,comp_index), value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! The following copies is_local%wrap%NStateImp(comp_index) to is_local%wrap%FBImp(comp_index,comp_index)
    call shr_nuopc_methods_FB_copy(is_local%wrap%FBImp(comp_index,comp_index), is_local%wrap%NStateImp(comp_index), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//trim(compname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)
  end subroutine med_connectors_post_generic

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2atm(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2atm)'
    !---------------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, compatm, 'atm', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)
  end subroutine med_connectors_prep_med2atm

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2ocn(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2ocn)'
    !---------------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, compocn, 'ocn', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)
  end subroutine med_connectors_prep_med2ocn

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2ice(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2ice)'
    !---------------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, compice, 'ice', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)
  end subroutine med_connectors_prep_med2ice

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2lnd(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2lnd)'
    !---------------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, complnd, 'lnd', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)
  end subroutine med_connectors_prep_med2lnd

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2rof(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer                       :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2rof)'
    !---------------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, comprof, 'rof', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)
  end subroutine med_connectors_prep_med2rof

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2wav(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer                       :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2wav)'
    !---------------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, compwav, 'wav', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)
  end subroutine med_connectors_prep_med2wav

  !-----------------------------------------------------------------------------

  subroutine med_connectors_prep_med2glc(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(med_connectors_prep_med2glc)'
    !---------------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_prep_generic(gcomp, compglc, 'glc', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)
  end subroutine med_connectors_prep_med2glc

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_atm2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_atm2med)'
    !---------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, compatm, 'atm', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_atm2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_ocn2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_ocn2med)'
    !---------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, compocn, 'ocn', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_ocn2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_ice2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_ice2med)'
    !---------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, compice, 'ice', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_ice2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_lnd2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_lnd2med)'
    !---------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, complnd, 'lnd', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_lnd2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_rof2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_rof2med)'
    !---------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, comprof, 'rof', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_rof2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_wav2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_wav2med)'
    !---------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, compwav, 'wav', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_wav2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_post_glc2med(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_connectors_post_glc2med)'
    !---------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call med_connectors_post_generic(gcomp, compglc, 'glc', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_post_glc2med

  !-----------------------------------------------------------------------------

  subroutine med_connectors_diagnose(State, cntr, string, rc)

    use ESMF                  , only : ESMF_State, ESMF_MAXSTR, ESMF_StateGet
    use NUOPC                 , only : NUOPC_Write
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_diagnose
    use med_constants_mod     , only : statewrite_flag => med_constants_statewrite_flag

    ! input/output variables
    type(ESMF_State), intent(in)    :: State
    integer         , intent(inout) :: cntr
    character(len=*), intent(in)    :: string
    integer         , intent(out)   :: rc

    ! local variables
    integer                        :: fieldCount
    character(ESMF_MAXSTR),pointer :: fieldnamelist(:)
    integer                        :: dbrc
    character(len=*),parameter     :: subname='(med_connectors_diagnose)'
    !---------------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//trim(string)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Obtain the field names in State - allocate memory which will be deallocated at the end
    allocate(fieldnamelist(fieldCount))
    call ESMF_StateGet(State, itemNameList=fieldnamelist, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call shr_nuopc_methods_State_diagnose(State, string=trim(subname)//trim(string), rc=rc)
    endif

    ! Write out the fields in State to netcdf files
    if (cntr > 0 .and. statewrite_flag) then
       call ESMF_LogWrite(trim(subname)//trim(string)//": writing out fields", ESMF_LOGMSG_INFO, rc=dbrc)
       call NUOPC_Write(State, &
            fieldnamelist(1:fieldCount), &
            "field_"//trim(string)//"_", timeslice=cntr, &
            overwrite=.true., relaxedFlag=.true., rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    deallocate(fieldnamelist)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//trim(string)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_connectors_diagnose

  !-----------------------------------------------------------------------------

end module med_connectors_mod
