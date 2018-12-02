module shr_nuopc_scalars_mod

  !----------------------------------------------------------------------------
  ! scalars
  !----------------------------------------------------------------------------

  implicit none
  public

  integer, parameter :: flds_scalar_index_nx                 = 1  ! initialization - set by components
  integer, parameter :: flds_scalar_index_ny                 = 2  ! initialization = set by components

  ! initialization - set by rof
  integer, parameter :: flds_scalar_index_rofice_present     = 3  ! does rof have iceberg coupling on
  integer, parameter :: flds_scalar_index_flood_present      = 4  ! does rof have flooding on

  ! initialization - set by ice
  integer, parameter :: flds_scalar_index_iceberg_prognostic = 5  ! does ice model support icebergs

  ! run time settings
  integer, parameter :: flds_scalar_index_precip_fact        = 6   ! runtime - set by ocn
  integer, parameter :: flds_scalar_index_nextsw_cday        = 7   ! runtime - set by atm
  integer, parameter :: flds_scalar_index_valid_glc_input    = 8   ! runtime - set by mediator
                                                                   ! does glc have is valid accumulated data being sent to it?

  ! total number of scalar data in scalar field in any import/export field bundle
  integer, parameter :: flds_scalar_num = 8

  ! name of scalar field 
  character(len=*) , parameter :: flds_scalar_name = "cpl_scalars"

end module shr_nuopc_scalars_mod
