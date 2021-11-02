module seq_drydep_mod

  use shr_drydep_mod, only: seq_drydep_setHCoeff=>shr_drydep_setHCoeff
  use shr_drydep_mod

  implicit none

  ! method specification
  character(len=*), parameter :: DD_XLND = 'xactive_lnd' ! dry-dep land
  character(len=*), parameter :: drydep_method = DD_XLND ! XLND is the only option now
  logical, protected :: lnd_drydep

contains

  subroutine seq_drydep_readnl(NLFilename, ID, shr_drydep_fields)

    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    integer         , intent(in)  :: ID         ! seq_comm ID
    character(len=*), intent(out) :: shr_drydep_fields

    call shr_drydep_readnl(NLFilename, ID, shr_drydep_fields)

    lnd_drydep = n_drydep>0

  end subroutine seq_drydep_readnl

end module seq_drydep_mod
