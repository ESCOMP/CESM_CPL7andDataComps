module seq_drydep_mod

  use shr_drydep_mod, only: seq_drydep_readnl=>shr_drydep_readnl
  use shr_drydep_mod, only: seq_drydep_setHCoeff=>shr_drydep_setHCoeff
  use shr_drydep_mod

  implicit none

  ! method specification
  character(len=*), parameter :: DD_XLND = 'xactive_lnd' ! dry-dep land
  character(len=*), parameter :: drydep_method = DD_XLND ! XLND is the only option now
  logical, parameter :: lnd_drydep = .true.

end module seq_drydep_mod
