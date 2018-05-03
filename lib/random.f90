function rng(aleat)

  use f90_kind

  integer(int32) :: aleat
  real(double) :: rng
  real(double), parameter :: huge=4294967295.d0

  aleat = ior(aleat,1)
  aleat = aleat + ishft(aleat,16) + ishft(aleat,1)

  rng = 0.5D0 + float(aleat)/huge

  return

end function rng
