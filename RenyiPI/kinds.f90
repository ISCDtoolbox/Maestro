! ***********************************************
!>\file kinds.f90
! ***********************************************


! ***********************************************
!>\brief Double precision kind
! ***********************************************

module kinds

! ***********************************************
!>\brief Double precision kind
! ***********************************************

  integer, parameter :: dp=kind(0.0D0)
  integer, parameter :: dp_size=8

  integer, parameter :: sg=kind(0.0)
  integer, parameter :: sg_size=4

end module kinds
