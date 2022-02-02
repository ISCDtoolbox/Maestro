subroutine scnd(t)

  implicit none
  real(8) :: t
  real(4) :: etime
  real(4), dimension (:), allocatable :: tarray

  allocate(tarray(2))

  t = etime(tarray)
  t = tarray(1)
   
  deallocate(tarray)

  return
end
