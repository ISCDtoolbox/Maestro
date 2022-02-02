subroutine  root2mat(mat,errnoise) 

  real*8 mat(2,2),z,a,b,c,aa,bb,cc
  integer errnoise
!   This subroutine computes the square root of a positive definite
!   symmetric 2x2 matrix A
!   matrix. The output is the unique 2x2 symmetric positive definite matrix B such
!   that B^2 =A. If the input is positive definite the output replaces
!   the input, otherwise the input matrix is unchanged.
  aa=mat(1,1)
  bb=mat(2,1)
  cc=mat(2,2)
  errnoise=0
  if(aa.lt.0.d0 .or. cc.lt.0.d0 .or. aa*cc-bb**2 .lt.0.d0) then
!         check if it is positive definite
     errnoise=2
     return
  endif 

  z=aa+cc+dsqrt(4.d0*(aa*cc-bb**2))

  if(z.le.0.d0) then
!         square root of zero is zero.
     return
  endif
  z=dsqrt(z)

  a=0.5d0*(z+(aa-cc)/z)
  b=bb/z
  c=0.5d0*(z-(aa-cc)/z)
         
  mat(1,1)=a
  mat(2,1)=b
  mat(1,2)=b
  mat(2,2)=c 

  return
end


