subroutine r3cross(x,y,z)
  implicit none
  double precision, intent(in) :: x(3)
  double precision, intent(in) :: y(3)
  double precision, intent(out) :: z(3)
  z(1)=x(2)*y(3)-x(3)*y(2)
  z(2)=x(3)*y(1)-x(1)*y(3)
  z(3)=x(1)*y(2)-x(2)*y(1)
  return
end subroutine


!!! this subroutine is taken from the paper:    ACS omega 2019, 4, 9271-9283
subroutine eckart(nat,ndim,x_eq,x_loc,amas,u_mat)

  implicit none
  integer i,nat,ndim,lworkd,info
  double precision amas(nat),x_eq(nat,ndim),x_loc(nat,ndim)
  double precision, allocatable :: cmm(:,:),eigv(:),qvec(:),work(:)
  double precision, allocatable :: csi1(:,:), csi2(:,:)
  double precision u_mat(3,3)
  
  allocate(cmm(4,4),eigv(4),qvec(4),csi1(nat,ndim),csi2(nat,ndim))
  
  do i=1,nat
    csi1(i,1)=x_eq(i,1)+x_loc(i,1);  csi2(i,1)=x_eq(i,1)-x_loc(i,1)
    csi1(i,2)=x_eq(i,2)+x_loc(i,2);  csi2(i,2)=x_eq(i,2)-x_loc(i,2)
    csi1(i,3)=x_eq(i,3)+x_loc(i,3);  csi2(i,3)=x_eq(i,3)-x_loc(i,3)
  enddo
  
  cmm=0.d0
  do i=1,nat
    cmm(1,1)=cmm(1,1)+amas(i)*(csi2(i,1)**2 + csi2(i,2)**2 + csi2(i,3)**2)
    cmm(1,2)=cmm(1,2)+amas(i)*(csi1(i,2)*csi2(i,3) - csi2(i,2)*csi1(i,3))
    cmm(1,3)=cmm(1,3)+amas(i)*(csi2(i,1)*csi1(i,3) - csi1(i,1)*csi2(i,3))
    cmm(1,4)=cmm(1,4)+amas(i)*(csi1(i,1)*csi2(i,2) - csi2(i,1)*csi1(i,2))
    cmm(2,2)=cmm(2,2)+amas(i)*(csi2(i,1)**2 + csi1(i,2)**2 + csi1(i,3)**2)
    cmm(2,3)=cmm(2,3)+amas(i)*(csi2(i,1)*csi2(i,2) - csi1(i,1)*csi1(i,2))
    cmm(2,4)=cmm(2,4)+amas(i)*(csi2(i,1)*csi2(i,3) - csi1(i,1)*csi1(i,3))
    cmm(3,3)=cmm(3,3)+amas(i)*(csi1(i,1)**2 + csi2(i,2)**2 + csi1(i,3)**2)
    cmm(3,4)=cmm(3,4)+amas(i)*(csi2(i,2)*csi2(i,3) - csi1(i,2)*csi1(i,3))
    cmm(4,4)=cmm(4,4)+amas(i)*(csi1(i,1)**2 + csi1(i,2)**2 + csi2(i,3)**2)
  end do
  
  eigv=0.d0
  lworkd=4*10
  allocate(work(lworkd))
  call dsyev('V','U',4,cmm,4,eigv,work,lworkd,info)
  
  qvec(:)=cmm(:,1)

  u_mat(1,1)=qvec(1)**2 + qvec(2)**2 - qvec(3)**2 - qvec(4)**2
  u_mat(1,2)=2.d0*(qvec(2)*qvec(3) + qvec(1)*qvec(4))
  u_mat(1,3)=2.d0*(qvec(2)*qvec(4) - qvec(1)*qvec(3))
  
  u_mat(2,1)=2.d0*(qvec(2)*qvec(3) - qvec(1)*qvec(4))
  u_mat(2,2)=qvec(1)**2 + qvec(3)**2 - qvec(2)**2 - qvec(4)**2
  u_mat(2,3)=2.d0*(qvec(3)*qvec(4) + qvec(1)*qvec(2))
  
  u_mat(3,1)=2.d0*(qvec(2)*qvec(4) + qvec(1)*qvec(3))
  u_mat(3,2)=2.d0*(qvec(3)*qvec(4) - qvec(1)*qvec(2))
  u_mat(3,3)=qvec(1)**2 + qvec(4)**2 - qvec(2)**2 - qvec(3)**2


  deallocate(cmm,eigv,qvec,work,csi1,csi2)
end subroutine eckart
