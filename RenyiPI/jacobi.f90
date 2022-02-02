subroutine jacobi(nat,ndim,x_loc,amas,x_new)
  
  implicit none
  integer i1,i2,nat,ndim
  double precision amas(nat)
  double precision x_loc(nat,ndim),x_new(nat,ndim)
  x_new=0.d0
  do i1=1,nat-1
    do i2=1,i1
      x_new(i1,:)=x_new(i1,:)+amas(i2)*x_loc(i2,:)/sum(amas(1:i1))
    end do
    x_new(i1,:)=x_new(i1,:)-x_loc(i1+1,:)
  end do
  do i1=1,ndim
    x_new(nat,i1)=sum(amas(:)*x_loc(:,i1))/sum(amas(:)) !0.d0
  end do
end subroutine jacobi


subroutine gmatrix(nat,ndim,amas,gmat)
  implicit none
  integer i1,i2,i3,i4,i5,i6,nat,ndim,cc1,cc2,cc3
  double precision, allocatable :: bvec(:,:)
  double precision gmat(nat*ndim,nat*ndim),amas(nat)
  
  allocate(bvec(nat*ndim,nat*ndim))
  bvec=0.d0

!!! ds_1/dx_(1:nat)  
  bvec(1,1)=1.d0; bvec(2,2)=1.d0; bvec(3,3)=1.d0
  bvec(1,4)=-1.d0; bvec(2,5)=-1.d0; bvec(3,6)=-1.d0

!!! ds_(2:nat-1)/dx_(1:nat)  
  cc1=3
  do i1=2,nat-1  
    do i2=1,ndim
      cc1=cc1+1
      cc2=mod(cc1,3)
      if (mod(cc1,3).eq.0) cc2=3
      do i3=1,i1
        bvec(cc1,cc2)=amas(i3)/sum(amas(1:i1))
        cc2=cc2+3
      end do
      bvec(cc1,cc1+3)=-1.d0
    end do
     
  end do
  
!!! ds_(nat)/dx_(1:nat)  
  do i1=1,nat  
      bvec(nat*ndim-2,1+(i1-1)*ndim)=amas(i1)/sum(amas(1:nat))
      bvec(nat*ndim-1,2+(i1-1)*ndim)=amas(i1)/sum(amas(1:nat))
      bvec(nat*ndim,3+(i1-1)*ndim)=amas(i1)/sum(amas(1:nat))
  end do

  do i1=1,nat*ndim
    do i2=1,nat*ndim
      write(*,*) i1,i2,bvec(i1,i2)
    end do
    write(*,*)
  end do
  
  gmat=0.d0
  cc1=0
  do i1=1,nat
    do i2=1,ndim
      cc1=cc1+1
      cc2=0
      do i3=1,nat
        do i4=1,ndim
          cc2=cc2+1
          cc3=0
          do i5=1,nat
            do i6=1,ndim 
              cc3=cc3+1
              gmat(cc1,cc2)=gmat(cc1,cc2)+1.d0/amas(i5)*bvec(cc1,cc3)*bvec(cc2,cc3)
            end do
          end do
        end do  
      end do
    end do
  end do
  
  
  deallocate(bvec)
end subroutine


subroutine align_vector_to_z(ndim,x_loc,rotmatv1)
  implicit none
  integer ndim
  double precision :: x_loc(ndim),normx,normr,rotmatv1(ndim,ndim)
  double precision, allocatable :: x_buf(:),zzz(:)
  
  allocate(x_buf(ndim),zzz(ndim))
  x_buf=0.d0
  zzz(3)=1.d0; zzz(1:2)=0.d0
  normx=sqrt(sum(x_loc(:)**2))
  x_buf(:)=x_loc(:)/normx
  normr=sum((x_buf(:)+zzz(:))**2)
  
  rotmatv1=0.d0
  rotmatv1(1,1)=2*(x_buf(1)+zzz(1))*(x_buf(1)+zzz(1))/normr-1.d0
  rotmatv1(1,2)=2*(x_buf(1)+zzz(1))*(x_buf(2)+zzz(2))/normr
  rotmatv1(1,3)=2*(x_buf(1)+zzz(1))*(x_buf(3)+zzz(3))/normr
  
  rotmatv1(2,1)=rotmatv1(1,2)
  rotmatv1(2,2)=2*(x_buf(2)+zzz(2))*(x_buf(2)+zzz(2))/normr-1.d0
  rotmatv1(2,3)=2*(x_buf(2)+zzz(2))*(x_buf(3)+zzz(3))/normr
  
  rotmatv1(3,1)=rotmatv1(1,3)
  rotmatv1(3,2)=rotmatv1(2,3)
  rotmatv1(3,3)=2*(x_buf(3)+zzz(3))*(x_buf(3)+zzz(3))/normr-1.d0
  
  deallocate(x_buf,zzz)
end subroutine


subroutine place_vector_to_xz(ndim,x_loc,rotmatp2)
  implicit none
  integer ndim
  double precision :: x_loc(ndim),angle,normx,rotmatp2(ndim,ndim)
  double precision, allocatable :: x_buf(:),yyy(:)
  
  allocate(x_buf(ndim),yyy(ndim))
  x_buf=0.d0 
  yyy(1)=0.d0; yyy(2)=1.d0; yyy(3)=0.d0
  normx=sqrt(sum(x_loc(:)**2))
  x_buf(:)=x_loc(:)/normx
  
  angle=abs(sum(x_buf(:)*yyy(:)))/(abs(sum(x_buf))*abs(sum(yyy)))
  angle=asin(angle)
  
  rotmatp2=0.d0
  rotmatp2(1,1)=cos(-angle); rotmatp2(1,2)=-sin(-angle)
  rotmatp2(2,1)=sin(-angle); rotmatp2(2,2)=cos(-angle)
  rotmatp2=1.d0
  
  deallocate(x_buf,yyy)
end subroutine
