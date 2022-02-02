! the dynamical matrix is the mass-reduced Fourier transform of the force constant matrix!
! but for molecules q=0 only !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dynmatrix(rpos_loc,yeswrite_loc)
  
  use md_variables 

  implicit none 
  integer :: i,j,k,l,ll,indi,indj,lworkd
  real(8) :: ep,ep_plus,ep_minus,deltasq
  real(8) :: ep_plxply,ep_plxmy,ep_mxply,ep_mxmy
  real(8) :: delta, zeta_rand
  real(8) :: drand1,drand2
  real(8) :: rpos_loc(ndimMD,n)
!  real(8), parameter :: delta=1.d-2
  real(8), allocatable :: dynmat_force(:,:),work(:)
  real(8), parameter :: dyn_cut = 10000.d0  ! cm^{-1} units
  logical :: yeswrite_loc

  delta=delta_harm
  deltasq=delta**2

  allocate(dynmat_force(n*ndimMD,n*ndimMD))

  do i=1,n 
     do l=1,ndimMD
        rtilde(i,l) = rpos_loc(l,i)
     enddo
  enddo
  
  call calcpot(ep,rtilde) 
        
      
! Computing dynamical matrix
! lower triangle
  
  do j=1,n
     do ll=1,ndimMD

        indj=(j-1)*ndimMD+ll

! diagonal part (equal particles, equal coordinates)

        rtilde(j,ll) = rpos_loc(ll,j) + delta
        call calcpot(ep_plus,rtilde)
 
        rtilde(j,ll) = rpos_loc(ll,j) - delta
        call calcpot(ep_minus,rtilde)
        
        rtilde(j,ll) = rpos_loc(ll,j)
        dynmat_force(indj,indj) = (ep_plus+ep_minus-2.d0*ep)/deltasq
        dynmat(indj,indj) =  dynmat_force(indj,indj)/amas(indx(j))

! add diagonal term to dynmat to avoid negative eigenvalues ????????????
        dynmat_force(indj,indj) =  dynmat_force(indj,indj) + dyn_cut*amas(indx(j))/ha2cm1**2
        dynmat(indj,indj) =  dynmat(indj,indj) + dyn_cut/ha2cm1**2

! off-diag part (equal particles)

        do l=ll+1,ndimMD
           indi=(j-1)*ndimMD+l 
           
           rtilde(j,l) = rpos_loc(l,j) + delta
           rtilde(j,ll) = rpos_loc(ll,j) + delta
           call calcpot(ep_plxply,rtilde)
                        
           rtilde(j,l) = rpos_loc(l,j) - delta
           call calcpot(ep_mxply,rtilde)
                 
           rtilde(j,ll) = rpos_loc(ll,j) - delta
           call calcpot(ep_mxmy,rtilde)
                 
           rtilde(j,l) = rpos_loc(l,j) + delta                 
           call calcpot(ep_plxmy,rtilde)
                 
                 
           dynmat_force(indi,indj) = (ep_plxply+ep_mxmy-ep_plxmy-ep_mxply)/4.d0/deltasq
           dynmat(indi,indj) = dynmat_force(indi,indj)/sqrt(amas(indx(j))*amas(indx(j))) 

           rtilde(j,l) = rpos_loc(l,j)
           rtilde(j,ll) = rpos_loc(ll,j)
        enddo

! off-diag part (unequal particles, unequal coordinates)
        do i=j+1,n
           do l=1,ndimMD
              indi=(i-1)*ndimMD+l 

              rtilde(i,l) = rpos_loc(l,i) + delta
              rtilde(j,ll) = rpos_loc(ll,j) + delta
              call calcpot(ep_plxply,rtilde)
                        
              rtilde(i,l) = rpos_loc(l,i) - delta
              call calcpot(ep_mxply,rtilde)
                 
              rtilde(j,ll) = rpos_loc(ll,j) - delta
              call calcpot(ep_mxmy,rtilde)
                 
              rtilde(i,l) = rpos_loc(l,i) + delta                 
              call calcpot(ep_plxmy,rtilde)
          
               
              dynmat_force(indi,indj) = (ep_plxply+ep_mxmy-ep_plxmy-ep_mxply)/4.d0/deltasq
              dynmat(indi,indj) = dynmat_force(indi,indj)/sqrt(amas(indx(i))*amas(indx(j)))    

              rtilde(i,l) = rpos_loc(l,i)
              rtilde(j,ll) = rpos_loc(ll,j)

           enddo
        enddo

     enddo
  enddo

! Write out the dynamical matrix
!  do indi=1,n*ndimMD
!     do indj=1,n*ndimMD
!        write(6,*) indi,indj,dynmat(indi,indj)
!     enddo
!  enddo
 
! Save force-constant matrix before it is destroyed by the diagonalization
  !dynmat_force0=dynmat_force

 
! Diagonalization of the dynamical matrix
  lworkd = 3*n*ndimMD
  
  if(allocated(work)) deallocate(work)
  allocate(work(lworkd))
  
  call dsyev('V','L', n*ndimMD, dynmat, n*ndimMD, dynmat_eig, work, lworkd, info)  


  if(info.ne.0) then
     write(6,*) 'some problem in dsyev from dynmatrix routine: info =',info
     stop
  endif

  dynmat_eig(:)=abs(dynmat_eig(:))   ! The matrix is positive definite

!  write(6,*) 'Dynamical matrix eigenvalues are'
!  do i=1,n*ndimMD
!     write(6,*) i,dynmat_eig(i)
!  enddo

!  if(yeswrite_loc) then
!     write(6,*) 'Harmonic frequencies in cm^-1'
!     do i=1,n*ndimMD
!        write(6,*) i,sqrt(dynmat_eig(i))*ha2cm1
!     enddo
!  endif

! Diagonalization of force-constant matrix  
! dynmat_force contains the eigenvectors after the diagonalization
  call dsyev('V','L', n*ndimMD, dynmat_force, n*ndimMD, dynmatforce_eig, work, lworkd, info)  
  deallocate(work)

  dynmatforce_eig(:)=abs(dynmatforce_eig(:))   ! The matrix is positive definite

  if(info.ne.0) then
     write(6,*) 'some problem in dsyev (force-constant matrix) from dynmatrix routine: info =',info
     stop
  endif

!!! sigmavar=sigmacov*sigmacov

! Noise on the CC forces
  do k=1,n*ndimMD
     call random_number(drand1)
     call random_number(drand2)
     zeta_rand=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2)
     do i=1,n*ndimMD
!        fk(k,kp_ion+i)=dynmat(i,k)*sqrt(dynmat_eig(k))
!        zeta_rand=dsqrt(-2.d0*dlog(1.d0-drand1()))*dcos(2.d0*pi*drand1())
!        fk(k,kp_ion+i)=dynmat_force(i,k)*zeta_rand*sqrt(dynmat_eig(k)*sigmavar)        

!        if (irun .eq. 2) then 
!           fk(k,i)=dynmat(i,k)*zeta_rand*sqrt(sigmavar*dynmat_eig(k))
!        elseif (irun .eq. 3) then
           fk(k,i)=dynmat_force(i,k)*zeta_rand*sqrt(sigmavar*dynmatforce_eig(k))
!        endif 
     enddo
!     write(6,*) k,fnoiseMD(k) 
  enddo
  

! Covariance matrix of the CC forces
!  covMD=0.d0 
!  do i=1,n*ndimMD
!     do j=1,n*ndimMD
!        do k=1,n*ndimMD 
!           covMD(i,j) = covMD(i,j) + dynmat_force(i,k)*dynmat_eig(k)*dynmat_force(j,k)
!        enddo 
!     enddo
!  enddo 

!  covMD=sigmavar*covMD 
     
!   Write out the dynamical matrix
!  do indi=1,n*ndimMD
!     do indj=1,n*ndimMD
!        write(6,*) indi,indj,covMD(indi,indj)
!     enddo
!  enddo


  deallocate(dynmat_force)

!  stop 

  return         
end subroutine dynmatrix
