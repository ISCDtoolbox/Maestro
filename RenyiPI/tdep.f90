!!! this subroutine calculates the TDEP of a molecule !!!

program tdep

  use md_variables
  implicit none
  
  integer :: n_time_steps,n_row_dmd,n_unkn  !! dmd stants for displacements Molecular Dynamics
  integer :: lworkd,indi,indj,ll
  integer :: i,j,cc,l
  integer, allocatable :: idxm(:,:)  !! contiene gli indici della matrice M_dmd per ogni time step
  double precision, allocatable :: M_dmd(:,:), S_svd(:), U_svd(:,:)
  double precision, allocatable :: S_UT_svd(:,:), UT_svd(:,:)
  double precision, allocatable :: VT_svd(:,:), VTT_svd(:,:)
  double precision, allocatable :: work(:)
  double precision, allocatable :: M_pinv(:,:),phi_anahr(:),dynmat_anahr_eig(:)
  double precision, allocatable :: dynmat_force_anahr(:,:),dynmat_anahr(:,:)
  
  double precision, allocatable :: f_tdep(:)  !! forces stored for tdep
  double precision, allocatable :: d_tdep(:)  !! displacements stored for tdep
 
  call input_tools_tom
 
  n_time_steps=12
  n_unkn = ndimMD*n*(ndimMD*n+1)/2
  n_row_dmd = n_time_steps*ndimMD*n
  write(*,*) n_row_dmd,n_time_steps,n,ndimMD
  lworkd=6*n_unkn
  
  OPEN(UNIT=11,FILE='output_files_dir/forces_tdep.dat',STATUS='OLD')
  OPEN(UNIT=12,FILE='output_files_dir/positions_tdep.dat',STATUS='OLD')

  OPEN(UNIT=14,FILE='output_files_dir/phonons_tdep.dat')
  
  
  allocate(M_dmd(n_row_dmd,n_unkn),S_svd(n_unkn),U_svd(n_row_dmd,n_row_dmd))
  allocate(VT_svd(n_unkn,n_unkn),S_UT_svd(n_unkn,n_row_dmd))
  allocate(VTT_svd(n_unkn,n_unkn),UT_svd(n_row_dmd,n_row_dmd))
  allocate(M_pinv(n_unkn,n_row_dmd),phi_anahr(n_unkn))
  allocate(dynmat_force_anahr(n*ndimMD,n*ndimMD),dynmat_anahr(n*ndimMD,n*ndimMD))
  allocate(dynmat_anahr_eig(n*ndimMD))
  
  allocate(work(lworkd))
  allocate(f_tdep(n_row_dmd),d_tdep(n_row_dmd))
  allocate(idxm(n*ndimMD,n*ndimMD))
  
  M_dmd=0.d0
  S_svd=0.d0
  U_svd=0.d0; UT_svd=0.d0; S_UT_svd=0.d0
  VT_svd=0.d0; VTT_svd=0.d0
  M_pinv=0.d0; phi_anahr=0.d0
  dynmat_force_anahr=0.d0; dynmat_anahr=0.d0
  dynmat_anahr_eig=0.d0


  idxm=0
  cc=0
  do i=1,n*ndimMD
    do j=1,i
      cc=cc+1
      idxm(i,j)=cc
      idxm(j,i)=idxm(i,j)
    end do
  end do


!!! read the input forces and displacements (with respect to the initial positions!!) !!!
  do i=1,n_time_steps
    read(11,*) f_tdep(1+(i-1)*ndimMD*n : i*ndimMD*n)
    read(12,*) d_tdep(1+(i-1)*ndimMD*n : i*ndimMD*n)
  end do
  
!!! construction of the M_dmd matrix which contains the displacements and reflect the symmetry of the force constant matrix !!!
  do i=1,n_time_steps
    do j=1,n*ndimMD
      do l=1,n*ndimMD!n_unkn
        cc=idxm(j,l)
        M_dmd((i-1)*n*ndimMD+j,cc)=d_tdep((i-1)*n*ndimMD+l)
      enddo
    enddo
  enddo

  call dgesvd('A','A',n_row_dmd,n_unkn,M_dmd,n_row_dmd,S_svd,U_svd,n_row_dmd, &
              VT_svd,n_unkn,work,lworkd,info)
              
!!! invert the diagonal S matrix coming from svd decomposition              
  do i=1,n_unkn
    if(S_svd(i).ne.0.d0) S_svd(i)=1.d0/S_svd(i)
  end do

!!! make the transpose of all the matrices coming from the SVD decomposition !!!   

!! V matrix ---> VTT=V !!  
  do i=1,n_unkn
    do j=1,n_unkn
      VTT_svd(j,i)=VT_svd(i,j)
    end do
  end do
!! U matrix
  do i=1,n_row_dmd
    do j=1,n_row_dmd
      UT_svd(j,i)=U_svd(i,j)
    end do
  end do

!!! questa è una prova per evitare una moltiplicazione tra matrici !!!
  do i=1,n_row_dmd
    if (i .le. n_unkn) then
      do j=1,n_row_dmd
        S_UT_svd(i,j)=UT_svd(i,j)*S_svd(j)
      end do
    end if
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
!!! this operation produces the matrix M_pinv, that is the pseudoinverse of the initial matrix M_dmd
  call dgemm('N','N',n_unkn,n_row_dmd,n_unkn,1.d0,VTT_svd,n_unkn,S_UT_svd,n_unkn,0.d0,M_pinv,n_unkn)


!!! finally I multiply the pseudoinverse M_pinv for the force vector f_tdep
!!! this is a (n_unkn,n_row_dmd)*(n_row_dmd) operation
  do i=1,n_unkn
    do j=1,n_row_dmd
        phi_anahr(i)=phi_anahr(i)+M_pinv(i,j)*f_tdep(j)
    end do
  end do  
  
  do i=1,n*ndimMD
    do j=1,n*ndimMD
      cc=idxm(i,j)
      dynmat_force_anahr(i,j)=phi_anahr(cc)
    end do
  end do
  
  
     
!! dividing the dynmat_anahr matrix for the sqrt of the masses     
  do j=1,n
    do ll=1,ndimMD
      indj=(j-1)*ndimMD+ll
      dynmat_anahr(indj,indj) =  dynmat_force_anahr(indj,indj)/amas(indx(j))
! off-diag part (equal particles)
      do l=ll+1,ndimMD
        indi=(j-1)*ndimMD+l 
        dynmat_anahr(indi,indj) = dynmat_force_anahr(indi,indj)/sqrt(amas(indx(j))*amas(indx(j))) 
      enddo
! off-diag part (unequal particles, unequal coordinates)
      do i=j+1,n
        do l=1,ndimMD
          indi=(i-1)*ndimMD+l 
          dynmat_anahr(indi,indj) = dynmat_force_anahr(indi,indj)/sqrt(amas(indx(i))*amas(indx(j)))    
        enddo
      enddo
    enddo
  enddo

! Write out the dynamical matrix
!  open(12,file='zz_dinmat')
!  do indi=1,n*ndimMD
!     do indj=1,n*ndimMD        
!        write(*,*) indi,indj,dynmat(indi,indj)   
!     enddo
!   enddo

! Diagonalization of the dynamical matrix
    lworkd = 3*n*ndimMD
  
    if(allocated(work)) deallocate(work)
    allocate(work(lworkd))
  
    call dsyev('V','L', n*ndimMD, dynmat_anahr, n*ndimMD, dynmat_anahr_eig, work, lworkd, info)  

    if(info.ne.0) then
      write(6,*) 'some problem in dsyev from dynmatrix routine: info =',info
      stop
    endif

    dynmat_anahr_eig(:)=abs(dynmat_anahr_eig(:))   ! The matrix is positive definite

    write(6,*) 'Harmonic frequencies in cm^-1'
    do i=7,n*ndimMD
       write(6,*) i,sqrt(dynmat_anahr_eig(i))*ha2cm1
       write(14,*) i,sqrt(dynmat_anahr_eig(i))*ha2cm1
    enddo
  
  
  close(11)
  close(12)
  close(14)
  
  deallocate(M_dmd,S_svd,U_svd,VT_svd,work)
  deallocate(f_tdep,d_tdep,idxm)
  deallocate(S_UT_svd,VTT_svd,UT_svd)
  deallocate(M_pinv,phi_anahr,dynmat_force_anahr,dynmat_anahr)
  deallocate(dynmat_anahr_eig)
  deallocate(nnucl,amas,ion_name,indx)
  
  return
end program tdep
