!!! if you want to sample the highest frequencies, 
!!! the following relation for delta_t (delta in our input file) must hold 
!!! 1/delta_t > 2*w_max  ....this is the Nyquist-Shannon theorem 

program vv_autocorr

  use md_variables
  use estimator
  implicit none

  integer ntau,n_therm,cc,nstep_vv
  integer i1,i2,i3,i4,lworkd
  double precision x1,x2,y1,y2,z1,z2,dtau,wmax,wmin,domega,norm
  double precision, allocatable :: v_file(:,:)
  double precision, allocatable :: Gtau_cl(:),Gtau_corr(:)
  double precision, allocatable :: mbuf(:,:),wr(:),wi(:),work(:),vl(:,:),vr(:,:)
  double precision, allocatable :: ft_Gtau_cl(:),ft_Gtau_corr(:)
  
  type(avg_matrix),allocatable :: Gtau_mm(:)
  type(avg_matrix) mnorm

  call input_tools_tom
  
  dtau=iprint*delt*timeau2fs
  
  n_therm=nstep/iprint
  nstep_vv=nstep*(nblocks-4)/iprint
  ntau=5000!(nstep*(nblocks-4)/iprint)
  
  write(*,*) '-----------------------------------------------------------------------------'
  write(*,*) 'Calculations of power spectrum by velocity-velocity correlation matrix method'
  write(*,*) '-----------------------------------------------------------------------------'
  write(*,*)
  write(*,*) '# of thermalizion steps ', n_therm
  write(*,*) '# of dynamical steps ', nstep_vv
  write(*,*) '# of discretized tau ', ntau
  write(*,*) 'dt (time step of the printed files) expressed in fs ', dtau
  write(*,*)

  OPEN(13,file='velocities.dat',STATUS='OLD',iostat=iflagerr)
  if(iflagerr.ne.0) then
    write(*,*) 'Opening velocities.dat error!'
    stop
  endif
  OPEN(UNIT=14,FILE='gtau.dat')  

  allocate(v_file(ndimMD*n,nstep_vv))
  allocate(Gtau_cl(ntau),ft_Gtau_cl(ntau))
  allocate(Gtau_corr(ntau),ft_Gtau_corr(ntau),mbuf(n*ndimMD,n*ndimMD))
  allocate(wr(n*ndimMD),wi(n*ndimMD),vr(n*ndimMD,n*ndimMD),vl(n*ndimMD,n*ndimMD))
  
  allocate(Gtau_mm(ntau))
  do i1=1,ntau
    call alloc(Gtau_mm(i1),n*ndimMD,n*ndimMD)
    call reset(Gtau_mm(i1))
  end do 
  call alloc(mnorm,n*ndimMD,n*ndimMD)
  call reset(mnorm)

!!! jump the first # n_therm thermalization steps !!! 
  do i1=1,n_therm
    read(13,*)
  end do
 
  do i1=1,nstep_vv
    read(13,'(400e15.7)') v_file(:,i1)
  end do

 
  Gtau_cl=0.d0
  Gtau_corr=0.d0

!!!-----------------------------------------------------------------!!!
!!! loov over tau to create the function g(tau)= <v(t_0)v(t_0+tau)>   
!!!-----------------------------------------------------------------!!!
  do i1=1,ntau   
     write(*,*) i1

!!! loop = integral (sum!) over t_0
     do i2=1,nstep_vv   
!!! normalization factor is <v(t_0)v(t_0)>: I can done this operation only at the first step       
        mbuf=0.d0
        if (i1 .eq. 1)  then
          do i3=1,n*ndimMD
            do i4=1,n*ndimMD
               mbuf(i3,i4) = v_file(i3,i2)*v_file(i4,i2)
            end do
          end do
          call push(mnorm,mbuf)          
        end if
        
        mbuf=0.d0
        if (i1-1+i2 .le. nstep_vv) then   
          do i3=1,n*ndimMD
            do i4=1,n*ndimMD
               mbuf(i3,i4) = v_file(i3,i2)*v_file(i4,(i1-1)+i2)
            end do
          end do
          call push(Gtau_mm(i1),mbuf)
        end if
        
     end do

  end do 
  
  do i1=1,ntau
    call calc(Gtau_mm(i1))
  end do
  call calc(mnorm)
  
!!! now I distinguish two cases: sum of the diagonal part only of the g(tau) matrix ---> g(tau)_clas
!!! and the sum over the eigenvalues of the whole matrix ---> g(tau)_corr  
  do i1=1,ntau
    mbuf=0.d0
    Gtau_mm(i1)%average(:,:)=Gtau_mm(i1)%average(:,:)/(n*ndimMD*mnorm%average(:,:))
    do i2=1,n*ndimMD
      Gtau_cl(i1)=Gtau_cl(i1)+Gtau_mm(i1)%average(i2,i2)
      do i3=1,n*ndimMD     
        mbuf(i2,i3)=Gtau_mm(i1)%average(i2,i3)
      end do
    end do
    
    lworkd = 3*n*ndimMD
    if(allocated(work)) deallocate(work)
    allocate(work(lworkd))
    call dgeev('N','N',n*ndimMD,mbuf,ndimMD*n,wr,wi,vl,ndimMD*n,vr,ndimMD*n,work,lworkd,info)
    Gtau_corr(i1)=Gtau_corr(i1)+sum(wr(:))    
  end do
  
  do i1=1,ntau
    write(14,*) i1*dtau,Gtau_cl(i1),Gtau_corr(i1)
  end do 
 
!!! compute the maximum and minimum frequencies and convert to cm^(-1)
  wmax=hbar2pi*(1.d0/(1.d0*dtau))*ev2cm1
  wmin=hbar2pi*(1.d0/(1.d0*ntau*dtau))*ev2cm1
  domega=(wmax-wmin)/ntau
 
  write(*,*) 'w (frequency) max is: ', wmax
  write(*,*) 'w (frequency) min is: ', wmin
  write(*,*) 'delta omega is: ', domega 
  write(*,*)
  
  call ft_write_spectrum('clas',Gtau_cl,ft_Gtau_cl,ntau,domega)
  call ft_write_spectrum('corr',Gtau_corr,ft_Gtau_corr,ntau,domega) 
  
  deallocate(v_file)
  deallocate(Gtau_cl,ft_Gtau_cl)
  deallocate(Gtau_corr,ft_Gtau_corr,mbuf)
  deallocate(wr,wi,work,vl,vr)
  
  do i1=1,ntau
    call free(Gtau_mm(i1))
  end do
  deallocate(Gtau_mm)
  call free(mnorm)
  
  deallocate(nnucl,amas,ion_name,indx)
  
  close(13)
  close(14)
 
end program vv_autocorr


subroutine ft_write_spectrum(clas_or_corr,input,output,nn,domega)
    
  implicit none
  character*4 :: clas_or_corr
  integer :: nn,iw,i,j
  double precision dw,omega,domega
  double precision, dimension(nn) :: input
  double precision, dimension(nn) :: output
  double precision gen_ker
  double precision icompl
  double precision, parameter :: pi=4.d0*atan(1.d0)
     
  do i=1,nn
    output(i)=0.d0
    do j=0,nn-1
      icompl=cos(2*pi*j*i/nn)
      output(i)=output(i)+input(j+1)*icompl*gen_ker(j,nn)
    end do
  end do

   open(9,file='vv_power_spectrum_'//clas_or_corr//'.dat')
     dw=domega
     do iw=20,nn/2+1
        omega=dw*dble(iw)
        write(9,*) omega,real(output(iw))
     end do
   close(9)

end subroutine ft_write_spectrum

!!! Jakson kernel to filter the spectrum
double precision function gen_ker(ixn,ixnn)
  implicit none
  integer ixn,ixnn
  double precision pi
  pi=4.d0*atan(1.d0)
  gen_ker=((ixnn-(ixn-1)+1)*cos(pi*(ixn-1)/(ixnn+1))+sin(pi*(ixn-1)/(ixnn+1))*&
           cos(pi/(ixnn+1))/sin(pi/(ixnn+1)))/(ixnn+1)
end function

