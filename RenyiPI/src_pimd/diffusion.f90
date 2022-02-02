program readdiff
  
  use md_variables, only : filen, lfnl, nbeadMD, nion, ndimMD 
  implicit none 

  integer :: i,j,k,icount,nbin,resto,ii,kk
  real(8)  :: var, mean, countj, wetot, diff
  integer  :: lbin,istepinit,istepmax
  real(8), dimension(:), allocatable ::  wei, meanp,meanp2,p
  real(8), dimension(:,:), allocatable:: pbuf,rion,rion_old
  logical :: ifex,yesquantum

  istepmax=0
  write(6,*) '  bin length, istepinit ?'  
  read(5,*) lbin,istepinit 
  if(istepinit.lt.0) then
     write(6,*) ' Max generation read ? '
     read(5,*) istepmax
     istepinit=-istepinit
  endif

! istepinit must be at least equal to 2
  istepinit=max(2,istepinit)
  
  var=0.d0
  countj=0.d0
  mean=0.d0

  yesquantum = .true.

  open(222,file='prefix.input',status='old',form='formatted')

  read(222,*) filen
  lfnl=index(filen,' ')-1
  write(6,*) 'name of the run: ',filen(1:lfnl)
  close(222)

! does this file exist?
  inquire(file=filen(1:lfnl)//'.rs',exist=ifex)
  if(ifex) then
     open(unit=100,file=filen(1:lfnl)//'.rs',status='old',form ='formatted')
     write(6,*) 'reading .rs file'
     read (100,*) nbeadMD, nion, ndimMD
     if(nbeadMD.eq.1) yesquantum=.false.
     close(100)
  else
     write(6,*) 'checkpoint file not found!'
     stop
  endif

  open(111,file='positions.dat',status='old')

  if(yesquantum) then
     write(6,*) 'Quantum molecular dynamics'
  else
     write(6,*) 'Classical molecular dynamics'
  endif

  allocate(rion(ndimMD,nion),rion_old(ndimMD,nion))

  if(istepmax.eq.0) then        
     i=0
     icount = 0
     do while(i.ge.0)
        read(111,*,end=100)   
        icount=icount+1
     enddo
100  continue
     istepmax=icount

  else
     icount=istepmax

  endif
      
  write(*,*) 'num lines files, tot step statistic =',icount , icount-istepinit
  if(icount-istepinit.le.0) then 
     write(*,*) 'ERROR istepinit > icount'
     stop
  endif
  nbin=(icount-istepinit)/lbin

  !write(*,*) 'bin', nbin
  if(mod(icount-istepinit,lbin).ne.0) then 
     nbin=nbin+1 
  endif
  
  allocate(pbuf(1,nbin),wei(nbin),meanp(1),meanp2(1),p(1))

  p=0.d0
  wei=lbin
  pbuf=0d0
  resto = (icount-istepinit)-(nbin-1)*lbin
  wei(nbin)=resto

  write(*,*) 'bins =', nbin
  write(*,*) 'length last bin=', resto 

  rewind(111)

  p=0.
  j=0
  k=1
  icount = 0
  do while(icount.lt.istepmax)

     read(111,*,end=200) ((rion(kk,ii),kk=1,ndimMD),ii=1,nion)
	
     icount=icount+1

     if(icount.ge.istepinit) then
 
        diff=sum((rion(:,:)-rion_old(:,:))**2)

        p(1) = p(1) + diff
        var=var+ diff**2
        countj=countj+1.d0
        mean=mean+diff
        j=j+1
        if ( k .lt. nbin .and. j .eq. lbin ) then            
           pbuf(:,k) = p(:)/wei(k)  
           p=0d0
           j=0
           k=k+1
        else if (k .eq. nbin .and. j .eq. resto) then             
           pbuf(:,k) = p(:)/wei(k)  
        endif
     endif
     rion_old=rion 
  enddo
200 continue     

  close(111)

  var=var/countj
  mean=mean/countj
  var=var-mean**2
  
  meanp=0.d0
  meanp2=0.d0
  wetot = 0.d0

  do i=1,nbin
     wetot= wetot+wei(i)
     meanp(:) = meanp(:) + pbuf(:,i)*wei(i)
     meanp2(:) = meanp2(:) + pbuf(:,i)*pbuf(:,i)*wei(i)  
  enddo
   
  meanp2 = meanp2/wetot   
  meanp = meanp/wetot

  meanp2(:)=dsqrt((meanp2(:) - meanp(:)**2)/(nbin-1))

  write(*,*) 'Diffusion atoms  +/-', meanp(1), meanp2(1)

  deallocate(pbuf,wei)
       
end program readdiff
