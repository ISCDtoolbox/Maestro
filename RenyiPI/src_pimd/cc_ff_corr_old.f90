program cc_ff_autocorr

  use estimator
  use md_variables
  implicit none

  integer n_therm,cc,len_file,nstep_after_therm
  integer i1,i2,i3,i4,i5,lworkd,skip,lbin
  integer cc1,cc2
  double precision, allocatable :: f_file(:,:), v_file(:,:), c_file(:,:), t_file(:)
  double precision, allocatable :: pp_loc(:,:),ccm1_loc(:,:), vvm1_loc(:,:), ff_loc(:,:)
  double precision, allocatable :: diag_mm(:,:),diag_mmbuf(:,:),mbuf(:,:)
  double precision, allocatable :: work(:),centre_of_mass_coord(:)
  double precision, allocatable :: x_new(:,:), p_new(:,:), f_new(:,:)
  double precision, allocatable :: eigv_corr_loc(:,:), eigv_corr_main(:,:)
  double precision, allocatable :: vecprov(:), vecprovbuf(:), checkpos(:)
  double precision, allocatable :: x_eq(:,:), x_loc(:,:), p_loc(:,:), f_loc(:,:), umat(:,:)
  
  double precision, allocatable :: S_svd(:), U_svd(:,:)
  double precision, allocatable :: VT_svd(:,:)
  double precision, allocatable :: proj_U_Sm1(:,:), proj_VT(:,:)
  
  type(avg_scalar) :: avg_temp
  type(esti_array) :: avg_x
  type(esti_matrix) :: eigv_corr
  type(avg_matrix) :: ccm1mat,ccm1matJN
  type(avg_matrix) :: ffmat,ffmatJN
  type(avg_matrix) :: ppmat,ppmatJN
  type(avg_matrix) :: vvm1mat,vvm1matJN

  call input_tools_tom
  skip=5
  n_therm=10*nstep/iprint
  nstep_after_therm=nstep*(nblocks-60)/iprint-n_therm
  lbin=1000
   
  write(*,*) '-----------------------------------------------------------------------------------------------------------'
  write(*,*) 'Calculations of phonons corrected by the <rr> and <ff> covariance matrices with Eckart rotations of vectors'
  write(*,*) '-----------------------------------------------------------------------------------------------------------'
  write(*,*)
  write(*,*) '# of thermalizion steps ', n_therm
  write(*,*) '# of dynamical steps ', nstep_after_therm
  write(*,*) 'Temperature in Kelvin: ', tempMD
  write(*,*)

  OPEN(UNIT=12,FILE='configinit',STATUS='OLD',iostat=iflagerr)
  if(iflagerr.ne.0) then
    write(*,*) 'Opening configinit error!'
    stop
  endif
  
  OPEN(UNIT=13,FILE='forces.dat',STATUS='OLD',iostat=iflagerr)
  if(iflagerr.ne.0) then
    write(*,*) 'Opening forces.dat error!'
    stop
  endif
  
  OPEN(UNIT=14,FILE='positions.dat',STATUS='OLD',iostat=iflagerr)
  if(iflagerr.ne.0) then
    write(*,*) 'Opening positions.dat error!'
    stop
  endif
  
  OPEN(UNIT=15,FILE='velocities.dat',STATUS='OLD',iostat=iflagerr)
  if(iflagerr.ne.0) then
    write(*,*) 'Opening velocities.dat error!'
    stop
  endif

  OPEN(UNIT=16,FILE='local_temp.dat',STATUS='OLD',iostat=iflagerr)
  if(iflagerr.ne.0) then
    write(*,*) 'Opening local_temp.dat error!'
    stop
  endif

  OPEN(UNIT=25,FILE='cc_ff_spectrum.dat')  

!!! compute the length of vectors which depends on how many configurations you want to skip
  len_file=0
  do i1=1,n_therm
    read(13,*)
  end do
  do i1=1,nstep_after_therm,skip
    len_file=len_file+1
    read(13,'(400e15.7)')
  end do
  rewind(13)
  write(*,*) 'test over file lenght: ',len_file,nstep_after_therm/skip

!!!---------------------------------------------------------------------------------------------
!!!-------------------------allocation and inizialization section-------------------------------
  allocate(f_file(ndimMD*n,len_file))
  allocate(c_file(ndimMD*n,len_file))
  allocate(v_file(ndimMD*n,len_file))
  allocate(t_file(len_file))

  allocate(diag_mm(n*ndimMD,n*ndimMD),diag_mmbuf(n*ndimMD,n*ndimMD),mbuf(n*ndimMD,n*ndimMD))
  allocate(x_new(n,ndimMD),p_new(n,ndimMD),f_new(n,ndimMD))
  allocate(eigv_corr_loc(n*ndimMD,4),eigv_corr_main(n*ndimMD,4)) 
  allocate(x_eq(n,ndimMD),x_loc(n,ndimMD),p_loc(n,ndimMD),f_loc(n,ndimMD),umat(ndimMD,ndimMD))
  allocate(vecprov(ndimMD),vecprovbuf(ndimMD),checkpos(n*ndimMD))
  allocate(pp_loc(ndimMD*n,ndimMD*n),ff_loc(ndimMD*n,ndimMD*n),centre_of_mass_coord(ndimMD))

  allocate(S_svd(n*ndimMD), U_svd(n*ndimMD,n*ndimMD))
  allocate(VT_svd(n*ndimMD,n*ndimMD))
  allocate(ccm1_loc(n*ndimMD,n*ndimMD),vvm1_loc(n*ndimMD,n*ndimMD))
  
  call alloc(ccm1mat,n*ndimMD,n*ndimMD)
  call alloc(vvm1mat,n*ndimMD,n*ndimMD)
  call alloc(ffmat,n*ndimMD,n*ndimMD)
  call alloc(avg_x,ndimMD*n)
  call alloc(ppmat,n*ndimMD,n*ndimMD)
  call reset(ccm1mat)
  call reset(vvm1mat)
  call reset(ffmat)
  call reset(avg_x)
  call reset(ppmat)
  call reset(avg_temp)
  
  f_file=0.d0
  c_file=0.d0
  v_file=0.d0
  t_file=0.d0
  
  diag_mm=0.d0; diag_mmbuf=0.d0; mbuf=0.d0
  x_new=0.d0; p_new=0.d0; f_new=0.d0;
  eigv_corr_loc=0.d0; eigv_corr_main=0.d0
  x_eq=0.d0; x_loc=0.d0; p_loc=0.d0; umat=0.d0
  vecprov=0.d0; vecprovbuf=0.d0; checkpos=0.d0
  pp_loc=0.d0; ff_loc=0.d0; centre_of_mass_coord=0.d0
  
  S_svd=0.d0; U_svd=0.d0
  VT_svd=0.d0; ccm1_loc=0.d0; vvm1_loc=0.d0
  
!!!-------------------------end of inizialization and allocation section------------------------
!!!---------------------------------------------------------------------------------------------

  
!!! diag_mm is the trivial mass matrix   
  do i1=1,n
    diag_mm(i1*3-2,i1*3-2)=amas(indx(i1))
    diag_mm(i1*3-1,i1*3-1)=amas(indx(i1))
    diag_mm(i1*3,i1*3)=amas(indx(i1))
  enddo

!!!---------------------------------------------------------------------------------------------
!!! ---------------------reading the file section----------------------------------------------- 
!!! jump the first # n_therm thermalization steps...  
  do i1=1,n_therm
    read(13,*)
    read(14,*)
    read(15,*)
    read(16,*)
  end do
!!! ...and read files and skip some configurations to have decorrelated configurations (this must be tested with autocorr.f90!)
  cc=0
  do i1=1,nstep_after_therm,skip
    cc=cc+1
    read(13,'(400e15.7)') f_file(:,cc)
    read(14,'(400e15.7)') c_file(:,cc)
    read(15,'(400e15.7)') v_file(:,cc)
    read(16,'(e15.7)') t_file(cc)
!!! transform velocities to momenta  
    do i2=1,n
      v_file(1+ndimMD*(i2-1):ndimMD*i2,cc)=v_file(1+ndimMD*(i2-1):ndimMD*i2,cc)*amas(indx(i2))
    end do
  end do
!!!---------------------end reading file section------------------------------------------------  
!!!---------------------------------------------------------------------------------------------

!!!---------------------------------------------------------------------------------------------
!!!------------creating the initial equilibrium position vector---------------------------------
  x_eq=0.d0 
  do i1=1,n
    read(12,*) x_eq(i1,1),x_eq(i1,2),x_eq(i1,3)
  enddo
  centre_of_mass_coord=0.d0
  do i2=1,n
    centre_of_mass_coord(1)=centre_of_mass_coord(1)+amas(i2)*x_eq(i2,1)
    centre_of_mass_coord(2)=centre_of_mass_coord(2)+amas(i2)*x_eq(i2,2)
    centre_of_mass_coord(3)=centre_of_mass_coord(3)+amas(i2)*x_eq(i2,3)
  end do
  centre_of_mass_coord(:)=centre_of_mass_coord(:)/sum(amas)  
  do i2=1,n
    x_eq(i2,1)=x_eq(i2,1)-centre_of_mass_coord(1)
    x_eq(i2,2)=x_eq(i2,2)-centre_of_mass_coord(2)
    x_eq(i2,3)=x_eq(i2,3)-centre_of_mass_coord(3)
  end do
!!!------------end the initialitazion of equilibrium position vector----------------------------   
!!!---------------------------------------------------------------------------------------------    

!  open(unit=81,file='prova_phys.xyz')
!  open(unit=82,file='prova_eck.xyz')
   
!!! --------------------------------------------------------------------------------------------    
!!! ---------------loop over files to create temperatures and coordinates averages--------------
  do i1=1,len_file
    
    call push(avg_temp,t_file(i1))    
    
    centre_of_mass_coord=0.d0
    do i2=1,n
      centre_of_mass_coord(1)=centre_of_mass_coord(1)+amas(i2)*c_file(1+(i2-1)*ndimMD,i1)
      centre_of_mass_coord(2)=centre_of_mass_coord(2)+amas(i2)*c_file(2+(i2-1)*ndimMD,i1)
      centre_of_mass_coord(3)=centre_of_mass_coord(3)+amas(i2)*c_file(3+(i2-1)*ndimMD,i1)
    end do
    centre_of_mass_coord(:)=centre_of_mass_coord(:)/sum(amas)
    
    x_loc=0.d0
    do i2=1,n
      x_loc(i2,1)=c_file(1+(i2-1)*ndimMD,i1)-centre_of_mass_coord(1)
      x_loc(i2,2)=c_file(2+(i2-1)*ndimMD,i1)-centre_of_mass_coord(2)
      x_loc(i2,3)=c_file(3+(i2-1)*ndimMD,i1)-centre_of_mass_coord(3)
    end do
    
    umat=0.d0
!!! compute the rotation matrix 'umat'...    
    call eckart(n,ndimMD,x_eq,x_loc,amas,umat)
    
    checkpos=0.d0
    do i2=1,n
      do i3=1,ndimMD
        do i4=1,ndimMD
          checkpos(i3+(i2-1)*ndimMD)=checkpos(i3+(i2-1)*ndimMD)+umat(i3,i4)*x_loc(i2,i4)
        end do
      end do
    end do
    
    call push(avg_x,checkpos)    
  
  end do
  call calc(avg_x)
  call calc(avg_temp)
!!! --------------------------------------------------------------------------------------------    
!!! ---------------end of calculation of temperature and coordinate averages--------------------


!!! --------------------------------------------------------------------------------------------    
!!! ---------------main loop: accumlate the <r_i r_j> matrix------------------------------------  
  do i1=1,len_file
    
    centre_of_mass_coord=0.d0
    do i2=1,n
      centre_of_mass_coord(1)=centre_of_mass_coord(1)+amas(i2)*c_file(1+(i2-1)*ndimMD,i1)
      centre_of_mass_coord(2)=centre_of_mass_coord(2)+amas(i2)*c_file(2+(i2-1)*ndimMD,i1)
      centre_of_mass_coord(3)=centre_of_mass_coord(3)+amas(i2)*c_file(3+(i2-1)*ndimMD,i1)
    end do
    centre_of_mass_coord(:)=centre_of_mass_coord(:)/sum(amas)
    
    x_loc=0.d0
    do i2=1,n
      x_loc(i2,1)=c_file(1+(i2-1)*ndimMD,i1)-centre_of_mass_coord(1)
      x_loc(i2,2)=c_file(2+(i2-1)*ndimMD,i1)-centre_of_mass_coord(2)
      x_loc(i2,3)=c_file(3+(i2-1)*ndimMD,i1)-centre_of_mass_coord(3)
    end do

    p_loc=0.d0
    do i2=1,n
      p_loc(i2,1)=v_file(1+(i2-1)*ndimMD,i1)
      p_loc(i2,2)=v_file(2+(i2-1)*ndimMD,i1)
      p_loc(i2,3)=v_file(3+(i2-1)*ndimMD,i1)
    end do

    f_loc=0.d0
    do i2=1,n
      f_loc(i2,1)=f_file(1+(i2-1)*ndimMD,i1)
      f_loc(i2,2)=f_file(2+(i2-1)*ndimMD,i1)
      f_loc(i2,3)=f_file(3+(i2-1)*ndimMD,i1)
    end do
   
    umat=0.d0
!!! compute the rotation matrix 'umat'...    
    call eckart(n,ndimMD,x_eq,x_loc,amas,umat)
    
!!! check the Eckart conditions...
!!! **********************************************************************
    checkpos=0.d0; x_new=0.d0; p_new=0.d0; f_new=0.d0
    do i2=1,n
      do i3=1,ndimMD
        do i4=1,ndimMD
          checkpos(i3+(i2-1)*ndimMD)=checkpos(i3+(i2-1)*ndimMD)+umat(i3,i4)*x_loc(i2,i4)
          x_new(i2,i3)=x_new(i2,i3)+umat(i3,i4)*x_loc(i2,i4)
          p_new(i2,i3)=p_new(i2,i3)+umat(i3,i4)*p_loc(i2,i4)
          f_new(i2,i3)=f_new(i2,i3)+umat(i3,i4)*f_loc(i2,i4)
        end do
      end do
    end do

!    write(81,*) 7
!    write(81,*) 'cartesian trajectory' 
!    write(82,*) 7
!    write(82,*) 'eckart trajectory' 
!    do i2=1,n
!      if(i2.le.2) write(81,*) 'O',x_loc(i2,:) 
!      if(i2.gt.2) write(81,*) 'H',x_loc(i2,:) 
!      if(i2.le.2) write(82,*) 'C',checkpos(1+(i2-1)*ndimMD:3+(i2-1)*ndimMD) 
!      if(i2.gt.2) write(82,*) 'He',checkpos(1+(i2-1)*ndimMD:3+(i2-1)*ndimMD) 
!    end do

    vecprov=0.d0
    do i2=1,n
      vecprovbuf=0.d0
      call r3cross(x_eq(i2,:),checkpos((1+(i2-1)*ndimMD):(3+(i2-1)*ndimMD)),vecprovbuf(:))
      vecprov(:)=vecprov(:)+amas(i2)*vecprovbuf(:)
    end do
    if (sqrt(sum(vecprov(:)**2)) .gt. 1.d-6) then
      write(*,*) 'you are not in the Eckart frame...error!!'
      stop
    end if
!!!************************************************************************        
    
!!! ...and rotate the coordinates to eliminate the translational and rotational degrees of freedom   
    ccm1_loc=0.d0; pp_loc=0.d0; ff_loc=0.d0; vvm1_loc=0.d0
    cc1=0
    do i2=1,n
      do i3=1,ndimMD
        cc1=cc1+1
        cc2=0
        do i4=1,n
          do i5=1,ndimMD
            cc2=cc2+1
              
              ccm1_loc(cc1,cc2)=(x_new(i2,i3)-avg_x%average(cc1))*(x_new(i4,i5)-avg_x%average(cc2))
              pp_loc(cc1,cc2)=p_new(i2,i3)*p_new(i4,i5)
              vvm1_loc(cc1,cc2)=p_new(i2,i3)/amas(i2)*p_new(i4,i5)/amas(i4)
              ff_loc(cc1,cc2)=f_new(i2,i3)*f_new(i4,i5)
          
          end do
        end do 
      end do
    end do

!!! accumulate the <r_i r_j> symmetric matrix, <p_i p_j> (generalized mass matrix)
!!! and <f_i f_j>   
    call push(ccm1mat,ccm1_loc)
    call push(ppmat,pp_loc)
    call push(vvm1mat,vvm1_loc)
    call push(ffmat,ff_loc)
    
  end do
  call calc(ccm1mat)
  call calc(ppmat)
  call calc(vvm1mat)
  call calc(ffmat)
!!! ---------------------------end main loop----------------------------------------------------
!!! --------------------------------------------------------------------------------------------
!  close(81)
!  close(82)
  
  write(*,*) 'average temperature: ',avg_temp%average

!!! average the <r_i r_j> and <p_i p_j> matrices over all the decorrelated steps
  ccm1_loc=0.d0
  pp_loc=0.d0
  ff_loc=0.d0
  vvm1_loc=0.d0
  ccm1_loc=ccm1mat%average
  pp_loc=ppmat%average
  ff_loc=ffmat%average
  vvm1_loc=vvm1mat%average

!!!---------------------------------------------------------------------------------------------
!!!---------------------------invert the <r_i r_j> matrix with an SVD decomposition-------------

  lworkd = 100*n*ndimMD
  if(allocated(work)) deallocate(work)
  allocate(work(lworkd))
  call dgesvd('A','A',n*ndimMD,n*ndimMD,ccm1_loc,n*ndimMD,S_svd,U_svd,n*ndimMD, &
              VT_svd,n*ndimMD,work,lworkd,info)
 
!!! invert the diagonal S matrix coming from svd decomposition              
  cc=0
  do i1=1,n*ndimMD
   ! if(S_svd(i1).ne.0.d0) S_svd(i1)=1.d0/S_svd(i1)
    if(S_svd(i1).gt.1.d-6) then
       S_svd(i1)=1.d0/S_svd(i1)
    else
      cc=cc+1
      S_svd(i1)=0.d0
    endif

  end do
  allocate(proj_U_Sm1(ndimMD*n,ndimMD*n-cc),proj_VT(ndimMD*n-cc,ndimMD*n))
  proj_U_Sm1=0.d0
  proj_VT=0.d0

!!! project over the subspaces with non-zero values... 
  do i1=1,n*ndimMD
     do i2=1,n*ndimMD
        if(i2 .le. ndimMD*n-cc) proj_U_Sm1(i1,i2)=U_svd(i1,i2)*S_svd(i2)
        if(i1 .le. ndimMD*n-cc) proj_VT(i1,i2)=VT_svd(i1,i2)
     end do
  end do
 
  ccm1_loc=0.d0
  call dgemm('N','N',n*ndimMD,n*ndimMD,ndimMD*n-cc,1.d0,proj_U_Sm1,n*ndimMD,&
             proj_VT,n*ndimMD-cc,0.d0,ccm1_loc,n*ndimMD)
  deallocate(proj_U_Sm1,proj_VT)
  
  !!!---------------------------------------------------------------------------------------------
!!!---------------------------invert the <v_i v_j> matrix with an SVD decomposition-------------
  S_svd=0.d0; U_svd=0.d0; VT_svd=0.d0
  lworkd = 100*n*ndimMD
  if(allocated(work)) deallocate(work)
  allocate(work(lworkd))
  call dgesvd('A','A',n*ndimMD,n*ndimMD,vvm1_loc,n*ndimMD,S_svd,U_svd,n*ndimMD, &
              VT_svd,n*ndimMD,work,lworkd,info)
  
!!! invert the diagonal S matrix coming from svd decomposition              
  cc=0
  do i1=1,n*ndimMD
    write(*,*) S_svd(i1)
    if(S_svd(i1).gt.1.d-10) then
       S_svd(i1)=1.d0/S_svd(i1)
    else
      cc=cc+1
      S_svd(i1)=0.d0
    endif

  end do
  allocate(proj_U_Sm1(ndimMD*n,ndimMD*n-cc),proj_VT(ndimMD*n-cc,ndimMD*n))
  proj_U_Sm1=0.d0
  proj_VT=0.d0

!!! project over the subspaces with non-zero values... 
  do i1=1,n*ndimMD
     do i2=1,n*ndimMD
        if(i2 .le. ndimMD*n-cc) proj_U_Sm1(i1,i2)=U_svd(i1,i2)*S_svd(i2)
        if(i1 .le. ndimMD*n-cc) proj_VT(i1,i2)=VT_svd(i1,i2)
     end do
  end do
 
  vvm1_loc=0.d0
  call dgemm('N','N',n*ndimMD,n*ndimMD,ndimMD*n-cc,1.d0,proj_U_Sm1,n*ndimMD,&
             proj_VT,n*ndimMD-cc,0.d0,vvm1_loc,n*ndimMD)
  deallocate(proj_U_Sm1,proj_VT)     
!!!-----------------------------end of matrix inversion-----------------------------------------
!!!---------------------------------------------------------------------------------------------

 
 ! do i1=1,ndimMD*n
 !    do i2=1,ndimMD*n
 !      write(*,*) i1,i2,ccm1_loc(i1,i2)
 !    end do
 !    write(*,*)
 ! end do

!!!---------------------------------------------------------------------------------------------  
!!!-------------------finally solve the generalized eigenvalue problem-------------------------- 	
  mbuf=ccm1_loc*(avg_temp%average/kbm1)
  diag_mmbuf=diag_mm
  call localise(diag_mmbuf,mbuf,eigv_corr_main(:,1),n*ndimMD,1)
  eigv_corr_main(:,1)=sqrt(abs(eigv_corr_main(:,1)))

  mbuf=ccm1_loc
  diag_mmbuf=vvm1_loc
  call dpotrf('U',n*ndimMD,diag_mmbuf,n*ndimMD,info)
  if(info .ne. 0) then
     write(*,*) 'matrix vvm1_loc is not positive definite..error!'
     stop
  end if
  call localise(vvm1_loc,mbuf,eigv_corr_main(:,2),n*ndimMD,1)
  eigv_corr_main(:,2)=sqrt(abs(eigv_corr_main(:,2)))
  
  mbuf=ff_loc/(avg_temp%average/kbm1)
  diag_mmbuf=diag_mm
  call localise(diag_mmbuf,mbuf,eigv_corr_main(:,3),n*ndimMD,1)
  eigv_corr_main(:,3)=sqrt(abs(eigv_corr_main(:,3)))
  
  mbuf=ff_loc
  diag_mmbuf=pp_loc
  call dpotrf('U',n*ndimMD,diag_mmbuf,n*ndimMD,info)
  if(info .ne. 0) then
     write(*,*) 'matrix pp_loc is not positive definite..error!'
     stop
  end if
  call localise(pp_loc,mbuf,eigv_corr_main(:,4),n*ndimMD,1)
  eigv_corr_main(:,4)=sqrt(abs(eigv_corr_main(:,4)))  
  
  
  write(6,*) 'Corrected Harmonic frequencies in cm^-1'
  do i1=1,n*ndimMD
     write(6,*) i1,eigv_corr_main(i1,:)*ha2cm1
  enddo


!!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!*********************************************************************************************
!!!---------------------------------------------------------------------------------------------
!!!-------------------------now estimate the error with Jacknife-------------------------------

  write(*,*) 
  write(*,*)
  write(*,*) 'Now estimate the error with Jacknife'


!!! now run the Jackknife method to compute the error
  call alloc(ccm1matJN,n*ndimMD,n*ndimMD)
  call alloc(vvm1matJN,n*ndimMD,n*ndimMD)
  call alloc(ppmatJN,n*ndimMD,n*ndimMD)
  call alloc(ffmatJN,n*ndimMD,n*ndimMD)
  call alloc(eigv_corr,n*ndimMD,4)
  
  call reset(ppmatJN)
  call reset(ffmatJN)
  call reset(ccm1matJN)
  call reset(vvm1matJN)
  call reset(eigv_corr)

  do i1=1,len_file
          
    centre_of_mass_coord=0.d0
    do i2=1,n
      centre_of_mass_coord(1)=centre_of_mass_coord(1)+amas(i2)*c_file(1+ndimMD*(i2-1),i1)
      centre_of_mass_coord(2)=centre_of_mass_coord(2)+amas(i2)*c_file(2+ndimMD*(i2-1),i1)
      centre_of_mass_coord(3)=centre_of_mass_coord(3)+amas(i2)*c_file(3+ndimMD*(i2-1),i1)
    end do
    centre_of_mass_coord(:)=centre_of_mass_coord(:)/sum(amas)
    
    x_loc=0.d0
    do i2=1,n
      x_loc(i2,1)=c_file(1+(i2-1)*ndimMD,i1)-centre_of_mass_coord(1)
      x_loc(i2,2)=c_file(2+(i2-1)*ndimMD,i1)-centre_of_mass_coord(2)
      x_loc(i2,3)=c_file(3+(i2-1)*ndimMD,i1)-centre_of_mass_coord(3)
    end do

    p_loc=0.d0
    do i2=1,n
      p_loc(i2,1)=v_file(1+(i2-1)*ndimMD,i1)
      p_loc(i2,2)=v_file(2+(i2-1)*ndimMD,i1)
      p_loc(i2,3)=v_file(3+(i2-1)*ndimMD,i1)
    end do

    f_loc=0.d0
    do i2=1,n
      f_loc(i2,1)=f_file(1+(i2-1)*ndimMD,i1)
      f_loc(i2,2)=f_file(2+(i2-1)*ndimMD,i1)
      f_loc(i2,3)=f_file(3+(i2-1)*ndimMD,i1)
    end do
    
    umat=0.d0
!!! compute the rotation matrix 'umat'...    
    call eckart(n,ndimMD,x_eq,x_loc,amas,umat)
    
    x_new=0.d0; p_new=0.d0; f_new=0.d0 
    do i2=1,n
      do i3=1,ndimMD
        do i4=1,ndimMD
          x_new(i2,i3)=x_new(i2,i3)+umat(i3,i4)*x_loc(i2,i4)
          p_new(i2,i3)=p_new(i2,i3)+umat(i3,i4)*p_loc(i2,i4)
          f_new(i2,i3)=f_new(i2,i3)+umat(i3,i4)*f_loc(i2,i4)
        end do
      end do
    end do
    
!!! store the local matrices <ccm1> and <pp> 
    ccm1_loc=0.d0; pp_loc=0.d0; vvm1_loc=0.d0; ff_loc=0.d0   
    cc1=0
    do i2=1,n
      do i3=1,ndimMD
        cc1=cc1+1
        cc2=0
        do i4=1,n
          do i5=1,ndimMD
            cc2=cc2+1
            ccm1_loc(cc1,cc2)=(x_new(i2,i3)-avg_x%average(cc1))*(x_new(i4,i5)-avg_x%average(cc2))
            pp_loc(cc1,cc2)=p_new(i2,i3)*p_new(i4,i5)
            vvm1_loc(cc1,cc2)=p_new(i2,i3)/amas(i2)*p_new(i4,i5)/amas(i4)
            ff_loc(cc1,cc2)=f_new(i2,i3)*f_new(i4,i5)
          end do
        end do 
      end do
    end do
    call push(ccm1matJN,ccm1_loc)
    call push(ppmatJN,pp_loc)
    call push(vvm1matJN,vvm1_loc)
    call push(ffmatJN,ff_loc)

    ccm1_loc=0.d0; pp_loc=0.d0; vvm1_loc=0.d0; ff_loc=0.d0
    if(ccm1matJN%num .eq. lbin) then
         ccm1_loc(:,:)=(ccm1mat%summation(:,:)-ccm1matJN%summation(:,:))/(ccm1mat%num-ccm1matJN%num)
         pp_loc(:,:)=(ppmat%summation(:,:)-ppmatJN%summation(:,:))/(ppmat%num-ppmatJN%num)
         ff_loc(:,:)=(ffmat%summation(:,:)-ffmatJN%summation(:,:))/(ffmat%num-ffmatJN%num)
         vvm1_loc(:,:)=(vvm1mat%summation(:,:)-vvm1matJN%summation(:,:))/(vvm1mat%num-vvm1matJN%num)
         
!!!---------------------------------------------------------------------------------------------
!!!---------------------------invert the <r_i r_j> matrix with an SVD decomposition-------------

         lworkd = 100*n*ndimMD
         if(allocated(work)) deallocate(work)
         allocate(work(lworkd))
         call dgesvd('A','A',n*ndimMD,n*ndimMD,ccm1_loc,n*ndimMD,S_svd,U_svd,n*ndimMD, &
                     VT_svd,n*ndimMD,work,lworkd,info)
 
!!! invert the diagonal S matrix coming from svd decomposition              
         cc=0
         do i2=1,n*ndimMD
   ! if(S_svd(i1).ne.0.d0) S_svd(i1)=1.d0/S_svd(i1)
           if(S_svd(i2).gt.1.d-10) then
              S_svd(i2)=1.d0/S_svd(i2)
           else
             cc=cc+1
             S_svd(i2)=0.d0
           endif
         end do
         allocate(proj_U_Sm1(ndimMD*n,ndimMD*n-cc),proj_VT(ndimMD*n-cc,ndimMD*n))
         proj_U_Sm1=0.d0
         proj_VT=0.d0

!!! project over the subspaces with non-zero values... 
         do i2=1,n*ndimMD
           do i3=1,n*ndimMD
             if(i3 .le. ndimMD*n-cc) proj_U_Sm1(i2,i3)=U_svd(i2,i3)*S_svd(i3)
             if(i2 .le. ndimMD*n-cc) proj_VT(i2,i3)=VT_svd(i2,i3)
           end do
         end do
 
         ccm1_loc=0.d0
         call dgemm('N','N',n*ndimMD,n*ndimMD,ndimMD*n-cc,1.d0,proj_U_Sm1,n*ndimMD,&
                    proj_VT,n*ndimMD-cc,0.d0,ccm1_loc,n*ndimMD)
         deallocate(proj_U_Sm1,proj_VT)
  
!!!---------------------------------------------------------------------------------------------
!!!---------------------------invert the <v_i v_j> matrix with an SVD decomposition-------------
         S_svd=0.d0; U_svd=0.d0; VT_svd=0.d0
         lworkd = 100*n*ndimMD
         if(allocated(work)) deallocate(work)
         allocate(work(lworkd))
         call dgesvd('A','A',n*ndimMD,n*ndimMD,vvm1_loc,n*ndimMD,S_svd,U_svd,n*ndimMD, &
                     VT_svd,n*ndimMD,work,lworkd,info)
  
!!! invert the diagonal S matrix coming from svd decomposition              
         cc=0
         do i2=1,n*ndimMD
   ! if(S_svd(i1).ne.0.d0) S_svd(i1)=1.d0/S_svd(i1)
           if(S_svd(i2).gt.1.d-10) then
             S_svd(i2)=1.d0/S_svd(i2)
           else
             cc=cc+1
             S_svd(i2)=0.d0
           endif
         end do
         allocate(proj_U_Sm1(ndimMD*n,ndimMD*n-cc),proj_VT(ndimMD*n-cc,ndimMD*n))
         proj_U_Sm1=0.d0
         proj_VT=0.d0

!!! project over the subspaces with non-zero values... 
         do i2=1,n*ndimMD
           do i3=1,n*ndimMD
             if(i3 .le. ndimMD*n-cc) proj_U_Sm1(i2,i3)=U_svd(i2,i3)*S_svd(i3)
             if(i2 .le. ndimMD*n-cc) proj_VT(i2,i3)=VT_svd(i2,i3)
           end do
         end do
 
         vvm1_loc=0.d0
         call dgemm('N','N',n*ndimMD,n*ndimMD,ndimMD*n-cc,1.d0,proj_U_Sm1,n*ndimMD,&
                    proj_VT,n*ndimMD-cc,0.d0,vvm1_loc,n*ndimMD)
         deallocate(proj_U_Sm1,proj_VT)     
!!!-----------------------------end of matrix inversion-----------------------------------------
!!!---------------------------------------------------------------------------------------------
         eigv_corr_loc=0.d0
         mbuf=ccm1_loc*(avg_temp%average/kbm1)
         diag_mmbuf=diag_mm

         call localise(diag_mmbuf,mbuf,eigv_corr_loc(:,1),n*ndimMD,1)
         eigv_corr_loc(:,1)=sqrt(abs(eigv_corr_loc(:,1)))

         mbuf=ccm1_loc
         diag_mmbuf=vvm1_loc
         call dpotrf('U',n*ndimMD,diag_mmbuf,n*ndimMD,info)
         if(info .ne. 0) then
           write(*,*) 'matrix vvm1_loc is not positive definite..error!'
           stop
         end if
         call localise(vvm1_loc,mbuf,eigv_corr_loc(:,2),n*ndimMD,1)
         eigv_corr_loc(:,2)=sqrt(abs(eigv_corr_loc(:,2)))
  
         mbuf=ff_loc/(avg_temp%average/kbm1)
         diag_mmbuf=diag_mm
         call localise(diag_mmbuf,mbuf,eigv_corr_loc(:,3),n*ndimMD,1)
         eigv_corr_loc(:,3)=sqrt(abs(eigv_corr_loc(:,3)))
  
         mbuf=ff_loc
         diag_mmbuf=pp_loc
         call dpotrf('U',n*ndimMD,diag_mmbuf,n*ndimMD,info)
         if(info .ne. 0) then
           write(*,*) 'matrix pp_loc is not positive definite..error!'
           stop
         end if
         call localise(pp_loc,mbuf,eigv_corr_loc(:,4),n*ndimMD,1)
         eigv_corr_loc(:,4)=sqrt(abs(eigv_corr_loc(:,4))) 
         
         call push(eigv_corr,eigv_corr_loc(:,:))
       
         call reset(ccm1matJN)
         call reset(vvm1matJN)
         call reset(ffmatJN)
         call reset(ppmatJN)
     endif
  end do
  
  call calc(eigv_corr)

  write(6,*) 'Corrected Harmonic frequencies in cm^-1'
  write(25,*) 'Corrected Harmonic frequencies in cm^-1'
  do i1=1,n*ndimMD
     write(6,*) i1,eigv_corr%average(i1,:)*ha2cm1,'plus or minus',eigv_corr%deviation(i1,:)*ha2cm1
     write(25,*) i1,tempMD,eigv_corr%average(i1,:)*ha2cm1,eigv_corr%deviation(i1,:)*ha2cm1
  enddo

 
 
  call free(eigv_corr)
  call free(avg_x)
  call free(ccm1mat)
  call free(ccm1matJN)
  call free(ppmat)
  call free(ppmatJN)
  call free(vvm1mat)
  call free(vvm1matJN)
  call free(ffmat)
  call free(ffmatJN)
  

  deallocate(f_file)
  deallocate(c_file)
  deallocate(v_file)
  deallocate(t_file)

  deallocate(diag_mm,diag_mmbuf,mbuf)
  deallocate(x_new,p_new,f_new)
  deallocate(eigv_corr_loc,eigv_corr_main) 
  deallocate(x_eq,x_loc,p_loc,f_loc,umat)
  deallocate(vecprov,vecprovbuf,checkpos)
  deallocate(pp_loc,ff_loc,centre_of_mass_coord)

  deallocate(S_svd,U_svd)
  deallocate(VT_svd)
  deallocate(ccm1_loc)

  deallocate(nnucl,amas,ion_name,indx)
  
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(25)
 
end program cc_ff_autocorr


subroutine localise(K0,K2,frequencies,ndim,gee_type)
    implicit none
    integer, intent(in) :: ndim,gee_type
    double precision, dimension(ndim,ndim), intent(in) :: K0,K2
    double precision, dimension(ndim), intent(out) :: frequencies

    double precision, dimension(ndim,ndim) :: K0_loc,K2_loc
    double precision, dimension(:), allocatable :: work
    double precision :: query
    integer :: lwork,info

    K0_loc=K0
    K2_loc=K2

    lwork=-1
    call dsygv(gee_type,'V','U',ndim,K2_loc,ndim,K0_loc,ndim,&
                         frequencies,query,lwork,info)

    lwork=floor(query)
    allocate(work(lwork))

    call dsygv(gee_type,'V','U',ndim,K2_loc,ndim,K0_loc,ndim,&
                         frequencies,work,lwork,info)
    deallocate(work)
    return
end subroutine localise
