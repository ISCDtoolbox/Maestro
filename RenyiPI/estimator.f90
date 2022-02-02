!#define _DEBUG
module estimator
  implicit none

  type :: avg_scalar
    integer :: num
    double precision :: summation, average
  end type avg_scalar

  type :: avg_array
    integer :: num,dimen
    double precision, allocatable :: summation(:), average(:)
  end type avg_array

  type :: avg_matrix
    integer :: num,dimen(2)
    double precision, allocatable :: summation(:,:), average(:,:)
  end type avg_matrix

  type :: esti_scalar
    integer :: num
    double precision :: summation, sumsqr, average, vari, deviation
  end type esti_scalar

  type :: esti_array
    integer :: num,dimen
    double precision, allocatable :: summation(:), sumsqr(:), average(:), vari(:), deviation(:)
  end type esti_array

  type :: esti_matrix
    integer :: num,dimen(2)
    double precision, allocatable :: summation(:,:), sumsqr(:,:), average(:,:), vari(:,:), deviation(:,:)
  end type esti_matrix

  type :: corr_basic
    integer :: num
    type(avg_array) :: estia,estib
    double precision, allocatable :: corrsum(:), corrfun(:)
  end type corr_basic

  type :: corr_advanced
    integer :: num
    type(esti_array) :: estia,estib
    double precision, allocatable :: corrsum(:), corrfun(:)
  end type corr_advanced

  interface alloc
     module procedure allocate_esti_matrix
     module procedure allocate_esti_array
     module procedure allocate_avg_array
     module procedure allocate_avg_matrix
     module procedure allocate_corr_advanced_array
     module procedure allocate_corr_basic_array
  end interface alloc

  interface free
     module procedure free_esti_matrix
     module procedure free_esti_array
     module procedure free_avg_array
     module procedure free_avg_matrix
     module procedure free_corr_advanced_array
     module procedure free_corr_basic_array
  end interface free

  interface reset
     module procedure reset_esti_scalar,reset_esti_array,reset_esti_matrix
     module procedure reset_avg_scalar,reset_avg_array,reset_avg_matrix
     module procedure reset_corr_advanced_array
     module procedure reset_corr_basic_array
  end interface reset

  interface push
     module procedure push_esti_scalar,push_esti_array,push_esti_matrix
     module procedure push_avg_scalar,push_avg_array,push_avg_matrix
     module procedure push_corr_advanced_array
     module procedure push_corr_basic_array
  end interface push

  interface calc
     module procedure calc_esti_scalar,calc_esti_array,calc_esti_matrix
     module procedure calc_avg_scalar,calc_avg_array,calc_avg_matrix
     module procedure calc_corr_advanced_array
     module procedure calc_corr_basic_array
  end interface calc

contains

!!!avg!!!
  subroutine reset_avg_scalar(estimator)
    implicit none
    type(avg_scalar) :: estimator

    estimator%num = 0
    estimator%summation = 0.d0
    estimator%average = 0.d0
  end subroutine reset_avg_scalar


  subroutine reset_avg_array(estimator)
    implicit none
    type(avg_array) :: estimator

#ifdef _DEBUG
    if(.not.allocated(estimator%summation)) then
      write(*,*) "ERROR! all the arrays for the estimator should be allocated! please call the allocate subroutine"
      stop
    endif
#endif
    estimator%num = 0
    estimator%summation = 0.d0
    estimator%average = 0.d0
  end subroutine reset_avg_array

  subroutine reset_avg_matrix(estimator)
    implicit none
    type(avg_matrix) :: estimator

#ifdef _DEBUG
    if(.not.allocated(estimator%summation)) then
      write(*,*) "ERROR! all the arrays for the estimator should be allocated! please call the allocate subroutine"
      stop
    endif
#endif
    estimator%num = 0
    estimator%summation = 0.d0
    estimator%average = 0.d0
  end subroutine reset_avg_matrix
 
  subroutine allocate_avg_array(estimator,ndim)
    implicit none
    type(avg_array) :: estimator
    integer, intent(in) :: ndim

    if(.not.allocated(estimator%summation)) then
      estimator%dimen=ndim
      allocate(estimator%summation(ndim),estimator%average(ndim))
      estimator%num = 0
      estimator%summation = 0.d0
      estimator%average = 0.d0
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't allocate the same array twice!"
#endif
    endif
  end subroutine allocate_avg_array

  subroutine allocate_avg_matrix(estimator,ndim1,ndim2)
    implicit none
    type(avg_matrix) :: estimator
    integer, intent(in) :: ndim1,ndim2

    if(.not.allocated(estimator%summation)) then
      estimator%dimen(1)=ndim1
      estimator%dimen(2)=ndim2
      allocate(estimator%summation(ndim1,ndim2),estimator%average(ndim1,ndim2))
      estimator%num = 0
      estimator%summation = 0.d0
      estimator%average = 0.d0
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't allocate the same array twice!"
#endif
    endif
  end subroutine allocate_avg_matrix

 
  subroutine push_avg_scalar(estimator,element)
    implicit none
    type(avg_scalar) :: estimator
    double precision :: element

    estimator%num = estimator%num + 1
    estimator%summation = estimator%summation + element
  end subroutine push_avg_scalar

  subroutine push_avg_array(estimator,element)
    implicit none
    type(avg_array) :: estimator
    double precision, intent(in) :: element(:)

    estimator%num = estimator%num + 1
    estimator%summation = estimator%summation + element
  end subroutine push_avg_array

  subroutine push_avg_matrix(estimator,element)
    implicit none
    type(avg_matrix) :: estimator
    double precision, intent(in) :: element(:,:)

    estimator%num = estimator%num + 1
    estimator%summation(:,:) = estimator%summation(:,:) + element(:,:)
  end subroutine push_avg_matrix

  subroutine calc_avg_scalar(estimator)
    implicit none
    type(avg_scalar) :: estimator

    if(estimator%num<=0) return
    estimator%average = estimator%summation / estimator%num
  end subroutine calc_avg_scalar

  subroutine calc_avg_array(estimator)
    implicit none
    type(avg_array) :: estimator

    if(estimator%num<=0) return
    estimator%average = estimator%summation / estimator%num
  end subroutine calc_avg_array

  subroutine calc_avg_matrix(estimator)
    implicit none
    type(avg_matrix) :: estimator

    if(estimator%num<=0) return
    estimator%average(:,:) = estimator%summation(:,:) / estimator%num
  end subroutine calc_avg_matrix

  subroutine free_avg_array(estimator)
    implicit none
    type(avg_array) :: estimator

    if(allocated(estimator%summation)) then
      deallocate(estimator%summation,estimator%average)
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't free the same array twice!"
#endif
    endif
  end subroutine free_avg_array

  subroutine free_avg_matrix(estimator)
    implicit none
    type(avg_matrix) :: estimator

    if(allocated(estimator%summation)) then
      deallocate(estimator%summation,estimator%average)
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't free the same array twice!"
#endif
    endif
  end subroutine free_avg_matrix

!!!esti!!!
  subroutine reset_esti_scalar(estimator)
    implicit none
    type(esti_scalar) :: estimator

    estimator%num = 0
    estimator%summation = 0.d0
    estimator%sumsqr = 0.d0
    estimator%average = 0.d0
    estimator%vari = 0.d0
    estimator%deviation =0.d0
  end subroutine reset_esti_scalar

  subroutine reset_esti_matrix(estimator)
    implicit none
    type(esti_matrix) :: estimator

    estimator%num = 0
    estimator%summation = 0.d0
    estimator%sumsqr = 0.d0
    estimator%average = 0.d0
    estimator%vari = 0.d0
    estimator%deviation =0.d0
  end subroutine reset_esti_matrix

  subroutine allocate_esti_array(estimator,ndim)
    implicit none
    type(esti_array) :: estimator
    integer, intent(in) :: ndim
    if(.not.allocated(estimator%summation)) then
      estimator%dimen=ndim
      allocate(estimator%summation(ndim),estimator%average(ndim),estimator%sumsqr(ndim),&
               estimator%vari(ndim),estimator%deviation(ndim))
      estimator%num = 0
      estimator%summation = 0.d0
      estimator%sumsqr = 0.d0
      estimator%average = 0.d0
      estimator%vari = 0.d0
      estimator%deviation =0.d0
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't allocate the same array twice!"
#endif
    endif
  end subroutine allocate_esti_array

  subroutine allocate_esti_matrix(estimator,ndim1,ndim2)
    implicit none
    type(esti_matrix) :: estimator
    integer, intent(in) :: ndim1,ndim2
    if(.not.allocated(estimator%summation)) then
      estimator%dimen(1)=ndim1
      estimator%dimen(2)=ndim2
      allocate(estimator%summation(ndim1,ndim2),estimator%average(ndim1,ndim2),estimator%sumsqr(ndim1,ndim2),&
               estimator%vari(ndim1,ndim2),estimator%deviation(ndim1,ndim2))
      estimator%num = 0
      estimator%summation = 0.d0
      estimator%sumsqr = 0.d0
      estimator%average = 0.d0
      estimator%vari = 0.d0
      estimator%deviation =0.d0
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't allocate the same array twice!"
#endif
    endif
  end subroutine allocate_esti_matrix


  subroutine reset_esti_array(estimator)
    implicit none
    type(esti_array) :: estimator

#ifdef _DEBUG
    if(.not.allocated(estimator%summation)) then
      write(*,*) "ERROR! all the arrays for the estimator should be allocated! please call the allocate subroutine"
      stop
    endif
#endif
    estimator%num = 0
    estimator%summation = 0.d0
    estimator%sumsqr = 0.d0
    estimator%average = 0.d0
    estimator%vari = 0.d0
    estimator%deviation =0.d0
  end subroutine reset_esti_array


  subroutine push_esti_scalar(estimator,element)
    implicit none
    type(esti_scalar) :: estimator
    double precision :: element

    estimator%num = estimator%num + 1
    estimator%summation = estimator%summation + element
    estimator%sumsqr = estimator%sumsqr + element * element
  end subroutine push_esti_scalar


  subroutine push_esti_array(estimator,element)
    implicit none
    type(esti_array) :: estimator
    double precision, intent(in) :: element(:)

    estimator%num = estimator%num + 1
    estimator%summation = estimator%summation + element
    estimator%sumsqr = estimator%sumsqr + element * element
  end subroutine push_esti_array

  subroutine push_esti_matrix(estimator,element)
    implicit none
    type(esti_matrix) :: estimator
    integer lda,ldb,ldc
    double precision, intent(in) :: element(:,:)
    double precision, allocatable :: element2(:,:)
    
    lda=max(estimator%dimen(1),1)
    ldb=max(estimator%dimen(2),1)
    ldc=max(estimator%dimen(1),1)
    
    allocate(element2(estimator%dimen(1),estimator%dimen(2)))
    element2=0.d0
    estimator%num = estimator%num + 1
    estimator%summation(:,:) = estimator%summation(:,:) + element(:,:)
!    call dgemm('N','N',estimator%dimen(1),estimator%dimen(2),estimator%dimen(2),1.d0,&
!                element,lda,element,ldb,0.d0,element2,ldc)
    estimator%sumsqr(:,:) = estimator%sumsqr(:,:) + element(:,:)**2
    deallocate(element2)
  end subroutine push_esti_matrix

  subroutine calc_esti_scalar(estimator)
    implicit none
    type(esti_scalar) :: estimator
    integer :: ii
    
    if(estimator%num<=0) return
    estimator%average = estimator%summation / estimator%num
    estimator%vari = estimator%sumsqr / estimator%num - estimator%average * estimator%average
    estimator%deviation = dsqrt(estimator%vari)
    if(estimator%vari<0.d0) write(6,*) "Warning : negative variation due to machine precision!"
  end subroutine calc_esti_scalar

  subroutine calc_esti_array(estimator)
    implicit none
    type(esti_array) :: estimator
    integer :: ii
    
    if(estimator%num<=0) return
    estimator%average = estimator%summation / estimator%num
    estimator%vari = estimator%sumsqr / estimator%num - estimator%average * estimator%average
    estimator%deviation = dsqrt(estimator%vari)
    do ii=1,estimator%dimen
      if(estimator%vari(ii)<0.d0) write(6,*) "Warning : negative variation due to machine precision!"
    enddo
  end subroutine calc_esti_array

  subroutine calc_esti_matrix(estimator)
    implicit none
    type(esti_matrix) :: estimator
    integer :: ii,jj,lda,ldb,ldc
    double precision, allocatable :: estav2(:,:)
    
    allocate(estav2(estimator%dimen(1),estimator%dimen(2)))
    if(estimator%num<=0) return
    estimator%average = estimator%summation / estimator%num

    lda=max(estimator%dimen(1),1)
    ldb=max(estimator%dimen(2),1)
    ldc=max(estimator%dimen(1),1)
    
    estav2=0.d0
    !call dgemm('N','N',estimator%dimen(1),estimator%dimen(2),estimator%dimen(2),1.d0 &
    !           ,estimator%average,lda,estimator%average,ldb,0.d0,estav2,ldc)
    
    estimator%vari(:,:) = estimator%sumsqr(:,:) / estimator%num - estimator%average(:,:)**2
    estimator%deviation = dsqrt(estimator%vari)
    do ii=1,estimator%dimen(1)
      do jj=1,estimator%dimen(2)
        if(estimator%vari(ii,jj)<0.d0) write(6,*) "Warning : negative variation due to machine precision!"
      enddo
    enddo
    deallocate(estav2)
  end subroutine calc_esti_matrix

  subroutine free_esti_array(estimator)
    implicit none
    type(esti_array) :: estimator

    if(allocated(estimator%summation)) then
      deallocate(estimator%summation,estimator%average,estimator%sumsqr,estimator%vari,estimator%deviation)
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't free the same array twice!"
#endif
    endif
  end subroutine free_esti_array

  subroutine free_esti_matrix(estimator)
    implicit none
    type(esti_matrix) :: estimator

    if(allocated(estimator%summation)) then
      deallocate(estimator%summation,estimator%average,estimator%sumsqr,estimator%vari,estimator%deviation)
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't free the same array twice!"
#endif
    endif
  end subroutine free_esti_matrix

!!!corr_basic!!!
  subroutine allocate_corr_basic_array(corr,ndima,ndimb)
    implicit none
    type(corr_basic) :: corr
    integer, intent(in) :: ndima,ndimb

    if(.not.allocated(corr%corrsum)) then
      call alloc(corr%estia,ndima)
      call alloc(corr%estib,ndimb)
      allocate(corr%corrsum(ndima*ndimb),corr%corrfun(ndima*ndimb))
      corr%num=0
      corr%corrsum=0.d0
      corr%corrfun=0.d0
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't allocate the same array twice!"
#endif
    endif
  end subroutine allocate_corr_basic_array


  subroutine reset_corr_basic_array(corr)
    implicit none
    type(corr_basic) :: corr

#ifdef _DEBUG
    if(.not.allocated(corr%corrsum)) then
      write(*,*) "ERROR! all the arrays for the estimator should be allocated! please call the allocate subroutine"
      stop
    endif
#endif
    call reset(corr%estia)
    call reset(corr%estib)
    corr%num=0
    corr%corrsum=0.d0
    corr%corrfun=0.d0
  end subroutine reset_corr_basic_array


  subroutine push_corr_basic_array(corr,elementa,elementb)
    implicit none
    type(corr_basic) :: corr
    double precision, intent(in) :: elementa(:), elementb(:)
    integer :: ii,jj

    corr%num = corr%num + 1
    call push(corr%estia,elementa)
    call push(corr%estib,elementb)
    call dgemm('N','N',corr%estia%dimen,corr%estib%dimen,1,1.d0,elementa,corr%estia%dimen,&
               elementb,1,1.d0,corr%corrsum,corr%estia%dimen)
  end subroutine push_corr_basic_array


  subroutine calc_corr_basic_array(corr)
    implicit none
    type(corr_basic) :: corr
    integer :: ii,jj

    if(corr%num<=0) return
    call calc(corr%estia)
    call calc(corr%estib)
    corr%corrfun=corr%corrsum / corr%num
    call dgemm('N','N',corr%estia%dimen,corr%estib%dimen,1,-1.d0,corr%estia%average,corr%estia%dimen,&
               corr%estib%average,1,1.d0,corr%corrfun,corr%estia%dimen)
  end subroutine calc_corr_basic_array


  subroutine free_corr_basic_array(corr)
    implicit none
    type(corr_basic) :: corr

    if(allocated(corr%corrsum)) then
      call free_avg_array(corr%estia)
      call free_avg_array(corr%estib)
      deallocate(corr%corrsum,corr%corrfun)
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't free the same array twice!"
#endif
    endif
  end subroutine free_corr_basic_array

!!!corr_advanced!!!
  subroutine allocate_corr_advanced_array(corr,ndima,ndimb)
    implicit none
    type(corr_advanced) :: corr
    integer, intent(in) :: ndima,ndimb

    if(.not.allocated(corr%corrsum)) then
      call alloc(corr%estia,ndima)
      call alloc(corr%estib,ndimb)
      allocate(corr%corrsum(ndima*ndimb),corr%corrfun(ndima*ndimb))
      corr%num=0
      corr%corrsum=0.d0
      corr%corrfun=0.d0
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't allocate the same array twice!"
#endif
    endif
  end subroutine allocate_corr_advanced_array


  subroutine reset_corr_advanced_array(corr)
    implicit none
    type(corr_advanced) :: corr

#ifdef _DEBUG
    if(.not.allocated(corr%corrsum)) then
      write(*,*) "ERROR! all the arrays for the estimator should be allocated! please call the allocate subroutine"
      stop
    endif
#endif
    call reset(corr%estia)
    call reset(corr%estib)
    corr%num=0
    corr%corrsum=0.d0
    corr%corrfun=0.d0
  end subroutine reset_corr_advanced_array


  subroutine push_corr_advanced_array(corr,elementa,elementb)
    implicit none
    type(corr_advanced) :: corr
    double precision, intent(in) :: elementa(:), elementb(:)
    integer :: ii,jj

    corr%num = corr%num + 1
    call push(corr%estia,elementa)
    call push(corr%estib,elementb)
    call dgemm('N','N',corr%estia%dimen,corr%estib%dimen,1,1.d0,elementa,corr%estia%dimen,&
               elementb,1,1.d0,corr%corrsum,corr%estia%dimen)
  end subroutine push_corr_advanced_array


  subroutine calc_corr_advanced_array(corr)
    implicit none
    type(corr_advanced) :: corr
    integer :: ii,jj

    if(corr%num<=0) return
    call calc(corr%estia)
    call calc(corr%estib)
    call dgemm('N','N',corr%estia%dimen,corr%estib%dimen,1,-1.d0,corr%estia%average,corr%estia%dimen,&
               corr%estib%average,1,1.d0,corr%corrfun,corr%estia%dimen)
  end subroutine calc_corr_advanced_array


  subroutine free_corr_advanced_array(corr)
    implicit none
    type(corr_advanced) :: corr

    if(allocated(corr%corrsum)) then
      call free_esti_array(corr%estia)
      call free_esti_array(corr%estib)
      deallocate(corr%corrsum,corr%corrfun)
#ifdef _DEBUG
    else
      write(*,*) "Warning! you can't free the same array twice!"
#endif
    endif
  end subroutine free_corr_advanced_array

end module estimator
