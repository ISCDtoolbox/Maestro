      SUBROUTINE DSYEV_MY( JOBZ, UPLO, N, A, LDA, W, INFO)
      implicit none
!
!  -- LAPACK driver routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, N,M,I
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, *), W( * )
      DOUBLE PRECISION, dimension(:,:), allocatable:: z
      double PRECISION, dimension(:), allocatable:: workp
      integer, dimension(:), allocatable:: iwork,ifail

      allocate(z(LDA,N),iwork(5*N),ifail(N),workp(8*N))

        CALL DSYEVX( JOBZ, 'A', UPLO, N, A, LDA, 0.d0, 0.d0, 1, 1,&
     &                   0.d0, M, W, Z, LDA, WORKP, 8*N, IWORK,&
     &                   IFAIL, INFO )


        do i=1,N
        A(1:N,i)=Z(1:N,i)
        enddo


      deallocate(z,iwork,ifail,workp)

      return
      END SUBROUTINE DSYEV_MY
