      subroutine dscalzero(n,zero,vet,m)
        implicit none
        integer n,m,i
        real*8 zero,vet(m,*)
        do i=1,n
        vet(1,i)=zero
        enddo
        return
      END
