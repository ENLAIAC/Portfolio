!------------------------------------------------------------------------
subroutine diagonalize_matrix(N,A,e)

      ! Diagonalize a square matrix

      implicit none

      ! Input variables

      integer,intent(in)            :: N
      double precision,intent(inout):: A(N,N)
      double precision,intent(out)  :: e(N)

      ! Local variables

      integer                       :: lwork,info
      integer                       :: i
      double precision,allocatable  :: work(:)


      ! Memory allocation

      allocate(work(3*N))
      lwork = size(work)

      call dsyev('V','U',N,A,N,e,work,lwork,info)

      if(info /= 0) then
        write(*,'(a)') 'Problem in diagonalize_matrix (dsyev)!!'
        stop
      endif

      do i = 1 , N
        if (abs(e(i)) < 1e-16) e(i) = 0
      end do

end subroutine diagonalize_matrix
!------------------------------------------------------------------------
