subroutine save_res(H,eigen,d,t,ratio,at1,at2)
  
  implicit none

  !INPUT VARIABLES
  integer, intent(in) :: d
  double precision, intent(in) :: H(d,d), eigen(d)
  double precision, intent(in) :: ratio
  double precision, intent(in) :: at1, at2
  character(len=4), intent(in) :: t
  

  !LOCAL VARIABLES
  character(len=256) :: folder, filename, dirpath, fullpath, which
  character(len=8) :: typ1, typ2, beta
  integer :: uw, i, j

  ! CASTING FROM DIGITS TO CHARACTER/STRING
  write(typ1,'(F8.3)') at1
  typ1=adjustl(typ1)
  write(typ2,'(F8.3)') at2
  typ2=adjustl(typ2)
  write(beta,'(F8.3)') ratio
  beta=adjustl(beta)
  write(folder, '(I0)') d
  if ( t .eq. 'L' ) then ;
    write(which,'(A)') "linear"
  else
    write(which, '(A)') "cyclic"
  end if

  ! CREATION OF THE FOLDER PATH 
  
  dirpath = "eigenvectors/"//trim(which)//"/"//trim(typ1)//"_"//trim(typ2)//"/"//trim(folder)//"_"//trim(beta)
  call execute_command_line('mkdir -p "'//trim(dirpath)//'"')
  
  ! GRAPHICAL BEGIN OF SAVING PROCEDURE
  write(*,*)
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!                    WAIT: WRITING THE EIGENVECTORS...                        !!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,*)

  ! PRINTING EIGENVALUES
  fullpath=trim(dirpath)//"/eigenvalues.txt"
  open(unit=11, file=fullpath, status='unknown')

  do i = 1, d
    write(11,*) (i-1)/dble(d-1), eigen(i) !writing
  end do
  close(11) 

  ! PRINTING EIGENVECTORS

  do i = 1 , d
    uw=i*20
    write(filename, '(A,I0)') 'eigenvector_',i
    fullpath = trim(dirpath)//"/"//trim(filename)
    open(unit=uw, file=trim(fullpath), status='unknown')
    do j=1, d
      write(uw,'(F10.5)') H(j,i)
    end do
    close(uw)
  end do

  write(*,'(1X,A,A)') "THE EIGENVECTORS HAVE BEEN COPIED IN THE FOLDER: ", trim(dirpath)
end subroutine save_res 
