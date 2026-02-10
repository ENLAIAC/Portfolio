subroutine HL_gap(t, ratio, d, eigen,at1,at2)
  
  implicit none

  ! INPUT VARIABLES BLOCK
  character(len=4), intent(in) :: t
  double precision, intent(in) :: ratio, at1, at2 !beta1/beta2 ratio absolute value
  integer, intent(in) :: d
  double precision, intent(in) :: eigen(d)

  ! LOCAL VARIABLES 
  integer :: n_el
  character (len=256) :: filename, fullpath, folder
  character (len=8) :: typ1, typ2
  
  n_el=d/2
  write(typ1,'(F8.3)') at1
  typ1=adjustl(typ1)
  write(typ2,'(F8.3)') at2
  typ2=adjustl(typ2)
  
  select case (t)
    case ('C')
      if ( at2 .ne. at1 .or. ratio .lt. 1.0d0 ) then ;
        write(*,*)
        write(*,'(1X,A,F12.6)') "THE VALUE OF THE HOMO LUMO GAP IS: ", eigen(n_el+1)-eigen(n_el)
        write(filename, '(A,F4.2)') 'gap_cycl_'//trim(typ1)//"_"//trim(typ2)//"_", ratio
        open(unit=15, file=filename, status='unknown', access='append')
        write(15,'(I4.4,F12.6)') d, eigen(n_el+1) - eigen(n_el)
     else if ( mod(d,4) .eq. 0) then;
        write(*,*)
        write(*,'(1X,A,F12.6)') "THE VALUE OF THE HOMO LUMO GAP IS: ", eigen(n_el+2)-eigen(n_el+1)
        write(filename, '(A,F4.2)') 'gap_cycl_'//trim(typ1)//"_"//trim(typ2)//"_", ratio
        open(unit=15, file=filename, status='unknown', access='append')
        write(15,'(I4.4,F12.6)') d, eigen(n_el+2) - eigen(n_el+1)
     else
        write(*,*)
        write(*,'(1X,A,F12.6)') "THE VALUE OF THE HOMO LUMO GAP IS: ", eigen(n_el+1)-eigen(n_el)
        write(filename, '(A,F4.2)') 'gap_cycl_'//trim(typ1)//"_"//trim(typ2)//"_", ratio
        open(unit=15, file=filename, status='unknown', access='append')
        write(15,'(I4.4,F12.6)') d, eigen(n_el+1) - eigen(n_el)
     end if
    case default
      write(*,*)
      write(*,'(1X,A,F12.6)') "THE VALUE OF THE HOMO LUMO GAP IS: ", eigen(n_el+1)-eigen(n_el)
      write(filename, '(A,F4.2)') 'gap_lin_'//trim(typ1)//"_"//trim(typ2)//"_", ratio
      open(unit=15, file=filename, status='unknown', access='append')
      write(15,'(I4.4,F12.6)') d, eigen(n_el+1) - eigen(n_el)
  end select
  close(15)
  write(*,'(1X,A,A)') "THE HOMO-LUMO GAP VALUE WAS APPEND IN THE FILE: ", filename 

end subroutine HL_gap
