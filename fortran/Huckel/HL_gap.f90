subroutine HL_gap(t, ratio, d, eigen,at1,at2)
  
  implicit none

  ! INPUT VARIABLES BLOCK

  character(len=4), intent(in) :: t
  double precision, intent(in) :: ratio, at1, at2 !ratio = beta1/beta2 ratio absolute value
  integer, intent(in) :: d 
  double precision, intent(in) :: eigen(d)

  ! LOCAL VARIABLES

  integer :: n_coup
  character (len=256) :: filename
  character (len=8) :: typ1, typ2

  ! VARIABLES INITIALIZATION

  n_coup=d/2 !amount of electron couples
  write(typ1,'(F8.3)') at1 ! converting the atom type to char
  typ1=adjustl(typ1) ! adjusting spaces in the just converted character variable
  write(typ2,'(F8.3)') at2 ! converting the atom type to char
  typ2=adjustl(typ2) !adjusting the length and the spaces in the variable

  ! HOMO-LUMO CALCULATION

  select case (t) 
    case ('C') ! if type is cyclic do:
      if ( at2 .ne. at1 .or. ratio .lt. 1.0d0 ) then ; ! case 1: cyclic chain with bond alternation
        write(*,*)
        write(*,'(1X,A,F12.6)') "THE VALUE OF THE HOMO LUMO GAP IS: ", eigen(n_coup+1)-eigen(n_coup) ! displaying the gap for case 1
        write(filename, '(A,F4.2)') 'gap_cycl_'//trim(typ1)//"_"//trim(typ2)//"_", ratio ! generating the file name where the gap
        open(unit=15, file=filename, status='unknown', access='append')                  ! info will be stored
        write(15,'(I4.4,F12.6)') d, eigen(n_coup+1) - eigen(n_coup)                      ! writing it in the just named file
      else if ( mod(d,4) .eq. 0) then; !case 2: cyclic with no bond alternation (whichever) and an antiaromatic structure
        write(*,*)
        write(*,'(1X,A,F12.6)') "THE VALUE OF THE HOMO LUMO GAP IS: ", eigen(n_coup+2)-eigen(n_coup+1) ! display the gap for case 2
        write(filename, '(A,F4.2)') 'gap_cycl_'//trim(typ1)//"_"//trim(typ2)//"_", ratio ! generating the name of the file
        open(unit=15, file=filename, status='unknown', access='append') ! opening the file unit
        write(15,'(I4.4,F12.6)') d, eigen(n_coup+2) - eigen(n_coup+1) ! writing the gap in the file
      else ! case 3 (last cyclic case): cyclic, no dimerization, and aromatic 
        write(*,*)
        write(*,'(1X,A,F12.6)') "THE VALUE OF THE HOMO LUMO GAP IS: ", eigen(n_coup+1)-eigen(n_coup) ! compute and display the gap
        write(filename, '(A,F4.2)') 'gap_cycl_'//trim(typ1)//"_"//trim(typ2)//"_", ratio ! name of the file to store the gap
        open(unit=15, file=filename, status='unknown', access='append') ! opening the corresponding unit
        write(15,'(I4.4,F12.6)') d, eigen(n_coup+1) - eigen(n_coup) ! writing the gap info in the file
      end if
    case default ! case 4: linear structure
        write(*,*)
        write(*,'(1X,A,F12.6)') "THE VALUE OF THE HOMO LUMO GAP IS: ", eigen(n_coup+1)-eigen(n_coup) ! computing the gap for case 4
        write(filename, '(A,F4.2)') 'gap_lin_'//trim(typ1)//"_"//trim(typ2)//"_", ratio ! new namefile before opening the unit
        open(unit=15, file=filename, status='unknown', access='append') ! opening the unit (file) where the gap will be stored
        write(15,'(I4.4,F12.6)') d, eigen(n_coup+1) - eigen(n_coup) ! writing
  end select
  close(15) ! closing unit(15) regardles which case was considered
  write(*,'(1X,A,A)') "THE HOMO-LUMO GAP VALUE WAS APPEND IN THE FILE: ", filename ! displaying to the user where the gap info was
                                                                                   ! stored
end subroutine HL_gap
