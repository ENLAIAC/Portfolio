subroutine matrix_fill (H, d, t, beta1, beta2,at1,at2)
  implicit none

  !INPUT VARIABLES
  integer, intent(in) :: d                  
  character, intent(in) :: t

  !LOCAL VARIABLES
  logical :: bond_alt, atom_alt !bond_alt = set to true to design a polyene whose bonds lenght is alternated
                                !atom_alt = set to true to designa polyene whose atoms are alternated
  integer :: i,j !counters to run through the matrix (reading, writing, filling)                                             

  !OPUTPUT VARIABLES
  double precision, intent(inout) :: H(d,d)  
  double precision, intent(out) :: beta1, beta2  
  double precision, intent(out) :: at1, at2      

  ! ASSIGNATION OF BOND LENGHTS                                                 
  beta1=-1.0d0
  beta2=-1.0d0

  ! BOND ALTERNATION CHOICE
  write(*,'(A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(A)') "!!!!!!!!!!         DEFINITION OF THE ATOM TYPES IN THE POLYENE CHAIN         !!!!!!!!!!!!"
  write(*,'(A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(A)') 
  write(*,'(A)') "INSTRUCTIONS:"
  write(*,'(A)') "- Enter TRUE if the bond lenght alternates across the polyene chain, enter FALSE otherwise. "
  write(*,'(A)') "- Any entry different from TRUE or FALSE will brake the code and the program crashes."
  write(*,'(A)')
  write(*,'(A)',advance='no') "   BOND ALTERNATION: "
  read(*,*) bond_alt
  if ( bond_alt .eqv. .True. ) then;
    write(*,'(A)',advance='no') "   ENTER THE LENGHT OF THE SECOND BOND: "
    read(*,*) beta2
  endif

  !ATOM ALTERNATIO CHOICE
  write(*,'(A)')
  write(*,'(A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(A)') "!!!!!!!!!         DEFINITION OF THE BOND LENGHTS IN THE POLYENE CHAIN         !!!!!!!!!!!"
  write(*,'(A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(A)')
  write(*,'(A)') "INSTRUCTIONS:"
  write(*,'(A)') "- Enter TRUE if the atom type alternates across the polyene chain, enter FALSE otherwise. "
  write(*,'(A)') "- Any entry different from TRUE or FALSE will brake the code and the program crashes."
  write(*,'(A)') 
  write(*,'(A)',advance='no') "   ATOM TYPE ALTERNATION: "
  read(*,*) atom_alt 
  write(*,*)
  write(*,'(A)',advance='no') "   ENTER THE VALUE OF THE FIRST ALPHA: "
  read(*,*) at1
  at2=at1
  if ( atom_alt .eqv. .True. ) then;
    write(*,'(A)',advance='no') "   ENTER THE VALUES OF THE SECOND ALPHA: "
    read(*,*) at2
  endif
  if ( at1 .eq. at2 ) then;
      write(*,'(A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,'(A)') "!!!!!!!!!!!!!!!!              THE SAME ATOM TYPE WAS USED                 !!!!!!!!!!!!!!!"
      write(*,'(A)') "!!!!!!!!!!!!!!!!   the dymer has the same atomic center on each monomer   !!!!!!!!!!!!!!!"
      write(*,'(A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*)
  end if
  

  ! FILLING PROCEDURE
  ! Diagonal filling
  H(:,:)=0.0d0
  do i=1, d
    if ( mod(i,2) .eq. 1 ) then ;
      H(i,i)=at1
    else
      H(i,i)=at2
    end if
  end do

  ! Off-diagonal elements filling:
  do i=1, d-1
    if ( mod(i,2) .eq. 1 ) then;
      H(i,i+1)=beta2
      H(i+1,i)=beta2
    else
      H(i,i+1)=beta1
      H(i+1,i)=beta1
    end if
  end do

    
  ! Periodic boundary conditions (cyclic case)
  if ( t .eq. 'C' ) then ;
          H(1,d)=beta1
          H(d,1)=beta1
  end if
    
  end subroutine matrix_fill
