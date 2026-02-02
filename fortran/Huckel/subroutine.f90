!The next subroutin is designed to filll the Huckel Matrix. The idea of designing a routine to perform this action 
!comes from the requirements of the program to match multiple conditions, i.e., to include dimerization, both in the
!form of bond alternation and atomic centers alternation (or even the both of them)
!Using a subroutine a overload of the memory is avoided because only the dimension of the matrix, the matrix itself and
!the type (cyclic or linear) are global variable, while the other ones are local. Since the bond alternation or atom
!alternation information are useful to fill the matrix, but they're value is explicitely expressed by the matrix form
!and shape, is not required to initialize them in the main program, rather they're used as local variable here and then
!their memory location gets deallocate when the subroutine concludes. In the case of this subroutine there isn't a 
!proper gain by no calling it, because the routine is required to be called to run the program, rather a slightly
!memory gain is accomplished by separating the procedure. 't' as the type variable is not required to be a global 
!variable, but do to the fact that cyclic polyenes must fulfill the "even number of centers condition" is required to
!declare the tyope before and then assess the dimension. Since the dimension initialization as well as the allocation
!of the matrix memory occurs in the main program, the definition of the type was required to be there to. The decision
!to declare H (huckel matrix) and d (its dimension), as globl variables derives from the fact that they are used in 
!almost any step of the program, it would be useful to locally allocate H memory to just fill the matrix and print it
!without any global function storing the result of the subroutine. About the dimension, the eigenvalues vector and the
!TPF both runs over loops that uses d, hence eclaring d in the filling subroutin would not generate any global variable
!able to be reused to control further step of the program, hence d and H have been declared as global variables

subroutine matrix_fill (H, d, t)
  !INPUT VARIABLES
  integer, intent(in) :: d                  ! d= the dimension of the matrix
  character, intent(in) :: t                !t = 'C' if the polyene is cyclic, 'L' if linear

  !LOCAL VARIABLES
  logical :: bond_alt, atom_alt !bond_alt = set to true to design a polyene whose bonds lenght is alternated
                                !atom_alt = set to true to designa polyene whose atoms are alternated
  integer :: at1, at2         !at* = stores the atom type eventially atom_alt is set to true
  double precision :: beta1=-1.0d0, beta2=-1.0d0 !beta* = bond lenghts. beta1 is already assigned, 
                                                 !beta2 is assigned by the user when bond_alt is True
  integer :: i,j !counters to run through the matrix (reading, writing, filling)                                             
  !OPUTPUT VARIABLES
  double precision, intent(inout) :: H(d,d) !H(:,:) = the Huckel Matrix 


  !Defining whether the polyene consist of a chain of dymer and which is the lenght of the second bond in the dymer 
  !unit. The following branch of code ask to the user whether or not the polyene has bond with different lenghts. 
  !If so, the program asks the user to provide the lenght of the second bond, while the first one is kept to -1 by 
  !default
  write(*,'(A)')
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

  !Defining if the dymer contains the same atomic center. As in the previous step the programs ask whether the polyene
  !has dymers as basic units with alternating atomic centers. If so, the programs keep asking the atomic type. 
  !By contrast to the previous branch of code, here, no assumptions are made about the default atomic type. When no
  !alternation is observed, then at1 and at2 are set to be equal, otherwise they're set to different values.
  !In both cases, the user decides the nature of the atomic center of the polyene
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
  read(*,*) atom_alt !Reading wheter the alternation is observed or not. If yes is stated by 'T', otherwise by 'F'
  write(*,*)
  if ( atom_alt .eqv. .False. ) then;
    write(*,'(A)') "----------------------------         ATOM TYPES         -------------------------------"
    write(*,*)
    write(*,'(A)') "          CARBON = 0           BORON = 1           NITROGEN = 2           OXYGEN = 3       "
    write(*,*)
    write(*,'(A)', advance='no') "   ENTER THE ATOM TYPE: "
    read(*,*) at1 !Reading the atomic type
    at2=at1 !Given the no alternation condition, at2 is set equal to at1
  endif
  if ( atom_alt .eqv. .TRUE. ) then;
    write(*,'(A)') "----------------------------         ATOM TYPES         -------------------------------"
    write(*,*)
    write(*,'(A)') "          CARBON = 0           BORON = 1           NITROGEN = 2           OXYGEN = 3       "
    write(*,*)
    write(*,'(A)', advance='no') "   ENTER THE ATOM TYPE FOF BOTH ATOMS :  "
    read(*,*) at1, at2 !Reading two different atom type, at1 and at2 are different in this case
    if ( at1 .eq. at2 ) then;
      write(*,'(A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,'(A)') "!!!!!!!!!!!!!!!!             THE SAME ATOM TYPE WAS ENTERED               !!!!!!!!!!!!!!!"
      write(*,'(A)') "!!!!!!!!!!!!!!!!   the dymer has the same atomic center on each monomer   !!!!!!!!!!!!!!!"
      write(*,'(A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*)
    end if
  end if
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FILLING PROCEDURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Here the main challenge is to properly fill the Huckel matrix without exceeding its dimension. To avoid any memory 
  ! error the filling procedure is split in 2 steps, the filling of the diagonal and then the assignment of 
  ! off-diagonal terms

  ! Diagonal filling:
  ! The reasoning here is pretty simple, element on the diagonal alternates (if atom_alt is true), otherwise they have
  ! the same value. To use a single do to fulfill both the alternating and not alternating case the do is designed as
  ! the atom types were different, if they're not then at1=at2 and the conditions does not affect the filling cause
  ! the atomic type are identical, otherwise the same conditions allows to alternate the center in the case of atom_alt
  ! being true
  H(:,:)=0.0d0
  do i=1, d
    if ( mod(i,2) .eq. 1 ) then ;
      H(i,i)=at1
    else
      H(i,i)=at2
    end if
  end do

  ! Off-diagonal elements filling:
  ! The reasoning here is based on the observation than any rows refers to an atom, hence in the matrix a linking is
  ! repeated two times. As an example, if atom 2 and 3 are linked, the H(2,3) and the H(3,2) element will have the same
  ! value. The loop runs over d-1 atoms. Indeed, for a cyclic polyene the (1,d) and (d,1) elements are filled outside
  ! of th filling loop if the ( t .eq. 'C' ) condition is matched, if is not then the matrix is linear and such 
  ! elements are supposed to be zero. The filling does not fills the (i-1) elements,  but rather double fills (i,i+1) 
  ! and its simmetric elements (i+1,i). Hence is not even required to get to the last row, cause the last center link
  ! is already done by the previous step at i=d-1 when the filled elements are (d-1,d) and (d,d-1). If the polyene is
  ! cyclic then the filling is accomplished as final step.
    do i=1, d-1
      if ( mod(i,2) .eq. 1 ) then;
        H(i,i+1)=beta2
        H(i+1,i)=beta2
      else
        H(i,i+1)=beta1
        H(i+1,i)=beta1
      end if
    end do

    ! Filling of (1,d) and (d,1) when the cyclic condition is met, if t .eq. 'L' the condition is not satisfid and the 
    ! program neglects the command in the if

    if ( t .eq. 'C' ) then ;
            H(1,d)=beta1
            H(d,1)=beta1
    end if
    
  end subroutine matrix_fill
