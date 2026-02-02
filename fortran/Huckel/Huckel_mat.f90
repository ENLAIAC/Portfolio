program Huckel
  ! INITIALIZIATION OF THE VARIABLES
  ! HL = Huckel matrix for linear polyenes; HC = Huckel matrix for cyclic polyenes, eigen = vector storing the
  ! eigenvectors; delta(:) = stores the differences in energy between successive eigenvalues in order to find the HOMO-LUMO gap
  ! beta = value of the beta distance; alpha = value of diagonal term; lambda = [...], mu = [...]
  ! i,j,k = counters to run through loops, d = lenght of the polyene, uw = unit to control the printing of
  ! eigenvectors; n_el = [...used later for TPS...]; filename = stores the name of the file that will host the
  ! eigenvalues
  implicit none
  double precision, ALLOCATABLE :: H(:,:), eigen(:), occupation(:)
  double precision :: beta1, beta2, alpha, lambda, mu
  integer  :: i,j,k,d, uw, n_el
  character(len=4) :: t
  character(len=256) :: filename, folder, fullpath

  write(*,*)
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "                             ___         ___   ___   ___   ___   ___   ___       "
  write(*,'(1X,A)') "!!!!   |   | |   | |    | / |    |      |   | |   | |   | |   | |   | |   | |\  /|   !!!!!!"
  write(*,'(1X,A)') "!!!!   |___| |   | |    |/  |__  |      |___| |___| |   | |     |___| |___| | \/ |   !!!!!!"
  write(*,'(1X,A)') "!!!!   |   | |   | |    |\  |    |      |     |\    |   | |  __ |\    |   | |    |   !!!!!!"
  write(*,'(1X,A)') "!!!!   |   | |___| |___ | \ |___ |___   |     | \   |___| |___| | \   |   | |    |   !!!!!!"
  write(*,'(1X,A)') "                                                                                             "
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,*)
  write(*,'(1X,A)') "                                         by CAZZANTI ELIA                                     "
  write(*,'(A)')
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!           BUILDING OF THE HUCKEL MATRIX           !!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(A)')

  ! Definition of the type of polyene by user with input error handling: if the user enters an unexpected entry
  ! the program doesn't stop the program, but arise an explicit error and keeps asking for a correct input

  do
    write(*,'(1X,A)',advance='no') "ENTER THE TYPE OF POLYENE:  'L' = LINEAR   |   'C' = CYCLIC  ->  "
    read(*,*) t
    if ( t .eq. 'C' .or. t .eq. 'L') then;
       exit
    endif
    write(*,*)
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!               INPUT ERROR: PLEASE PROVIDE A CORRECT INPUT              !!!!!"
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)
  end do

  ! Definition of the size of the polyene chain. For a cyclic chain the number of atoms must be even, the next branch of
  ! of code handles any odd entry for cyclic chains. The program keeps asking for a correct length until the users entry
  ! satisfy the even number condition
  do
    write(*,'(1X,A)',advance='no') "ENTER THE DIMENSION OF THE POLYENE (N):  "
    read(*,*) d
    if ( t .eq. 'C' .and. mod(d,2) .eq. 0) then ;
            exit;
    else if ( t .eq. 'L') then ;
            exit;
    endif

    write(*,*)
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,'(1X,A)') "!!!!!!!!                 SIZE ERROR: PLEASE PROVIDE A CORRECT INPUT                   !!!!!"
    write(*,'(1X,A)') "!!!!!!!!        NOTE: FOR A CYCLIC POLYENE THE DIMENSION MUST BE EVEN                 !!!!!"
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)
  end do


  ! Memory allocation of Huckel matrix (dxd) and eigenvalue vector (d) 
  allocate(H(d,d), eigen(d))

  ! 'Matrix_fill' ask to the user some additional features of the system like if dymers alternate bond lenght or the atomic center
  ! type and fills the Huckel's matrix accordingly. Look at 'subroutine.f90' for more details.
  call matrix_fill(H, d, t, beta1, beta2)

  !Open a file that is 'replaced' any time the program is run and stores the Huckel's matrix generated  by 'Matrix_fill'. The
  !writing unit is open in the main program, then the file is filled with a for running through the Huckel's matrix and the unit is
  !closed.
  open(unit=10, file="Huckel_matrix", status='replace')
  do i=1, d
    write(10,*) (H(i,j), j=1, d)
  end do

  close(10)

  ! DIAGONALIZATON  PHASE:
  ! Using a externl subroutine that exploits the fortran buil-in function dsyev, the Huckel Matrix is diagonalized and overwritten,
  ! hence the matrix returned by the subroutine does not mantain anymore its initali form, but rather contains the d (dimension)
  ! eigenvectors. In this part of code the eigenvalues will be displayed in files, one for each eigenvector, and printed in the
  ! working folder. Each file shows a reference to the eigenvectors it refers to. This eigenvectors are the columns of the
  ! diagonlized matris. A second printing procedure is performed over the eigenvalue vector: inside a file 'eigenvalues.txt' the
  ! eigenvalues are printed following their normalized index. This way of storing them helps the representation and allows to
  ! restrict the x-axis range from 0 to 1. The x value represent the principal quantum number of the eigenstate associated to the
  ! eigenvalue, normalizing n over the largest n allows to display the eigenvalue energy in the [0,1] interval.

  call diagonalize_matrix(d, H, eigen) !Diagonalization

  open (unit=11, file='eigenvalues.txt', status='unknown') !Eigenvalues unit opening

  do i = 1, d
    write(11,*) (i-1)/dble(d-1), eigen(i) !writing
  end do
  
  n_el=d/2 ! given 'n_el' is an integer, diving by 2 a odd number floors the result to the lower number, for even number the ratio
           ! is exact

  write(*,*)
  write(*,'(1X,A,F12.6)') "THE VALUE OF THE HOMO LUMO GAP IS: ", abs(eigen(n_el)-eigen(n_el+1))
  close(11) !Closing eigenvalue unit

  ! The following procedure may be tricky. Each iteration a different value of uw is assigned (where uw is the writing unit). In
  ! addition the content of filename is updated according to the changing 'i'. Hence at the firs iteration filename is
  ! "eigenvalue_1", at the second iteration becomes "eigenvalue_2", and so on up to d (dimension of the matrix). In each file is
  ! stored a column of the diagonalized matrix (i.e. a eigenvector), running over i shifts the column of reference by one each
  ! iteration.

  write(*,*)
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!                    WAIT: WRITING THE EIGENVECTORS...                        !!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,*)
  
  write(folder,'(I0)') d
  call execute_command_line("mkdir -p eigenvectors/"//trim(folder))

  do i = 1 , d
    uw=i*20 
    write(filename, '(A,I0)') 'eigenvector_',i
    fullpath = "eigenvectors/"//trim(folder)//"/"//trim(filename)
    open(unit=uw, file=fullpath, status='unknown')
    do j=1, d
      write(uw,'(F10.5)') H(j,i)
    end do
    close(uw)
  end do

  fullpath = "eigenvectors/"//trim(folder)//"/"
  write(*,'(1X,A,A)') "THE EIGENVECTORS HAVE BEEN COPIED IN THE FOLDER: ", trim(fullpath)
  ! TPS Calculation
  ! The TPS here will be calculated using spatial orbitals (i.e., neglecting spin) and evetually exploitng occupation to keep track
  ! of the double counting. To keep trck of closed and open shell system simultaneously the number of electrons is computed as it
  ! folows. As a second step a occupation array is filled to keep track of doubly and singly occupied orbitals. For a even number of
  ! electrons, the orbitals obtained by the diagonalization (spatial orbitals) are doubly occupied, otherwise when d is odd the last
  ! orbital will be singly occupied (i.e., provides a single contribution to the TPS). 
  ! Three nested loops are used to run through the TPS calculation formula.
  ! According to the provided formula, (see the report), the computation of the TPS follows by multypling term-wise the coefficients
  ! of AO used to expand the MOs obtained by the diagonalizaion of the Huckel matrix. Here, some assumptions have been done as the
  ! orthonormalization of the basis functions (AO) and the diagonal form of the position operator matrix, as well as the fact that
  ! the position are not centered (the center of the chain on the zero).
  !

  n_el=(d+1)/2 ! The formula to compute the number of electron follows from the requirement of disrtinguishing between even and odd
               ! number of electrons. The method holds on the assumption that any atom contribute with a single electron regardless
               ! the atom type. Hence for a even number of centers the number of electrons is the exact half and the system consists
               ! of a closed-shell system with d/2 doubly occupied orbitals. Given the integer kind of 'n_el' for a odd number of
               ! electrons, dividing by two returns an incorrect value (e.g. n_el=3 for d=7). To avoid this wrong counting the '+1'
               ! term is used. Doing so, for even number of electrons the division by two will still return the correct value as 
               ! integer value rounded to the floor, while for an odd amount of electrons the 'n_el' increase by one. Indeed, with
               ! an odd number of electrons the last spatial orbital is singly occupied, an occupation array keeps track of the
               ! occupation of the spatial orbitals. 
 
              
              
              
  allocate(occupation(n_el))
  occupation(:)=2.00d0
 
  if ( mod(d,2) .ne. 0 ) then;
    occupation(n_el)=1.0d0 
  end if

  lambda=0.0d0
  do i=1, n_el
    do j=n_el+1, d
      mu=0
      do k=1, d
        mu= mu + H(k,i)*H(k,j)*dble(k)
      end do
      mu=occupation(i)*(mu*mu)
      lambda = lambda + mu
    end do
  end do
  
  write(*,'(A)')
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!                 TPS PRINTING...                 !!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(A)')

  write(filename, '(A,F4.2)') 'TPS_lin_', abs(beta1/beta2)
  open(unit=15, file=filename, status='unknown', access='append')
  write(15,'(I4.4,F12.6)') d, lambda/dble(d) 
  close(15)
  write(*,'(1X,A,A)') "THE TPS VALUE WAS APPEND IN THE FILE: ", filename
  deallocate(H,eigen,occupation) !deallocating the memory reserved to H and eigen to avoid memory leakage
end program Huckel

