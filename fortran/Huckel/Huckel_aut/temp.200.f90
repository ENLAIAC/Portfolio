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
  double precision :: beta1, beta2, lambda, mu, at1, at2
  integer  :: i,j,k,d, uw, n_el
  character(len=4) :: t
  character(len=256) :: filename, folder, fullpath, typ
  double precision, parameter :: thresh=0.95

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

  !do
  !  write(*,'(1X,A)',advance='no') "ENTER THE TYPE OF POLYENE:  'L' = LINEAR   |   'C' = CYCLIC  ->  "
  !  read(*,*) t
  !  if ( t .eq. 'C' .or. t .eq. 'L') then;
  !     exit
  !  endif
  !  write(*,*)
  !  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  !  write(*,'(1X,A)') "!!!!!!!!!!!!!!!               INPUT ERROR: PLEASE PROVIDE A CORRECT INPUT              !!!!!"
  !  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  !  write(*,*)
  !end do

  t='C'


  ! Definition of the size of the polyene chain. For a cyclic chain the number of atoms must be even, the next branch of
  ! of code handles any odd entry for cyclic chains. The program keeps asking for a correct length until the users entry
  ! satisfy the even number condition
  !do
  !  write(*,'(1X,A)',advance='no') "ENTER THE DIMENSION OF THE POLYENE (N):  "
  !  read(*,*) d
  !  if ( mod(d,2) .eq. 0) then ;
  !          exit;
  !  endif
  !  write(*,*)
  !  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  !  write(*,'(1X,A)') "!!!!!!!!                 SIZE ERROR: PLEASE PROVIDE A CORRECT INPUT                   !!!!!"
  !  write(*,'(1X,A)') "!!!!!!!!        NOTE: FOR A CYCLIC POLYENE THE DIMENSION MUST BE EVEN                 !!!!!"
  !  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  !  write(*,*)
  !end do

  d=200

  ! MEMORY ALLOCATION

  allocate(H(d,d), eigen(d))

  ! MATRIX FILLING
  
  call matrix_fill(H, d, t, beta1, beta2, at1, at2)

  open(unit=10, file="Huckel_matrix.txt", status='replace')
  do i=1, d
    write(10,*) (H(i,j), j=1, d)
  end do

  close(10)
  
  ! HUCKEL MATRIX DIAGONALIZATION

  call diagonalize_matrix(d, H, eigen)

  open (unit=11, file='eigenvalues.txt', status='unknown')

  do i = 1, d
    write(11,*) (i-1)/dble(d-1), eigen(i) !writing
  end do
  close(11) 

  call HL_gap(t,abs(beta1/beta2),d,eigen,at1,at2)

  ! EIGENVALUES AND EIGENFUNCTIONS PRINTINGS

  call save_res(H,eigen,d,t,abs(beta1/beta2),at1,at2)

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

  if ( t .eq. 'C' ) then ;
    write(*,*)
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,'(1X,A)') "!!!!!!        THE COMPUTATION OF THE TPS FOR CYCLIC POLIENES IS NOT REQUIRED         !!!!!"
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)
    stop 1
  end if

  write(*,'(A)')
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!               TPS CALCULATION...                !!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(A)')
  

  n_el=d/2
  lambda=0.0d0
  do i=1, n_el
    do j=n_el+1, d
      mu=0
      do k=1, d
        mu= mu + H(k,i)*H(k,j)*dble(k)
      end do
      lambda = lambda + mu*mu
    end do
  end do
  
  write(*,*)
  write(filename, '(A,F4.2)') 'TPS_lin_', abs(beta1/beta2)
  open(unit=16, file=filename, status='unknown', access='append')
  write(16,'(I4.4,F12.6)') d, lambda/dble(d) 
  close(16)
  
  write(*,'(A)')
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!                 TPS PRINTING...                 !!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(A)')

  write(*,'(1X,A,F12.6)') "   THE TPS VALUE IS: ", lambda/dble(d)
  write(*,*)
  write(*,'(1X,A,A)') "    THE TPS VALUE WAS APPEND IN THE FILE: ", filename
  deallocate(H,eigen) !deallocating the memory reserved to H and eigen to avoid memory leakage
end program Huckel

