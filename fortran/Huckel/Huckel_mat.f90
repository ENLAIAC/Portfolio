program Huckel
  ! INITIALIZIATION OF THE VARIABLES
  implicit none
  double precision, ALLOCATABLE :: H(:,:), eigen(:)
  double precision :: beta1, beta2, lambda, mu, at1, at2
  integer  :: i,j,k,d, uw, n_el
  character(len=4) :: t
  character(len=256) :: filename

  ! GRAPHICAL INTERFACE
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
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!           BUILDING THE HUCKEL MATRIX            !!!!!!!!!!!!!!!!!!!"
  write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write(*,'(A)')

  do
    write(*,'(1X,A)',advance='no') "ENTER THE TYPE OF POLYENE:  'L' = LINEAR   |   'C' = CYCLIC  ->  "
    read(*,*) t
    if ( t .eq. 'C' .or. t .eq. 'L') then;
       exit
    endif
    write(*,*)
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!        UNEXPECTED INPUT: PLEASE PROVIDE A CORRECT INPUT         !!!!!!!!!!"
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)
  end do

  do
    write(*,'(1X,A)',advance='no') "ENTER THE DIMENSION OF THE POLYENE (N):  "
    read(*,*) d
    if ( mod(d,2) .eq. 0) then ;
            exit;
    endif
    write(*,*)
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,'(1X,A)') "!!!!!!!!                 SIZE ERROR: PLEASE PROVIDE A CORRECT INPUT                   !!!!!"
    write(*,'(1X,A)') "!!!!!!!!        NOTE: FOR A CYCLIC POLYENE THE DIMENSION MUST BE EVEN                 !!!!!"
    write(*,'(1X,A)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(*,*)
  end do


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

  ! HOMO-LUMO GAP CALCULATION

  call HL_gap(t,abs(beta1/beta2),d,eigen,at1,at2)

  ! EIGENVALUES AND EIGENFUNCTIONS PRINTINGS

  call save_res(H,eigen,d,t,abs(beta1/beta2),at1,at2)

  ! TPS Calculation
  
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
  deallocate(H,eigen) 
end program Huckel

