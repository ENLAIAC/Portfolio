:
HERE=$PWD
G90=gfortran
for j in 1.0 #1.01 1.25 2.00 3.00 100.0
do
  sed -e "s/beta2=-1.0d0/beta2=-$j/g" < subroutine.f90 > sub.$j.f90
  for i in 10 20 30 40 50 100 150 200 250 350 500 1000
  do
    sed -e "s/d=n/d=$i/g" < Huckel_mat.f90 > temp.$i.f90	
    $G90 temp.$i.f90 sub.$j.f90 diag.f90 save_res.f90 HL_gap.f90 -llapack -o Huckel_mat
    ./Huckel_mat
    mv eigenvalues.txt eigenvectors/$i/
    #rm temp.* sub.*
  done
done
