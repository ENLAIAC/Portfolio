# Huckel program

# Goal:

---

The goal of the program is to provide a simple implementation of the Huckel approach for the description of polyenes. It allows to include dimerizations by alternating  bonds or atom types (or both simultaneously). The program itself builds the Huckel matrix from scratch, following the users preferences. It asks for the length of the chain, the type (linear or cyclic), if dimerization is required and which type, letting the user to decide for the parameters. Given the distribution of eigenvalues’ energy depends by the relative difference of bond length ($\Delta\beta$), a value of beta is kept fixed, while only the second one can be arbitrarly changed by the user. The values of alternating alpha can be both set by the user according the experimental parameters. 

The program diagonalizes the Huckel matrix created by the user and derive the HOMO-LUMO gap according to the distribution of the eigenstates. It exploits a subroutine that selectively computes such a gap depending on the polyene’s parameters such as its type (linear or cyclic) and the amount of electrons (i.e., the amount of carbons). The analyses of the HOMO-LUMO gap provides a first hint about the metallic/insulator behaviour of the molecule.

For linear polyenes the Total Position Spread tensor (TPS) is computed too. The value of the TPS assumes a scalar form for 1D polyenes, while its trends against the length of the chain expresses how such a property change as the dimension of the polyene increases. Both HOMO-LUMO gap and the TPS value (for linear molecules) are printed to the screen and explicitely stated each run. Furthermore, the program automatically handles the storage of such a data. For each run eigenvalues, eigenvectors, HOMO-LUMO gap values and TPS values are stored in personalized/folders (see here). The folder contains a second folder `Huckel_aut/` that contains an automatized version of the program used to generate plots and iteratively run the program (see here).

# Description:

---

The description of the program will pass by each  file that can be find in this folder. When the file is a fortran file, the describe its content and explain its usage and goal. For building files (e.g., *makefile*), postprocessing of the results or data files a brief explanation is provided.

## Main program (`Huckel_mat.f90`)

### Variables declaration

```fortran
! INITIALIZIATION OF THE VARIABLES
```

This part of the program set the Local variables used in the main program, any local variable exlusively used in the subroutines is not displayed here.

**Variables:**

- `H(:,:)` : Huckel matrix.
- `beta1, beta2` : $\beta_1,\;\beta_2$ in atomic units (Hartree), respectively.
- `mu` : updates lambda over the atomic orbitals iterations.
- `at1, at2` : $\alpha_1,\;\alpha_2$ in atomic units (Hartree), respectively.

- `eigen(:)`: Eigenvalues vector.
- `lambda` : non-normalized TPS.
- `n_el`: number of electrons
- `uw` : writing unit
- `i,j,k` : counters.
- `t` : chain type (linear of cyclic)
- `filename` : stores the name of the saving files

---

### Type and dimension settings (by the user)

```fortran
! GRAPHICAL INTERFACE
```

This section is devoted to provide a first graphical interface of the program. This is the first interactive part of the program and it aims to store the parameters to build the Huckel matrix, here type and the length of the chain are collected, while dimerization parameters are gathered in the *filling routine*. It initializes it and defines the program stating its title as well as the author. Secondly, it asks to the user for the preferred type of the polyene. The programs allows the user to enters only two input values: `'L'` for linear, `'C'` for cyclic. The program handles entering an unexpected alphabetical input through a `do loop`, when the `( t .eq. 'C' .or. t .eq. 'L’` ) is met, the loop ends and the programs continues to the next interacting part. When an unexpected value for the polyene type is entered an error message is provided stating the provided input is unexpected and asks for a correct input.

The second interacting part asks the user for the dimension of the polyene. The program expects a even number, an approach identical to the one for the type with an appropriate condition is used to handle unexpect odd dimensions. As for the type, the program returns an explanation of the source of error and guides the user to provide a correct value. If an alphabetical (or non-numerical) the programs stops and arises an error “*Bad integer for item 1 in list input*”. `d` is delcared as integers, assigning it a non-integer value generates an error.

---

### Memory allocation of variables

```fortran
! MEMORY ALLOCATION
```

The Huckel matrix and the eigenvalues vector dimensions are assigned. `Huckel(d,d)` is a $d\times d$ matrix, while `eigenv(d)` store $d$ eigenvalues. The matrix is filled as the next step by the *filling subroutine*, while `eigenv` is filled as a consequence of the diagonalization in the *diagonalization subroutine.*

---

### Calling the filling function & printing Huckel matrix

```fortran
! MATRIX FILLING
```

In this branch of code the *flling procedure* is called to build the Huckel matrix according to the parameters collected in the [**type and dimension settings section](https://www.notion.so/Huckel-program-302cc1b49e848045b18bd7aad96e671a?pvs=21).** Once the matrix is built and returned to the main code, a `do loop` is used to pass through its elements and print it line by line in a file called `Huckel_matrix.txt`. Since in fortran each time `write` is called the file line is incremented by one, to mantain the matrix format an implicit do is used inside an outer explicit loop:

```fortran
do i=1, d
    write(10,*) (H(i,j), j=1, d)
  end do
```

---

### Huckel matrix diagonalization: subroutine call

```fortran
! HUCKEL MATRIX DIAGONALIZATION
```

This branch of code only calls the *matrix diagonalization* subroutine

---

### HOMO-LUMO calculation: subroutine call

```fortran
! HOMO-LUMO GAP CALCULATION   
```

This branch of code only calls the subroutine implemented to compute the *HOMO-LUMO gap* case by case.

---

### Storing eigenvalues and eigenvectors: subroutine call

```fortran
! EIGENVALUES AND EIGENFUNCTIONS PRINTINGS
```

In this section the program calls a fucntion devoted to save and store eigenvalues and eigenvectors in the appropriate folder, whose name is defined depending on the parameters collected and the h