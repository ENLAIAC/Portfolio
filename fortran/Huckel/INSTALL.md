# Installation and Usage — Hückel Polyene Fortran Program

This document explains how to compile and run the **Hückel/tight-binding polyene** code (`Huckel_mat.f90` + modules), including required dependencies and expected outputs.

## 1. Requirements

### Compilers
- A modern Fortran compiler (tested setup uses **gfortran**).

### Linear algebra library
- **LAPACK** is required (the program links with `-llapack` in the provided `makefile`).

> Note: on some systems LAPACK may require BLAS as well (often pulled automatically by the package manager). If you see linker errors about BLAS symbols, link with `-lblas` or use a combined implementation like OpenBLAS.

## 2. Install dependencies

### Ubuntu / Debian
```bash
sudo apt-get update
sudo apt-get install -y gfortran make liblapack-dev
```

### Fedora
```bash
sudo dnf install -y gcc-gfortran make lapack-devel
```

### macOS (Homebrew)
```bash
brew install gcc make lapack
```

If the linker cannot find LAPACK on macOS, you may need:
```bash
brew info lapack
```
and then add the suggested `-L...` and `-I...` paths (or link against OpenBLAS).

### HPC clusters (modules)
Typical example (names vary by cluster):
```bash
module load gcc
module load lapack
```

## 3. Build

From the directory containing the sources and `makefile`:
```bash
make
```

This produces the executable:
```bash
./Huckel_mat
```

To remove compiled artifacts:
```bash
make clean
```

## 4. Run (interactive)

Run:
```bash
./Huckel_mat
```

The program will prompt for:
- system type: `L` (linear) or `C` (cyclic)
- chain length `d` (even values enforced by the current driver)
- bond alternation flag (if enabled, you provide a second coupling)
- site alternation flag (if enabled, you provide two on-site energies)

## 5. Run (non-interactive / batch)

You can feed inputs via stdin. Example (values are illustrative; follow the prompt order shown by the program):

```bash
printf "L\n100\nF\nF\n0.0\n" | ./Huckel_mat
```

For parameter scans, it is recommended to drive the executable from a shell/Python script and collect outputs.

## 6. Outputs (what files to expect)

The program writes a mix of “single-run” files and “scan-friendly” appended datasets.

### Top-level outputs (current working directory)
- `Huckel_matrix.txt`  
  The Hamiltonian matrix printed row-wise.
- `eigenvalues.txt`  
  Eigenvalue list (with a normalized index column).

### Organized outputs (auto-created folders)
The routine `save_res.f90` creates directories under `eigenvectors/`:

- `eigenvectors/linear/<alpha1>_<alpha2>/<d>_<ratio>/eigenvalues.txt`
- `eigenvectors/linear/<alpha1>_<alpha2>/<d>_<ratio>/eigenvector_i`
- and analogous paths for `cyclic/`

Folders are created automatically via `mkdir -p`, so you do not need to pre-create them.

### Dataset-style outputs (appended lines)
- `gap_lin_...` / `gap_cyc_...` (written by `HL_gap.f90`)  
  Each run appends a line: `d   gap`
- `TPS_lin_...` (written by `Huckel_mat.f90`, linear only)  
  Each run appends a line: `d   TPS_norm` where the code currently writes `lambda/d`.

## 7. Troubleshooting

### Linker errors: “undefined reference to dsy ev / lapack symbols”
Install the development LAPACK package and rebuild. If needed, modify `LDFLAGS` in the makefile to include BLAS:
```make
LDFLAGS = -llapack -lblas
```

### “cannot find -llapack”
Your LAPACK is not in the default library search path. Install `liblapack-dev` (Linux) or use Homebrew paths on macOS, then add:
- `-L/path/to/lapack/lib` to `LDFLAGS`.

### Output directories not created
The code calls `execute_command_line('mkdir -p ...')`. If you run on a restricted environment, ensure:
- you have permission to create directories in the working directory,
- `execute_command_line` is supported by your compiler (gfortran supports it).

## 8. Reproducibility recommendation (optional)
For scan campaigns, keep:
- a copy of the exact executable used,
- a run log listing parameter grids,
- a single consolidated `summary.tsv` (if you later add it) to simplify plotting.

