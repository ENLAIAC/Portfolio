# Portfolio
Hi, I'm Elia, and I'm a Computational chemist. In this repository I'll store any project, program, or script I've done during my academic journey. The goal of this directory is to provide suggestions, help, and a free space to anyone who feels lost in the field of computational chemistry.
During my studies, I've always felt like I was missing something. I was an applied student and most of the time I preferred to sacrifice my free time to study and get good marks, and most of the time I accomplished it. Despite the good grades, I kept feeling behind, as everyone was proceeding, and I was stuck in the mud. Later, talking with colleagues and professors, I realized that it wasn't me, but everyone felt the same. And it wasn't my fault either, but the teaching system wasn't providing the right information for us to learn what requirements a computational chemist should fulfill. 
This repository won't replace your classes or any online courses provided by some IT guy, but rather will drive you through some of the most common tasks in computational chemistry. About that, what is the most efficient way to explain something about coding? Most of the time, the answer is coding, so let's delve into this and let's see what we get.

## Repository structure

| Folder | Contents |
|--------|----------|
| `fortran/` | Hückel MO theory, Hartree-Fock implementation |
| `python/` | Data extraction and analysis scripts for HPC workflows |
| `bash/` | Utility and job submission scripts for HPC clusters |
| `C/` | C programs (in progress) |
| `ML/` | Machine learning projects (in progress) |

## Python scripts

- **`energy_grep.py`** — reads CP2K energy output files from replica directories, converts energies from a.u. to kcal/mol, and generates energy distribution histograms
- **`extract.py`** — extracts individual configurations from a multi-frame XYZ trajectory file with atom count validation and debug logging

## Bash scripts

- **`check_convergence.sh`** — batch SCF convergence checker for CP2K labeling jobs; classifies runs as converged/failed/requires-manual-check and writes summary reports
- **`check_data.sh`** — validates the presence and consistency of expected data files across labeling dataset directories
- **`check_distance.sh`** — scans candidate configurations for atomic overlaps below a distance threshold
- **`energy_chart.sh`** — consolidates total energies from multiple CP2K replica directories into a single file for plotting
- **`energy_grep.sh`** — extracts energy data from LAMMPS/DeePMD run logs across multiple simulation directories
- **`generate_xyz.sh`** — converts LAMMPS trajectories and other formats to XYZ
- **`job_lammps-deepmd.sh`** — SLURM submission wrapper for LAMMPS+DeePMD exploration runs
- **`python_submission.sh`** — SLURM wrapper for submitting Python scripts to the cluster queue
