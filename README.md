# Fluor-FOS
Open-source software (Fluor-FOS) for fluorescent nanocomposite multi-layer optical modeling. The program can handle a variable number of layers, scattering particles, and fluorescent particles, all with distinct material and optical properties. This program is equipped with a user-friendly interface that allows users to insert particles and media, then arrange them in a multi-layer structure with adjustable structural parameters.

## Authors and References
+ Khalid Alhammadi: alhammak@purdue.edu
+ Daniel Carne: dcarne@purdue.edu
+ Xiulin Ruan: ruan@purdue.edu

Cite: *Fluor-FOS: Open-Source Code for Optical Modeling of Multilayer Nanocomposite Media with Fluorescent Inclusions* (In revision)
## Requirements and Installation

The software is written in Python 3.11 or later and is compatible with any operating system (Windows, macOS, Linux).

### Python Dependencies

The following packages are required:

| Package | Version | Purpose |
|---------|---------|---------|
| `numpy` | 1.23.5 | Array operations and numerical computing |
| `scipy` | 1.14.1 | Special functions and numerical integration |
| `numba` | 0.56.4 | JIT compilation for parallelized Monte Carlo |
| `matplotlib` | 3.7.1 | Spectral property visualization and color display |

### Option 1: pip (VS Code or any standard Python environment)

Open a terminal in your working directory and run:
```bash
pip install -r requirements.txt
```

This will install all packages listed in `requirements.txt`. If any package fails to install, please ensure you have the required system dependencies installed.

### Option 2: conda (Anaconda/Miniconda)

1. Create a new conda environment:
```bash
conda create --name fluorfos python=3.11
```

2. Activate the environment:
```bash
conda activate fluorfos
```

3. Install the required packages using pip inside the conda environment:
```bash
pip install -r requirements.txt
```

> **Note:** `pip` is used inside the conda environment rather than `conda install` because packages such as `numba` are more reliably installed via `pip` than through conda's default channels.

## Repository Structure

The following source files are provided in the repository. Upon executing `a_Main_fluor_gh_f.py`, all other modules are automatically called internally — **do not run the other Python files directly**.

| File | Description |
|------|-------------|
| `a_Main_fluor_gh_f.py` | Main script — entry point for execution |
| `a_MieTheory3fluor2_gh.py` | Mie theory calculations and effective medium corrections |
| `a_montefluor_gh_f.py` | Parallelized Monte Carlo photon transport algorithm |
| `Integration_fluor_gh.py` | Solar spectrum integration |
| `interpolatefluor_gh.py` | Spectral data interpolation |
| `color_post_gh.py` | Color post-processing and visualization |

All source files and material-property input files must be placed in the **same working directory**. Material-property input files include:
- Optical properties of matrix and filler materials
- Fluorescent emission profiles
- Quantum yield properties

---

## Running the Program

Execute the main script by running:
```bash
python a_Main_fluor_gh_f.py
```

Upon execution, the user is prompted to provide:
- The main input file (`input.txt`)
- Ensure that the optical‑property files for the matrix materials, filler materials, fluorescent emission profiles, and quantum yield are located in the same directory.

Example input files for three test cases and one uncertainty‑analysis case are provided in the ‘examples’ folder of the repository.

---

## Fluor-FOS Functionalities
The main functionalities our program offers can be summarized as follows:

1) **Mie theory calculations for spherical particle inclusions.** For non-spherical inclusions, users can input pre-calculated spectral attenuation properties (λ, μₛ, μₐ, and g) obtained from other numerical methods such as FEM (finite element method) and DDA (discrete dipole approximation). Other cases, such as dye inclusions, are discussed in the article.

2) **Multi-layer fluorescent particle modeling.** Users can insert multiple fluorescent particle inclusions in multi-layer structures with various incident light sources for creating new colors or spectral shifting purposes in 1D systems. Such strategies are relevant for radiative cooling and energy harvesting applications. The software has the potential to be expanded to 2D and 3D systems in future developments.

3) **Uncertainty quantification.** Users can incorporate uncertain optical properties by providing the standard deviation along with optical properties and specifying the desired number of runs for smooth intervals. The uncertainty analysis is based on a derived equation combining normal distribution with random number generation via the fundamental principles of Monte Carlo methods.

4) **Parallel processing capabilities.** The software utilizes parallel processing for photon launching using the Numba package, despite the computational challenges associated with 2D tensor calculations for reflectance, transmittance, and absorptance.
   
5) **Different types of absorptances.** The output provides detailed spectral and solar‑weighted values for fluorescent absorptance—quantified either as absorbed photons or radiative relaxation—as well as losses associated with the fluorescent process, including quantum‑yield limitations and Stokes‑shift energy loss. It also distinguishes between ordinary reflectance and fluorescent reflectance (i.e., fluorescence emitted toward the source side), both of which are reported and used in the plotted results.

6) **Substrate modeling.** Metallic or polymeric substrates can be included in the simulations.

## Output
The Fluor-FOS output for a given nanocomposite multi-layer medium is primarily directed toward thermal and visual performance:

1) **Spectral and solar‑weighted optical properties.** The spectral and solar-weighted radiosity (reflectance and fluorescence), transmittance, and different types of absorptance are provided, tunable with the wavelength range specified by the user (ensuring the wavelength range covers both excitation and emission).

2) **Color prediction.** The predicted color is calculated using CIE standards given the spectral radiosity and incident spectral power.


## Input File Format

Recommendations for Input File.txt Format:
- Iclude the material optical‑property text files next to the particle and medium keywords. The keywords before the ‘:’ must remain unchanged, as they are recognized by the program. An exception is the numbering of particle, matrix, quantum‑yield, and emission entries (e.g., particle 1:, particle 2:, matrix 1:, matrix 2:, qy1:, qy2:, emit1:, emit2:), which may be incremented or adjusted as needed. The number following emit and qy must correspond to the particle number of the fluorescent pigment (e.g., emit1 and qy1 refer to particle 1). All text appearing after the ‘:’ may be modified by the user. Here is an example input file:

 ```
  MC 
  output: aaa_case_uncertainty_results

  particle 1: bsr.txt  #
  particle 2: a_clas_red.txt
  particle 3: y2o3.txt
  matrix 1: acrylic.txt

  light: AM1555.txt

  emit1:bsr_emis.txt # 
  excit_start_end1: 325,700

  emit2:cals_emis.txt # emi2 is emission for particle 2
  excit_start_end2: 325, 700

  qy1: bsr_qy.txt  # 

  qy2: cals_qy.txt

  photons: 01
  Start: 310
  End: 1000
  Interval: 3

  Sim: 1
  Upper: 1
  Lower: 1
  Layer 1
  Matrix 1
  T: 200
  particle 1
  fluor: 1
  D: 5
  VF: 20
  std: 0
  particle 3
  fluor: 0
  D: 0.6
  VF: 10
  std: 0
```
- Another example when the user includes uncertainty analysis:

```
  MC
  output: aaa_case_uncertainty_results

  particle 1: bsr.txt  #  BSR.txt   a_green_phosphor   BSR_check_black  a_CLAS_red_check_black
  particle 2: a_clas_red.txt
  particle 3: y2o3.txt
  matrix 1: acrylic.txt

  light: AM1555.txt

  emit1:bsr_emis.txt #   emi1 is emission for particle 1
  excit_start_end1: 325,700

  emit2:cals_emis.txt # emi2 is emission for particle 2
  excit_start_end2: 325, 700

  qy1: bsr_qy.txt  # BSR_QY.txt ba2SiO4Eu2_QY.txt

  qy2: cals_qy.txt

  photons: 01
  Start: 310
  End: 1000
  Interval: 3
  Number_of_sims_uncertainty: 30

  Sim: 1
  Upper: 1
  Lower: 1
  Layer 1
  Matrix 1
  T: 200
  particle 1
  fluor: 1
  D: 5
  VF: 20
  std: 0
  particle 3
  fluor: 0
  D: 0.6
  VF: 10
  std: 0
```
For multiple simulations, an input file for such purpose can be made in this format: 

```
MC
output: aaa_case_valid_get_hub

particle 1:  bsr_increases.txt  # 
particle 2: a_clas_red.txt
particle 3: Y2O3.txt
matrix 1: acrylic.txt
matrix 2: silicone.txt

light: AM1555.txt

emit1: BSR_emis.txt # 
excit_start_end1: 312,900

emit2:cals_emis.txt # emi2 is emission for particle 2
excit_start_end2: 312, 900

qy1: cals_qy.txt  #

qy2: cals_qy.txt


photons: 10000
Start: 310
End: 1000
Interval: 5
Number_of_sims_uncertainty: 60


Sim: 1
Upper: 1
Lower: 1
Layer 1
Matrix 1
T: 500
particle 1
fluor: 1
D: 5
VF: 1
std: 0.0
particle 2
fluor: 1
D: 10
VF: 10
std: 0.0
particle 3
fluor: 0
D: 0.4
VF: 50
std: 0.0

Sim: 2
Upper: 1
Lower: 1
Layer 1
Matrix 1
T: 500
particle 1
fluor: 1
D: 5
VF: 10
std: 0.0
particle 2
fluor: 1
D: 10
VF: 10
std: 0.0
particle 3
fluor: 0
D: 0.4
VF: 50
std: 0.0

Sim: 3
Upper: 1
Lower: 1
Layer 1
Matrix 1
T: 500
particle 1
fluor: 1
D: 5
VF: 20
std: 0.0
particle 2
fluor: 1
D: 10
VF: 10
std: 0.0
particle 3
fluor: 0
D: 0.4
VF: 50
std: 0.0
```
