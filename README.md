# Free Vibration of an Euler–Bernoulli Beam with Lumped Mass and Local Spring Support

**Analytical solution (transfer matrix method) vs. Finite Element validation**

Authors: Wang Zichuan, Lyu Tianxiang, Fang Sizhe  
Institution: School of Mechanical Engineering, Tianjin University

---

## Overview

This project investigates the free vibration of a uniform Euler–Bernoulli beam carrying:

- A **pin support** at an intermediate point (point B)
- A **lumped mass** at point C
- A **vertical linear spring** at point D
- **Free ends** at both extremities (A and E)

The beam is initially at rest, then an **initial transverse velocity** $v_0$ is applied to the lumped mass at C.

> **The central question**: how do local discrete elements (a support, a mass, a spring) alter the continuous beam's natural frequencies, mode shapes, and dynamic response?

### Structural Model

![Structural model](assets/image-20260616001638017.png)

The beam has total length $4l$, with equal spans $AB = BC = CD = DE = l$.

---

## Methodology

Two independent approaches are developed and cross-validated:

### Analytical: Transfer Matrix Method

1. **Governing equation**: Euler–Bernoulli beam theory with separation of variables
2. **Non-dimensionalization**: coordinate $\xi = x/l$, frequency parameter $\beta = \lambda l$
3. **Jump conditions**: shear force jumps at C (inertia of lumped mass) and D (spring reaction)
4. **State vector**: $\mathbf{z}(\xi) = [Y, Y', Y'', Y''']^\mathrm{T}$ unified across continuous segments and discrete components
5. **Characteristic equation**: $\det\mathbf{A}(\beta) = 0$ → eigenvalues $\beta_n$ → natural frequencies $\omega_n = \frac{\beta_n^2}{l^2}\sqrt{\frac{EJ}{\rho S}}$

### Numerical: Finite Element Method

- **Element**: 2-node Euler–Bernoulli Hermite beam element (2 DOFs per node: $w$, $\theta$)
- **Mesh**: 80 elements, 81 nodes
- **Constraint handling**: B pin by elimination, C mass added to mass matrix, D spring added to stiffness matrix
- **Eigenvalue problem**: $\mathbf{K}\boldsymbol{\Phi}_n = \omega_n^2 \mathbf{M}\boldsymbol{\Phi}_n$
- **Response**: modal superposition with 24 modes

### Parameters

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Mass ratio | $\alpha = m / (\rho S l)$ | 0.5 |
| Spring stiffness | $\kappa = k l^3 / (EJ)$ | 20.0 |

---

## Results

### Natural Frequencies

Excellent agreement for the first 11 modes:

| Mode | Analytical $\bar{\omega}_n$ | FEM $\bar{\omega}_n$ | Error (%) |
|:----:|:--------------------------:|:--------------------:|:---------:|
| 1 | 1.29291120 | 1.29291121 | $6.6\times10^{-7}$ |
| 2 | 1.95408261 | 1.95408264 | $1.6\times10^{-6}$ |
| 3 | 3.17948221 | 3.17948233 | $3.7\times10^{-6}$ |
| 4 | 6.93522035 | 6.93522170 | $2.0\times10^{-5}$ |
| ... | ... | ... | ... |
| 12 | 73.64864210 | 73.67590498 | $3.7\times10^{-2}$ |

> Mode 12 marks a turning point where spatial discretization error becomes non-negligible.

### Mode Shapes

The first four mode shapes show near-perfect overlap between analytical and FEM solutions, with $L_2$ relative errors at the $10^{-9}$ to $10^{-8}$ level.

### Dynamic Response

- **Displacement** at point C: RMS error $1.27\times10^{-3}$
- **Velocity** at point C: RMS error $1.01\times10^{-1}$ (more sensitive to high-order truncation)
- No phase drift observed over the full time window, confirming frequency consistency.

---

## Workflow

```
Fortran FEM solver → CSV export → Python post-processing → SciencePlots visualization
```

### Key output files

| File | Content |
|------|---------|
| `outputs/fem_frequencies.csv` | FEM natural frequencies |
| `outputs/fem_modes.csv` | FEM mode shapes |
| `outputs/fem_c_response.csv` | Dynamic response at point C |
| `outputs/frequency_comparison.csv` | Analytical vs. FEM comparison table |

### Visualization gallery

| Figure | Description |
|--------|-------------|
| `assets/q3_visualization/q3_characteristic_roots.png` | Characteristic determinant scan |
| `assets/q3_visualization/q3_mode_shapes.png` | Analytical mode shapes (modes 1–4) |
| `assets/q3_visualization/q3_c_point_response.png` | Analytical response at point C |
| `assets/q3_visualization/q3_spacetime_response.png` | Full-field spacetime response |
| `assets/q3_comparison/comparison_frequency_spectrum.png` | Frequency bar chart |
| `assets/q3_comparison/comparison_frequency_error.png` | Frequency relative error |
| `assets/q3_comparison/comparison_mode_overlay.png` | Mode shape overlay |
| `assets/q3_comparison/comparison_c_response.png` | C-point response comparison |

---

## Repository Structure

```
.
├── assets/                          # Figures and visualizations
│   ├── image-20260616001638017.png  # Structural model diagram
│   ├── q3_visualization/            # Analytical solution plots
│   └── q3_comparison/               # Analytical vs. FEM comparison plots
├── outputs/                         # CSV data from FEM solver
├── paper/                           # Report and documentation
├── slides/                          # Presentation slides (Beamer)
├── slides_new/                      # Revised presentation (24 pages)
├── fem_q3_solver.f90                # Fortran FEM solver
├── fem_q3_solver.exe                # Compiled solver
├── postprocess_fem_q3.py            # FEM data post-processing
├── compare_analytic_fem_q3.py       # Comparison and plotting
├── visualize_q3.py                  # Analytical visualization
├── AGENTS.md                        # Development workflow notes
└── 结构动力学q3.md                    # Problem description (Chinese)
```

---

## Dependencies

- **Fortran**: any modern Fortran compiler (tested with `gfortran` 13+ via MSYS2)
- **Python**: 3.10+
  - `numpy`, `scipy`
  - `matplotlib` + `SciencePlots`
  - `pandas`
- **LaTeX** (for slides): `xelatex` with `ctexbeamer`, `newtxmath`

---

## How to Reproduce

```bash
# 1. Run FEM solver
gfortran fem_q3_solver.f90 -o fem_q3_solver
./fem_q3_solver

# 2. Post-process FEM data
python postprocess_fem_q3.py

# 3. Generate comparison figures
python compare_analytic_fem_q3.py

# 4. Build presentation slides
cd slides_new && xelatex q3_slides.tex
```

---

## License

This project is open-sourced for academic reference and reproducibility.
