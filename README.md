# Lecture: Flux Balance Analysis (CHEME 5430/5450, Spring 2026)

This repository contains lecture and example notebooks on metabolic engineering and flux balance analysis (FBA) for CHEME 5430/5450. Notebooks are written in [Julia](https://julialang.org) and use [Jupyter](https://jupyter.org) via [IJulia](https://github.com/JuliaLang/IJulia.jl).

---

## Notebooks

### Lecture
| Notebook | Description |
|---|---|
| `CHEME-5430-Lecture-FluxBalanceAnalysis-Spring-2026.ipynb` | Introduction to metabolic engineering and FBA. Covers the stoichiometric matrix, the FBA linear program, and the general and simplified flux bounds model. |
| `CHEME-5430-Advanced-Derivation-FluxBalanceAnalysis-Spring-2026.ipynb` | Derives the FBA steady-state mass balance constraint from open species mole balances. Shows how the Palsson constraints follow from steady-state, constant-volume, and no-transport assumptions, and explains the role of exchange reactions. |

### Examples
| Notebook | Description |
|---|---|
| `CHEME-5430-Example-SVD-StoichiometricMatrix-Spring-2026.ipynb` | Downloads genome-scale stoichiometric matrices from [BiGG Models](http://bigg.ucsd.edu/) and performs singular value decomposition (SVD). Extracts conservation relations from the left null space ($\mathbf{U}_0$) and steady-state flux patterns from the right null space ($\mathbf{V}_0$). |
| `CHEME-5450-Example-Solution-UreaCycle-S2026.ipynb` | Problem Set 2 (PS2). Constructs and solves an FBA problem for the urea cycle in HL-60 cells. Estimates reaction reversibility using [eQuilibrator](https://equilibrator.weizmann.ac.il) and maximum reaction velocities from [BRENDA](https://www.brenda-enzymes.org/). Solves the LP using [GLPK](https://www.gnu.org/software/glpk/). |

---

## Setup

### Prerequisites
- [Julia v1.10+](https://julialang.org/downloads/)
- [IJulia](https://github.com/JuliaLang/IJulia.jl) for Jupyter notebook support

### Installation

1. Clone this repository:
    ```bash
    git clone https://github.com/varnerlab/Lecture-5430-FluxBalanceAnalysis.git
    cd Lecture-5430-FluxBalanceAnalysis
    ```

2. Start Julia and install IJulia (first time only):
    ```julia
    using Pkg
    Pkg.add("IJulia")
    ```

3. Launch Jupyter:
    ```julia
    using IJulia
    notebook()
    ```

4. Open any notebook. The `Include.jl` file will automatically install and load all required packages on first run.

---

## Dependencies
Key packages used across notebooks:

| Package | Purpose |
|---|---|
| `JuMP.jl` | Formulating the FBA linear program |
| `GLPK.jl` | Solving the LP |
| `DataFrames.jl` | Tabular data for flux results |
| `PrettyTables.jl` | Formatted table display |
| `LinearAlgebra` | SVD and matrix operations |
| `Plots.jl` | Visualization |
| `HTTP.jl` / `JSON.jl` | BiGG Models API access |
