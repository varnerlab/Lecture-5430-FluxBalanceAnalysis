# Lecture: Flux Balance Analysis (CHEME 5430/5450, Spring 2026)

**Repository:** [https://github.com/varnerlab/Lecture-5430-FluxBalanceAnalysis](https://github.com/varnerlab/Lecture-5430-FluxBalanceAnalysis.git)

This repository contains lecture and example notebooks on metabolic engineering and flux balance analysis (FBA) for CHEME 5430/5450. Notebooks are written in [Julia](https://julialang.org) and use [Jupyter](https://jupyter.org) via [IJulia](https://github.com/JuliaLang/IJulia.jl).

---

## Overview

Cells are chemical factories: they consume raw materials, run hundreds of simultaneous reactions, and produce everything from energy carriers to complex macromolecules. Flux Balance Analysis (FBA) is a mathematical framework for analyzing the flow of metabolites through these reaction networks. Rather than tracking every concentration over time, FBA asks a simpler question: given the network structure and bounds on reaction rates, what steady-state flux distribution optimizes a chosen objective (e.g., growth rate or product yield)?

These notebooks build the FBA framework from first principles. We start with the biochemistry of metabolic networks, encode that biochemistry as a stoichiometric matrix, and then formulate and solve the FBA linear program. Along the way, we examine genome-scale models from public databases and work through a concrete FBA calculation on the urea cycle.

---

## Suggested Order of Study

Work through the materials in this order:

1. **Start with the lecture notebook** — introduces metabolic engineering, the stoichiometric matrix, and the FBA linear program.
2. **Read the advanced derivation** — derives the steady-state mass balance constraint from first principles, so you understand *where* $\mathbf{S}\hat{\mathbf{v}} = \mathbf{0}$ comes from.
3. **Work through the SVD example** — explores the structure of real genome-scale stoichiometric matrices using singular value decomposition.
4. **Complete the urea cycle example** — constructs and solves a full FBA problem for a small but realistic metabolic network.

---

## Notebooks

### Lecture

| Notebook | Description |
|---|---|
| `CHEME-5430-Lecture-FluxBalanceAnalysis-Spring-2026.ipynb` | Introduction to metabolic engineering and FBA. Covers the stoichiometric matrix, the FBA linear program, and the general and simplified flux bounds model. |
| `CHEME-5430-Advanced-Derivation-FluxBalanceAnalysis-Spring-2026.ipynb` | Derives the FBA steady-state constraint $\mathbf{S}\hat{\mathbf{v}} = \mathbf{0}$ from open species mole balances. Shows how the Palsson constraints follow from steady-state, constant-volume, and no-transport assumptions, and explains the role of exchange reactions. |

### Examples

| Notebook | Description |
|---|---|
| `CHEME-5430-Example-SVD-StoichiometricMatrix-Spring-2026.ipynb` | Downloads genome-scale stoichiometric matrices from [BiGG Models](http://bigg.ucsd.edu/) and performs singular value decomposition (SVD). Extracts conservation relations from the left null space ($\mathbf{U}_0$) and steady-state flux patterns from the right null space ($\mathbf{V}_0$). |
| `CHEME-5450-Example-Solution-UreaCycle-S2026.ipynb` | Constructs and solves an FBA problem for the urea cycle in HL-60 cells. Estimates reaction reversibility using [eQuilibrator](https://equilibrator.weizmann.ac.il) and maximum reaction velocities from [BRENDA](https://www.brenda-enzymes.org/). Solves the LP using [GLPK](https://www.gnu.org/software/glpk/). |

### Static Markdown (Pre-executed)
The `md/` folder contains pre-executed, read-only snapshots of all notebooks with outputs populated. These are useful for quickly reviewing results without running Julia.

---

## Key Concepts

These notebooks cover the following core ideas:

- **Metabolic networks and the stoichiometric matrix $\mathbf{S}$**: Biochemical reactions are encoded as a matrix $\mathbf{S}\in\mathbb{R}^{|\mathcal{M}|\times|\mathcal{R}|}$, where rows are metabolites and columns are reactions. The sign of each entry encodes whether a metabolite is consumed ($\sigma_{ij}<0$), produced ($\sigma_{ij}>0$), or uninvolved ($\sigma_{ij}=0$).

- **The FBA linear program**: At steady state, $\mathbf{S}\hat{\mathbf{v}} = \mathbf{0}$. FBA finds the flux vector $\hat{\mathbf{v}}$ that satisfies this constraint and flux bounds $\mathcal{L}\leq\hat{\mathbf{v}}\leq\mathcal{U}$ while maximizing (or minimizing) a linear objective $\mathbf{c}^{\top}\hat{\mathbf{v}}$.

- **Flux bounds**: Each reaction flux is bounded by a general model involving maximum velocity $v_{\max}$, reversibility, allosteric regulation, and substrate saturation. A simplified bounds model is often used in practice.

- **Null space structure**: The null spaces of $\mathbf{S}$ reveal important biological structure. The right null space ($\mathbf{V}_0$) spans feasible steady-state flux distributions; the left null space ($\mathbf{U}_0$) encodes conservation relations (e.g., conserved moieties).

- **Exchange reactions**: Because real cells exchange material with their surroundings, FBA models include hypothetical exchange reactions that allow the steady-state constraint to represent an open system.

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

> **Note for first-time users:** Package installation on first run can take several minutes. This is normal — Julia compiles packages on first use. Subsequent runs will be faster.

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

---

## External Resources

| Resource | Description |
|---|---|
| [BiGG Models](http://bigg.ucsd.edu/) | Genome-scale metabolic models for hundreds of organisms |
| [KEGG](https://www.genome.jp/kegg/) | Curated metabolic pathway maps and reaction data |
| [eQuilibrator](https://equilibrator.weizmann.ac.il) | Thermodynamic data for estimating reaction reversibility |
| [BRENDA](https://www.brenda-enzymes.org/) | Enzyme kinetics data, including $v_{\max}$ values |
| [Orth et al. (2010)](https://pubmed.ncbi.nlm.nih.gov/20212490/) | Foundational FBA primer ("What is flux balance analysis?") |
| [Palsson Lab](https://systemsbiology.ucsd.edu/) | Source of the genome-scale modeling framework used throughout |
