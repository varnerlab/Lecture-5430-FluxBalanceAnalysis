# Derivation: Flux Balance Analysis (FBA) Formulation
In this advanced lecture, we will derive the flux balance analysis (FBA) formulation. We'll start by introducing material balances, and then we'll use these balances to derive the FBA formulation, and show what assumptions we need to make to arrive at [the standard Palsson FBA formulation](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3108565/).

> __Learning Objectives:__
> 
> By the end of this derivation, you should be able to:
>
> * __Derive species concentration balances:__ Write dynamic species concentration balances for a system with reactions and transport streams, starting from open mole balance equations.
> * __Apply steady-state simplifications:__ Identify and apply the three assumptions (steady-state, constant volume, no physical transport) that reduce the general concentration balance to the standard FBA constraint.
> * __Explain exchange reactions:__ Describe how hypothetical exchange reactions allow a steady-state model to represent an open system that exchanges material with its surroundings.

Let's get started!
___

## Theory: Material Balances
Suppose we have [a system with abstract volume $V$](https://github.com/varnerlab/CHEME-5450-Lectures-Spring-2025/blob/main/lectures/week-5/L5b/docs/figs/Fig-System-Schematic.pdf), e.g., the physical volume inside a single cell, the physical volume in a test tube, or the mass of cells in a reactor, etc. Inside this system, we have reaction set $\mathcal{R}$, metabolite set $\mathcal{M}$, and a stream set $\mathcal{S}$ that connects the system and the surroundings. 

> __Material Balances:__
> 
> A material balance equation for species $i\in\mathcal{M}$ in a system with volume $V$ has four terms:
> $$
> \begin{equation}
> \text{Accumulation} = \text{Generation} + \text{Transport In} - \text{Transport Out}
> \end{equation}
> $$
> The terms of the material balance are defined as:
> * The __accumulation__ term is the rate of change of species $i$ in the system, the __generation__ term is the rate of production (consumption) of species $i$ by chemical reactions in the system, and the __transport__ terms describe the rate of physical (convection or passive diffusion) or logical transport of species $i$ into (from) the system.

Let's look at dynamic concentration balance equations.

## Dynamic species concentration balances
When describing systems with chemical reactions, we write reaction rate expressions in terms of concentration, e.g., mole per unit volume basis. The number of moles $n_{i}$ (units: `mol`) of species $i$ in a system is described by an _open species mole balance equation_:
$$
\begin{equation}
\sum_{s\in\mathcal{S}}\nu_{s}\dot{n}_{i,s} + \dot{n}_{G,i} = \frac{dn_{i}}{dt}
\qquad\forall{i}\in\mathcal{M}
\end{equation}
$$
However, we can re-write the number of moles of species $i$ as $n_{i} = C_{i}V$ for $i\in\mathcal{M}$
where $C_{i}$ is the concentration of species $i$ (units: `mole per volume`), and $V$ (units: `volume`) is the volume of the system. The species mole balance can be rewritten in concentration units as:
$$
\begin{equation}
\sum_{s\in\mathcal{S}}d_{s}C_{i,s}\dot{V}_{s} + \dot{C}_{G,i}V = \frac{d}{dt}\left(C_{i}V\right)\qquad\forall{i}\in\mathcal{M}
\end{equation}
$$
where $\dot{V}_{s}$ denotes the volumetric flow rate for stream $s$ (units: `volume/time`), $C_{i,s}$ denotes the concentration of species $i$ in stream $s$ (units: `concentration`), and $\dot{C}_{G,i}$ is the rate of generation of species $i$ by chemical reaction (units: `concentration/time`). 

The generation terms for species $i$ in the concentration balance can be written as:
$$
\begin{equation}
\dot{C}_{G,i}V = \sum_{j\in\mathcal{R}}\sigma_{ij}\hat{v}_{j}V\qquad\forall{i}\in\mathcal{M}
\end{equation}
$$
where $\sigma_{ij}$ denotes the stoichiometric coefficient of species $i$ in reaction $j$ (units: `dimensionless`), 
and $\hat{v}_j $ denotes the rate of the jth chemical reaction per unit volume (units: `concentration/volume-time`), and $V$ denotes the volume of the system (units: `volume`). Putting these ideas together gives the _species concentration balance_:
$$
\begin{equation}
\sum_{s\in\mathcal{S}}d_{s}C_{i,s}\dot{V}_{s} + \sum_{j\in\mathcal{R}}\sigma_{ij}\hat{v}_{j}V = \frac{d}{dt}\left(C_{i}V\right)\qquad\forall{i\in\mathcal{M}}
\end{equation}
$$
Note: the transport terms are shown as physical terms, but we could also write logical flow terms.

## Simplified Species Constraints
In the [FBA primer](https://pubmed.ncbi.nlm.nih.gov/20212490/) the material balance constraints were written as $\mathbf{S}\hat{\mathbf{v}} = 0$. This is a simplification; in reality, the material balance constraints are more complex. But let's see how we get there.
When doing FBA, we often will make three assumptions:
* __Steady-state, constant volume__: The biological system (or at least part of it) is in a steady state, and in whole-cell models, the volume of the culture is constant. However, this is not the case in fed-batch cultures; many industrial biotechnology processes operate in fed-batch mode. 
* __Specific units__: The volume basis for the _intracellular_ species concentrations is in specific units, i.e., per unit cell mass measured in grams dry weight (units: `gDW`). For cell-free systems, the volume basis is the volume of the reactor.
* __No transport or dilution terms__: We will assume that there are no physical transport or dilution terms in the material balance equations. This is a simplification, but it is often used in FBA. For example, this assumption is not the case for continuous cell-free systems.

### Palsson constraints
Let the volume of our system be written in specific units, i.e., $V=B\bar{V}$, where $B$ is the biomass concentration (units: `gDW/L`) and $\bar{V}$ is the volume of the culture (units: `L`). The material balance constraints can be simplified by assuming the species in our system are in a steady state. But there is more to this story. Let's expand the accumulation terms:
$$
\begin{align*}
\frac{d}{dt}\left(C_{i}V\right) &= \frac{d}{dt}\left(C_{i}B\bar{V}\right)\\
&= B\bar{V}\underbrace{\left(\frac{dC_{i}}{dt}\right)}_{\text{steady state}\,=\,0} + C_{i}B\underbrace{\left(\frac{d\bar{V}}{dt}\right)}_{\text{steady state\,=\,0}} + C_{i}\bar{V}\left(\frac{dB}{dt}\right)\\
C_{i}\bar{V}\left(\frac{dB}{dt}\right) & = \underbrace{\sum_{s\in\mathcal{S}}d_{s}C_{i,s}\dot{V}_{s}}_{\text{no transport\,=\,0}} + \sum_{j\in\mathcal{R}}\sigma_{ij}\hat{v}_{j}V\\
C_{i}\bar{V}\left(\frac{dB}{dt}\right) & = \sum_{j\in\mathcal{R}}\sigma_{ij}\hat{v}_{j}B\bar{V}\\
C_{i}\underbrace{\left[\frac{1}{B}\left(\frac{dB}{dt}\right)\right]}_{\text{specific growth rate $\mu$}} & = \sum_{j\in\mathcal{R}}\sigma_{ij}\hat{v}_{j}\\
C_{i}\mu & = \sum_{j\in\mathcal{R}}\sigma_{ij}\hat{v}_{j}\\
\sum_{j\in\mathcal{R}}\sigma_{ij}\hat{v}_{j} - \underbrace{C_{i}\mu}_{\text{small}\,\ll{1}} & = 0\\
\sum_{j\in\mathcal{R}}\sigma_{ij}\hat{v}_{j} & = 0\quad\forall{i\in\mathcal{M}}\quad\blacksquare
\end{align*}
$$

### Exchange reactions
Great! So then, does everything in our system have to be at a steady state? Not exactly.  
* We can think of the system we are studying as an open system, i.e., it can exchange material with the surroundings. Thus, while the components of the system are in a steady state, the universe (system + surroundings) as a whole is not.
* This is a subtle point, but it is essential to understand. The exchange of material with the surroundings is captured in the context of our three assumptions by writing _hypothetical reactions_ that exchange material with the surroundings. We call these reactions [the _exchange reactions_](https://github.com/varnerlab/CHEME-5450-Lectures-Spring-2025/blob/main/lectures/week-5/L5c/docs/figs/Fig-ExchangeReactions.png); this figure was reproduced from [Bordbar et al, 2014](https://pubmed.ncbi.nlm.nih.gov/24987116/).

___

## Summary
The standard FBA mass balance constraint follows from the dynamic species concentration balance after applying steady-state, constant-volume, and no-transport assumptions.

> __Key Takeaways:__
>
> * **Material balance structure:** Dynamic species concentration balances contain four terms — accumulation, generation by reaction, transport in, and transport out — connecting physical intuition to the mathematical FBA formulation.
> * **Palsson constraints:** The standard FBA mass balance constraint follows from assuming steady-state intracellular concentrations, constant system volume, and no physical transport terms, with growth dilution treated as negligible.
> * **Exchange reactions:** Exchange reactions capture material transfer with the surroundings, allowing a closed steady-state formulation to represent an open biological system.

These derivations connect the general material balance framework to the simplified FBA formulation used in practice.
___
