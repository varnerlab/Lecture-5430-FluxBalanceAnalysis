"""
Abstract base type for all BiGG database API endpoint models.
Subtypes represent specific API endpoints and carry any required parameters.
"""
abstract type AbstractBiggEndpointModel end

"""
Abstract base type for all flux calculation models.
Subtypes carry the data required to formulate and solve a flux balance problem.
"""
abstract type AbstractFluxCalculationModel end

"""
    MyBiggModelsEndpointModel <: AbstractBiggEndpointModel

Endpoint model for the BiGG `/api/v2/models` listing endpoint.
No parameters are required; instantiate with `MyBiggModelsEndpointModel()`.
"""
struct MyBiggModelsEndpointModel <: AbstractBiggEndpointModel

    # methods -
    MyBiggModelsEndpointModel() = new();
end

"""
    MyBiggModelsDownloadModelEndpointModel <: AbstractBiggEndpointModel

Endpoint model for downloading a specific BiGG model via `/api/v2/models/<bigg_id>/download`.

### Fields
- `bigg_id::String`: the BiGG model identifier (e.g., `"iJO1366"`).
"""
mutable struct MyBiggModelsDownloadModelEndpointModel <: AbstractBiggEndpointModel

    # data -
    bigg_id::String

    # methods -
    MyBiggModelsDownloadModelEndpointModel() = new();
end

"""
    MyPrimalFluxBalanceAnalysisCalculationModel <: AbstractFluxCalculationModel

Data model for the primal flux balance analysis (FBA) problem.
Populate all fields, then pass to `solve` to obtain the optimal flux distribution.

### Fields
- `S::Array{Float64,2}`: stoichiometric matrix (species × reactions).
- `fluxbounds::Array{Float64,2}`: flux bounds array (reactions × 2), where column 1 is the lower bound and column 2 is the upper bound.
- `objective::Array{Float64,1}`: objective function coefficients (length = number of reactions).
- `species::Array{String,1}`: species names/ids ordered to match the rows of `S`.
- `reactions::Array{String,1}`: reaction names/ids ordered to match the columns of `S`.
"""
mutable struct MyPrimalFluxBalanceAnalysisCalculationModel <: AbstractFluxCalculationModel

    # data -
    S::Array{Float64,2}; # stoichiometric matrix
    fluxbounds::Array{Float64,2}; # flux bounds
    objective::Array{Float64,1}; # objective function coefficients
    species::Array{String,1}; # species names/ids
    reactions::Array{String,1}; # reaction names/ids

    # methods -
    MyPrimalFluxBalanceAnalysisCalculationModel() = new();
end