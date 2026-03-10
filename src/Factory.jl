


# --- PUBLIC METHODS BELOW HERE -------------------------------------------------------------------------------- #
"""
    build(base::String, model::MyBiggModelsEndpointModel; apiversion::String = "v2") -> String

Build the URL string for the BiGG models listing endpoint.

### Arguments
- `base::String`: the base URL of the BiGG REST API (e.g., `"http://bigg.ucsd.edu"`).
- `model::MyBiggModelsEndpointModel`: the endpoint model instance.
- `apiversion::String`: the API version string (default `"v2"`).

### Returns
- `String`: the fully constructed URL for the models listing endpoint.
"""
function build(base::String, model::MyBiggModelsEndpointModel; apiversion::String = "v2")::String
    
    # TODO: implement this function, and remove the throw statement
    # throw(ArgumentError("build(base::String, model::MyWeatherGridPointEndpointModel) not implemented yet!"));

    # build the URL string -
    url_string = "$(base)/api/$(apiversion)/models";

    # return the URL string -
    return url_string;
end

"""
    build(base::String, model::MyBiggModelsDownloadModelEndpointModel; apiversion::String = "v2") -> String

Build the URL string for downloading a specific BiGG model.

### Arguments
- `base::String`: the base URL of the BiGG REST API (e.g., `"http://bigg.ucsd.edu"`).
- `model::MyBiggModelsDownloadModelEndpointModel`: the endpoint model instance. Must have `bigg_id` set.
- `apiversion::String`: the API version string (default `"v2"`).

### Returns
- `String`: the fully constructed URL for the model download endpoint.
"""
function build(base::String, model::MyBiggModelsDownloadModelEndpointModel; apiversion::String = "v2")::String

    # get data -
    bigg_id = model.bigg_id;

    # build the URL string -
    url_string = "$(base)/api/$(apiversion)/models/$(bigg_id)/download";

    # return the URL string -
    return url_string;
end

"""
    build(modeltype::Type{MyPrimalFluxBalanceAnalysisCalculationModel}, data::NamedTuple) -> MyPrimalFluxBalanceAnalysisCalculationModel

Construct a `MyPrimalFluxBalanceAnalysisCalculationModel` from a `NamedTuple` of data.

### Arguments
- `modeltype::Type{MyPrimalFluxBalanceAnalysisCalculationModel}`: the model type to construct.
- `data::NamedTuple`: a named tuple with fields `S`, `fluxbounds`, `objective`, `species`, and `reactions`.

### Returns
- `MyPrimalFluxBalanceAnalysisCalculationModel`: a fully populated FBA model ready to pass to `solve`.
"""
function build(modeltype::Type{MyPrimalFluxBalanceAnalysisCalculationModel},
    data::NamedTuple)::MyPrimalFluxBalanceAnalysisCalculationModel

    # get data -
    S = data.S;
    fluxbounds = data.fluxbounds;
    objective = data.objective;
    species = data.species;
    reactions = data.reactions;
 
    # build an empty model -
    model = modeltype();
 
    # add data to the model -
    model.S = S;
    model.fluxbounds = fluxbounds;
    model.objective = objective;
    model.species = species;
    model.reactions = reactions;
     
     # return the model -
     return model;
end
# --- PUBLIC METHODS ABOVE HERE -------------------------------------------------------------------------------- #