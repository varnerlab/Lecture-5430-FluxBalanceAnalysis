# --- PRIVATE METHODS BELOW HERE ------------------------------------------------------------------------------- #
"""
Private default handler that parses a raw JSON response string from the BiGG API.
Dispatches to the appropriate JSON.parse call based on the endpoint model type.
Returns a parsed Julia object (typically a `Dict`), or `nothing` if the model type is unrecognized.
"""
function _default_handler_process_bigg_response(model::Type{T},
    response::String) where T <: AbstractBiggEndpointModel

    # initialize -
    type_handler_dict = Dict{Any,Function}()

    # hardcode the response handler -
    type_handler_dict[MyBiggModelsEndpointModel] = (x::String) -> JSON.parse(x) # default handler
    type_handler_dict[MyBiggModelsDownloadModelEndpointModel] = (x::String) -> JSON.parse(x) # default handler

    # lookup the function to handle the response -
    if (haskey(type_handler_dict, model) == true)
        handler_function = type_handler_dict[model]
        return handler_function(response);
    end

    # default: return nothing
    return nothing
end
# --- PRIVATE METHODS ABOVE HERE ------------------------------------------------------------------------------- #

# --- PUBLIC METHODS BELOW HERE -------------------------------------------------------------------------------- #
"""
    process_forecast_response_dataframe(model::Type{T}, response::String) -> DataFrame where T <: MyBiggModelsEndpointModel

Process the JSON response string from a BiGG API endpoint and return a `DataFrame`.

### Arguments
- `model::Type{T}`: the endpoint model type. Must be a subtype of `MyBiggModelsEndpointModel`.
- `response::String`: the raw JSON response string from the BiGG API.

### Returns
- `DataFrame`: a `DataFrame` populated with the response data.
"""
function process_forecast_response_dataframe(model::Type{T}, response::String)::DataFrame where T <: MyBiggModelsEndpointModel

    # initialize -
    dataframe = DataFrame();

    # parse the response -
    response_dict = JSON.parse(response);
    for data ∈ response_dict["properties"]["periods"]
        
        # put data into a row. This is an example of a NamedTuple
        # Notice: we are selecting the keys we want from the JSON response, and returning them in the order we want.
        # This is a form of data transformation and cleaning. We only give the user the data they need.
        row = (
            startTime = data["startTime"],
            endTime = data["endTime"],
            isDayTime = data["isDaytime"],
            temperature = data["temperature"],
            temperatureUnit = data["temperatureUnit"],
            windSpeed = data["windSpeed"],
            windDirection = data["windDirection"],
            shortForecast = data["shortForecast"]
        );
        push!(dataframe, row);
    end
 
    # return the dataframe -
    return dataframe
end
# --- PUBLIC METHODS ABOVE HERE -------------------------------------------------------------------------------- #