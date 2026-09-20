import JLD2

const RESULT_SCHEMA_VERSION = 1

_string_key_dict(values::AbstractDict) =
    Dict{String,Any}(string(key) => value for (key, value) in pairs(values))

"""Create a portable result containing only plain dictionaries and arrays."""
function result_record(;
    metadata::AbstractDict,
    cluster::AbstractDict,
    parameters::AbstractDict,
    sector::AbstractDict,
    solver::AbstractDict,
    spectrum::AbstractDict,
    observables::AbstractDict=Dict{String,Any}(),
    provenance::AbstractDict=Dict{String,Any}(),
)
    record = Dict{String,Any}(
        "schema_version" => RESULT_SCHEMA_VERSION,
        "metadata" => _string_key_dict(metadata),
        "cluster" => _string_key_dict(cluster),
        "parameters" => _string_key_dict(parameters),
        "sector" => _string_key_dict(sector),
        "solver" => _string_key_dict(solver),
        "spectrum" => _string_key_dict(spectrum),
        "observables" => _string_key_dict(observables),
        "provenance" => _string_key_dict(provenance),
    )
    validate_result(record)
    return record
end

const _REQUIRED_RESULT_KEYS = (
    "schema_version",
    "metadata",
    "cluster",
    "parameters",
    "sector",
    "solver",
    "spectrum",
    "observables",
    "provenance",
)

function validate_result(record::AbstractDict)
    missing = [key for key in _REQUIRED_RESULT_KEYS if !haskey(record, key)]
    isempty(missing) ||
        throw(ArgumentError("result is missing keys: $(join(missing, ", "))"))
    record["schema_version"] == RESULT_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported result schema"))
    for key in _REQUIRED_RESULT_KEYS[2:end]
        record[key] isa AbstractDict ||
            throw(ArgumentError("result field '$key' must be a dictionary"))
    end
    spectrum = record["spectrum"]
    energies = spectrum["energy"]
    levels = spectrum["level_index"]
    length(energies) == length(levels) ||
        throw(ArgumentError("spectrum arrays have inconsistent lengths"))
    if haskey(spectrum, "energy_minus_block_ground")
        length(spectrum["energy_minus_block_ground"]) == length(energies) ||
            throw(ArgumentError("block-gap array has inconsistent length"))
    end
    if haskey(spectrum, "energy_minus_ground")
        length(spectrum["energy_minus_ground"]) == length(energies) ||
            throw(ArgumentError("gap array has inconsistent length"))
    end
    return true
end

"""Atomically save one `(sector,U/t)` JLD2 record."""
function save_result(path::AbstractString, record::AbstractDict; overwrite::Bool=false)
    endswith(lowercase(path), ".jld2") ||
        throw(ArgumentError("result path must end in .jld2"))
    validate_result(record)
    destination = abspath(path)
    mkpath(dirname(destination))
    isfile(destination) && !overwrite &&
        throw(ArgumentError("result already exists: $destination"))
    temporary = tempname(dirname(destination)) * ".jld2"
    try
        JLD2.jldopen(temporary, "w") do file
            file["result"] = Dict{String,Any}(record)
        end
        mv(temporary, destination; force=overwrite)
    catch
        isfile(temporary) && rm(temporary)
        rethrow()
    end
    return destination
end

function load_result(path::AbstractString)
    record = JLD2.load(path, "result")
    validate_result(record)
    return record
end
