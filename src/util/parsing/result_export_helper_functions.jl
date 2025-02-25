"""
    get_outage_pras_resultspec(
        template::RATemplate,
    )

Get ResultSpec to be used by PRAS.assess() to get asset availability based on the 
template.

# Arguments

  - `template::RATemplate`: PRAS problem template

# Returns

  - PRAS ResultSpec
"""

function get_outage_pras_resultspec(template::RATemplate)
    # We have to change this once we define a formulation for Lines
    resultsspecs = PRASCore.Results.ResultSpec[]
    push!(resultsspecs, ShortfallSamples())

    for model in template.device_models
        result_spec = get_pras_resultspec(typeof(get_formulation(model)))
        if !(result_spec in resultsspecs)
            push!(resultsspecs, result_spec)
        end
    end

    return resultsspecs
end

"""
Get PRAS ResultSpec for RAFormulation
"""
function get_pras_resultspec(::Type{GeneratorPRAS})
    return GeneratorAvailability()
end

"""
Get PRAS ResultSpec for RAFormulation
"""
function get_pras_resultspec(::Type{B}) where {B <: StoragePRAS}
    return StorageAvailability()
end

"""
Get PRAS ResultSpec for RAFormulation
"""
function get_pras_resultspec(::Type{B}) where {B <: GeneratorStoragePRAS}
    return GeneratorStorageAvailability()
end

"""
Get PRAS ResultSpec for AbstractRAFormulation
"""
function get_pras_resultspec(::Type{B}) where {B <: AbstractRAFormulation}
    return error("PRAS ResultSpec not defined for $(B)")
end

"""
Get DeviceRAModel for PRAS AbstractAvailabilityResult
"""
function get_device_ramodel(
    ::Type{Union{Nothing, PRASCore.Results.GeneratorAvailabilityResult}},
)
    return (model=GeneratorPRAS, key=:generators)
end

"""
Get DeviceRAModel for PRAS AbstractAvailabilityResult
"""
function get_device_ramodel(
    ::Type{Union{Nothing, PRASCore.Results.StorageAvailabilityResult}},
)
    return (model=StoragePRAS, key=:storages)
end

"""
Get DeviceRAModel for PRAS AbstractAvailabilityResult
"""
function get_device_ramodel(
    ::Type{Union{Nothing, PRASCore.Results.GeneratorStorageAvailabilityResult}},
)
    return (model=GeneratorStoragePRAS, key=:generatorstorages)
end

"""
Get SPIOutageResult fname for PRAS ResultSpec
"""
function get_outage_result_fname(::Type{GeneratorAvailability})
    return :gen_availability
end

"""
Get SPIOutageResult fname for PRAS ResultSpec
"""
function get_outage_result_fname(::Type{StorageAvailability})
    return :stor_availability
end

"""
Get SPIOutageResult fname for PRAS ResultSpec
"""
function get_outage_result_fname(::Type{GeneratorStorageAvailability})
    return :gen_stor_availability
end

"""
Get SPIOutageResult fname for PRAS ResultSpec
"""
function get_outage_result_fname(::Type{R}) where {R <: PRASCore.Results.ResultSpec}
    return error("SPIOutageResult currently is not defined for $(R)")
end

"""
    SPIOutageResult(; shortfall_samples, gen_availability, stor_availability, gen_stor_availability)

# Arguments
$(TYPEDFIELDS)

SPIOutageResult is used to parse Tuple{Vararg{PRAS.PRASCore.Results.Result}} and add structure to it.
"""
mutable struct SPIOutageResult
    "Shortfall Sample Result"
    shortfall_samples::Union{Nothing, PRASCore.Results.ShortfallSamplesResult}
    "Generator Availability Result"
    gen_availability::Union{Nothing, PRASCore.Results.GeneratorAvailabilityResult}
    "Storage Availability Result"
    stor_availability::Union{Nothing, PRASCore.Results.StorageAvailabilityResult}
    "GeneratorStorage Availability Result"
    gen_stor_availability::Union{
        Nothing,
        PRASCore.Results.GeneratorStorageAvailabilityResult,
    }

    function SPIOutageResult(;
        shortfall_samples=nothing,
        gen_availability=nothing,
        stor_availability=nothing,
        gen_stor_availability=nothing,
    )
        return new(
            shortfall_samples,
            gen_availability,
            stor_availability,
            gen_stor_availability,
        )
    end
end

function SPIOutageResult(results::T) where {T <: Tuple{Vararg{PRASCore.Results.Result}}}
    num_fields = fieldcount(SPIOutageResult)
    outage_result = SPIOutageResult()
    for i in 1:num_fields
        fname = fieldname(SPIOutageResult, i)
        ftype = fieldtype(SPIOutageResult, i)
        if ~(isempty(filter(x -> typeof(x) <: ftype, results)))
            setproperty!(
                outage_result,
                fname,
                first(filter(x -> typeof(x) <: ftype, results)),
            )
        end
    end
    return outage_result
end
