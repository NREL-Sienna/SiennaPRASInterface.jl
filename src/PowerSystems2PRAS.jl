"""
Add default data to a system from `OUTAGE_INFO_FILE` (ERCOT data).
"""
function add_default_data!(sys::PSY.System, outage_info_file=OUTAGE_INFO_FILE)
    @warn "No forced outage data available in the Sienna/Data PowerSystems System. Using generic outage data ..."
    df_outage = DataFrames.DataFrame(
        CSV.File(
            outage_info_file,
            types=Dict(:tech => String, :PrimeMovers => String, :ThermalFuels => String),
            missingstring="NA",
        ),
    )

    outage_values = outage_data[]
    for row in eachrow(df_outage)
        if ismissing(row.ThermalFuels)
            push!(
                outage_values,
                outage_data(
                    PSY.PrimeMovers(row.PrimeMovers),
                    row.ThermalFuels,
                    row.NameplateLimit_MW,
                    (row.FOR / 100),
                    row.MTTR,
                ),
            )
        else
            push!(
                outage_values,
                outage_data(
                    PSY.PrimeMovers(row.PrimeMovers),
                    PSY.ThermalFuels(row.ThermalFuels),
                    row.NameplateLimit_MW,
                    (row.FOR / 100),
                    row.MTTR,
                ),
            )
        end
    end

    # Add min capacity fields to outage_data objects
    add_min_capacity!(outage_values)
    # Adding generic data to components in the System
    for outage_val in outage_values
        λ, μ = rate_to_probability(outage_val.FOR, outage_val.MTTR)
        transition_data = PSY.GeometricDistributionForcedOutage(;
            mean_time_to_recovery=outage_val.MTTR,
            outage_transition_probability=λ,
        )

        comps = if ismissing(outage_val.fuel)
            PSY.get_components(
                x -> (
                    PSY.get_prime_mover_type(x) == outage_val.prime_mover &&
                    outage_val.min_capacity <=
                    PSY.get_max_active_power(x) <
                    outage_val.max_capacity
                ),
                PSY.Generator,
                sys,
            )
        else
            PSY.get_components(
                x -> (
                    PSY.get_prime_mover_type(x) == outage_val.prime_mover &&
                    PSY.get_fuel(x) == outage_val.fuel &&
                    outage_val.min_capacity <=
                    PSY.get_max_active_power(x) <
                    outage_val.max_capacity
                ),
                PSY.ThermalGen,
                sys,
            )
        end

        for comp in comps
            PSY.add_supplemental_attribute!(sys, comp, transition_data)
        end
    end
end

"""
    $(TYPEDSIGNATURES)

Extract region load as a matrix of Int64 values.

Uses to all StaticLoad objects.

Note that behavior is not controlled by a formulation.
"""
function get_region_loads(sys::PSY.System, s2p_meta::S2P_metadata, regions)
    region_load = Array{Int64, 2}(undef, length(regions), s2p_meta.N)

    # FIXME We should make this work for all ElectricLoads
    for (idx, region) in enumerate(regions)
        reg_load_comps =
            get_available_components_in_aggregation_topology(PSY.StaticLoad, sys, region)
        if (length(reg_load_comps) > 0)
            region_load[idx, :] =
                floor.(
                    Int,
                    sum(
                        get_ts_values.(
                            get_first_ts.(
                                PSY.get_time_series_multiple.(
                                    reg_load_comps,
                                    s2p_meta.filter_func,
                                    name="max_active_power",
                                )
                            )
                        ) .* PSY.get_max_active_power.(reg_load_comps),
                    ),
                ) # Any issues with using the first of time_series_multiple?
        else
            region_load[idx, :] = zeros(Int64, s2p_meta.N)
        end
    end
    return region_load
end

# Generator must not have HydroEenergyReservoir or have 0 max active power or be a hybrid system
"""
    $(TYPEDSIGNATURES)

Extraction of generators using formulation dictionary to create a list of generators
and appropriate indices for PRAS. Note that objects with 0 max active power are excluded.
"""
function get_generator_region_indices(
    sys::PSY.System,
    s2p_meta::S2P_metadata,
    regions,
    component_to_formulation::Dict{PSY.Device, PRASGenerator},
)
    gens = Array{PSY.Device}[]
    start_id = Array{Int64}(undef, length(regions))
    region_gen_idxs = Array{UnitRange{Int64}, 1}(undef, length(regions))

    reg_wind_gens = []
    reg_pv_gens = []
    for (idx, region) in enumerate(regions)
        reg_ren_comps = filter(
            x -> haskey(component_to_formulation, x),
            get_available_components_in_aggregation_topology(PSY.RenewableGen, sys, region),
        )
        wind_gs = filter(
            x -> (
                (PSY.get_prime_mover_type(x) == PSY.PrimeMovers.WT) &&
                get_lump_renewable_generation(component_to_formulation[x])
            ),
            reg_ren_comps,
        )
        pv_gs = filter(
            x -> (
                (PSY.get_prime_mover_type(x) == PSY.PrimeMovers.PVe) &&
                get_lump_renewable_generation(component_to_formulation[x])
            ),
            reg_ren_comps,
        )
        gs = filter(
            x -> (
                haskey(component_to_formulation, x) &&
                !get_lump_renewable_generation(component_to_formulation[x]) &&
                !(iszero(PSY.get_max_active_power(x))) &&
                PSY.IS.get_uuid(x) ∉ s2p_meta.hs_uuids
            ),
            get_available_components_in_aggregation_topology(PSY.Generator, sys, region),
        )
        push!(gens, gs)
        push!(reg_wind_gens, wind_gs)
        push!(reg_pv_gens, pv_gs)

        if (idx == 1)
            start_id[idx] = 1
        else
            if (length(reg_wind_gens[idx - 1]) > 0 && length(reg_pv_gens[idx - 1]) > 0)
                start_id[idx] = start_id[idx - 1] + length(gens[idx - 1]) + 2
            elseif (length(reg_wind_gens[idx - 1]) > 0 || length(reg_pv_gens[idx - 1]) > 0)
                start_id[idx] = start_id[idx - 1] + length(gens[idx - 1]) + 1
            else
                start_id[idx] = start_id[idx - 1] + length(gens[idx - 1])
            end
        end

        if (length(reg_wind_gens[idx]) > 0 && length(reg_pv_gens[idx]) > 0)
            region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx]) + 2)
        elseif (length(reg_wind_gens[idx]) > 0 || length(reg_pv_gens[idx]) > 0)
            region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx]) + 1)
        else
            region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx]))
        end
    end
    for (gen, region, reg_wind_gen, reg_pv_gen) in
        zip(gens, regions, reg_wind_gens, reg_pv_gens)
        if (length(reg_wind_gen) > 0)
            # Wind
            temp_lumped_wind_gen = PSY.RenewableDispatch(nothing)
            PSY.set_name!(temp_lumped_wind_gen, "Lumped_Wind_" * PSY.get_name(region))
            PSY.set_prime_mover_type!(temp_lumped_wind_gen, PSY.PrimeMovers.WT)
            ext = PSY.get_ext(temp_lumped_wind_gen)
            ext["region_gens"] = reg_wind_gen
            push!(gen, temp_lumped_wind_gen)
        end
        if (length(reg_pv_gen) > 0)
            # PV
            temp_lumped_pv_gen = PSY.RenewableDispatch(nothing)
            PSY.set_name!(temp_lumped_pv_gen, "Lumped_PV_" * PSY.get_name(region))
            PSY.set_prime_mover_type!(temp_lumped_pv_gen, PSY.PrimeMovers.PVe)
            ext = PSY.get_ext(temp_lumped_pv_gen)
            ext["region_gens"] = reg_pv_gen
            push!(gen, temp_lumped_pv_gen)
        end
    end
    gen = reduce(vcat, gens)
    return gen, region_gen_idxs
end

"""
    $(TYPEDSIGNATURES)

Extraction of storage devices using formulation dictionary to create a list of storage
devices and appropriate indices for PRAS.
"""
function get_storage_region_indices(
    sys::PSY.System,
    s2p_meta::S2P_metadata,
    regions,
    component_to_formulation::Dict{PSY.Device, StoragePRAS},
)
    stors = Array{PSY.Device}[]
    start_id = Array{Int64}(undef, length(regions))
    region_stor_idxs = Array{UnitRange{Int64}, 1}(undef, length(regions))

    for (idx, region) in enumerate(regions)
        reg_stor_comps =
            get_available_components_in_aggregation_topology(PSY.Storage, sys, region)
        push!(
            stors,
            filter(
                x ->
                    haskey(component_to_formulation, x) &&
                        PSY.IS.get_uuid(x) ∉ s2p_meta.hs_uuids,
                reg_stor_comps,
            ),
        )
        idx == 1 ? start_id[idx] = 1 :
        start_id[idx] = start_id[idx - 1] + length(stors[idx - 1])
        region_stor_idxs[idx] = range(start_id[idx], length=length(stors[idx]))
    end
    return reduce(vcat, stors), region_stor_idxs
end

"""
    $(TYPEDSIGNATURES)

Extract components with a generator storage formulation.
"""
function get_gen_storage_region_indices(
    sys::PSY.System,
    regions,
    component_to_formulation::Dict{PSY.Device, PRASGeneratorStorage},
)
    gen_stors = Array{PSY.Device}[]
    start_id = Array{Int64}(undef, length(regions))
    region_genstor_idxs = Array{UnitRange{Int64}, 1}(undef, length(regions))

    for (idx, region) in enumerate(regions)
        reg_gen_stor_comps =
            get_available_components_in_aggregation_topology(PSY.Generator, sys, region)
        gs = filter(x -> haskey(component_to_formulation, x), reg_gen_stor_comps)
        push!(gen_stors, gs)
        idx == 1 ? start_id[idx] = 1 :
        start_id[idx] = start_id[idx - 1] + length(gen_stors[idx - 1])
        region_genstor_idxs[idx] = range(start_id[idx], length=length(gen_stors[idx]))
    end
    return reduce(vcat, gen_stors), region_genstor_idxs
end

"""
Turn a time series into an Array of floored ints
"""
function get_pras_array_from_timseries(device::PSY.Device, filter_func, name, multiplier)
    return floor.(
        Int,
        get_ts_values(
            get_first_ts(PSY.get_time_series_multiple(device, filter_func, name=name)),
        ) * multiplier,
    )
end

"""
    $(TYPEDSIGNATURES)

Apply PRASGenerator to process all generators objects
into rows in PRAS matrices:
- Capacity, λ, μ

Negative max active power will translate into zeros for time series data.
"""
function process_generators(
    gen::Array{PSY.Device},
    s2p_meta::S2P_metadata,
    component_to_formulation::Dict{PSY.Device, PRASGenerator},
)
    if (length(gen) == 0)
        gen_names = String[]
    else
        gen_names = PSY.get_name.(gen)
    end

    gen_categories = get_generator_category.(gen)
    n_gen = length(gen_names)

    gen_cap_array = Matrix{Int64}(undef, n_gen, s2p_meta.N)
    λ_gen = Matrix{Float64}(undef, n_gen, s2p_meta.N)
    μ_gen = Matrix{Float64}(undef, n_gen, s2p_meta.N)

    #FIXME This should use a component map instead.
    for (idx, g) in enumerate(gen)
        if haskey(PSY.get_ext(g), "region_gens")
            reg_gens_DA = PSY.get_ext(g)["region_gens"]
            gen_cap_array[idx, :] =
                round.(
                    Int,
                    sum([
                        get_ts_values(
                            get_first_ts(
                                PSY.get_time_series_multiple(
                                    reg_gen,
                                    s2p_meta.filter_func,
                                    name=get_max_active_power(
                                        component_to_formulation[reg_gen],
                                    ),
                                ),
                            ),
                        ) * PSY.get_max_active_power(reg_gen) for reg_gen in reg_gens_DA
                    ]),
                )
        else
            if (
                PSY.has_time_series(g) && (
                    get_max_active_power(component_to_formulation[g]) in
                    PSY.get_name.(PSY.get_time_series_multiple(g, s2p_meta.filter_func))
                )
            )
                gen_cap_array[idx, :] = get_pras_array_from_timseries(
                    g,
                    s2p_meta.filter_func,
                    get_max_active_power(component_to_formulation[g]),
                    PSY.get_max_active_power(g),
                )
                if !(all(gen_cap_array[idx, :] .>= 0))
                    @warn "There are negative values in max active time series data for $(PSY.get_name(g)) of type $(gen_categories[idx]) is negative. Using zeros for time series data."
                    gen_cap_array[idx, :] = zeros(Int, s2p_meta.N)
                end
            else
                if (PSY.get_max_active_power(g) > 0)
                    gen_cap_array[idx, :] =
                        fill.(floor.(Int, PSY.get_max_active_power(g)), 1, s2p_meta.N)
                else
                    @warn "Max active power for $(PSY.get_name(g)) of type $(gen_categories[idx]) is negative. Using zeros for time series data."
                    gen_cap_array[idx, :] = zeros(Int, s2p_meta.N) # to handle components with negative active power (usually UNAVAIALABLE)
                end
            end
        end

        λ_gen[idx, :], μ_gen[idx, :] = get_outage_time_series_data(g, s2p_meta)
    end

    return PRASCore.Generators{
        s2p_meta.N,
        s2p_meta.pras_timestep,
        s2p_meta.pras_resolution,
        PRASCore.MW,
    }(
        gen_names,
        get_generator_category.(gen),
        gen_cap_array,
        λ_gen,
        μ_gen,
    )
end

function assign_to_stor_matrices!(
    ::EnergyReservoirLossless,
    s::PSY.Device,
    s2p_meta::S2P_metadata,
    charge_cap_array,
    discharge_cap_array,
    energy_cap_array,
    chrg_eff_array,
    dischrg_eff_array,
)
    fill!(charge_cap_array, floor(Int, PSY.get_input_active_power_limits(s).max))
    fill!(discharge_cap_array, floor(Int, PSY.get_output_active_power_limits(s).max))
    fill!(
        energy_cap_array,
        floor(Int, PSY.get_storage_level_limits(s).max * PSY.get_storage_capacity(s)),
    )
    fill!(chrg_eff_array, PSY.get_efficiency(s).in)
    fill!(dischrg_eff_array, PSY.get_efficiency(s).out)
end

"""
    $(TYPEDSIGNATURES)

Apply StoragePRAS to process all storage objects
"""
function process_storage(
    stor::Array{PSY.Device},
    s2p_meta::S2P_metadata,
    component_to_formulation::Dict{PSY.Device, StoragePRAS},
)
    if (length(stor) == 0)
        stor_names = String[]
    else
        stor_names = PSY.get_name.(stor)
    end

    n_stor = length(stor_names)

    stor_charge_cap_array = Matrix{Int64}(undef, n_stor, s2p_meta.N)
    stor_discharge_cap_array = Matrix{Int64}(undef, n_stor, s2p_meta.N)
    stor_energy_cap_array = Matrix{Int64}(undef, n_stor, s2p_meta.N)
    stor_chrg_eff_array = Matrix{Float64}(undef, n_stor, s2p_meta.N)
    stor_dischrg_eff_array = Matrix{Float64}(undef, n_stor, s2p_meta.N)
    λ_stor = Matrix{Float64}(undef, n_stor, s2p_meta.N)
    μ_stor = Matrix{Float64}(undef, n_stor, s2p_meta.N)

    for (idx, s) in enumerate(stor)
        assign_to_stor_matrices!(
            component_to_formulation[g],
            s,
            s2p_meta,
            view(stor_charge_cap_array, idx, :),
            view(stor_discharge_cap_array, idx, :),
            view(Array{Int64}(undef, s2p_meta.N), :),  # Empty inflow array since Storage has no inflow
            view(stor_energy_cap_array, idx, :),
            view(Array{Int64}(undef, s2p_meta.N), :),  # Empty gridinj array since Storage has no grid injection
        )

        λ_stor[idx, :], μ_stor[idx, :] = get_outage_time_series_data(s, s2p_meta)
    end

    stor_cryovr_eff = ones(n_stor, s2p_meta.N)   # Not currently available/ defined in PowerSystems

    return PRASCore.Storages{
        s2p_meta.N,
        s2p_meta.pras_timestep,
        s2p_meta.pras_resolution,
        PRASCore.MW,
        PRASCore.MWh,
    }(
        stor_names,
        get_generator_category.(stor),
        stor_charge_cap_array,
        stor_discharge_cap_array,
        stor_energy_cap_array,
        stor_chrg_eff_array,
        stor_dischrg_eff_array,
        stor_cryovr_eff,
        λ_stor,
        μ_stor,
    )
end

"""
    $(TYPEDSIGNATURES)

Apply HybridSystem Formulation to fill in a row of a PRAS Matrix.
Views should be passed in for all arrays.
"""
function assign_to_gen_stor_matrices!(
    formulation::HybridSystemPRAS,
    g_s::PSY.Device,
    s2p_meta::S2P_metadata,
    charge_cap_array,
    discharge_cap_array,
    inflow_array,
    energy_cap_array,
    gridinj_cap_array,
)
    fill!(
        charge_cap_array,
        floor(Int, PSY.get_input_active_power_limits(PSY.get_storage(g_s)).max),
    )
    fill!(
        discharge_cap_array,
        floor(Int, PSY.get_output_active_power_limits(PSY.get_storage(g_s)).max),
    )
    fill!(
        energy_cap_array,
        floor(Int, PSY.get_state_of_charge_limits(PSY.get_storage(g_s)).max),
    )
    fill!(gridinj_cap_array, floor(Int, PSY.get_output_active_power_limits(g_s).max))

    if (
        PSY.has_time_series(PSY.get_renewable_unit(g_s)) && (
            get_max_active_power(formulation) in
            PSY.get_name.(
                PSY.get_time_series_multiple(
                    PSY.get_renewable_unit(g_s),
                    s2p_meta.filter_func,
                )
            )
        )
    )
        inflow_array .= get_pras_array_from_timseries(
            PSY.get_renewable_unit(g_s),
            s2p_meta.filter_func,
            get_max_active_power(formulation),
            PSY.get_max_active_power(PSY.get_renewable_unit(g_s)),
        )
    else
        fill!(
            inflow_array,
            floor(Int, PSY.get_max_active_power(PSY.get_renewable_unit(g_s))),
        )
    end
end

"""
    $(TYPEDSIGNATURES)

Apply HydroEnergyReservoir Formulation to fill in a row of a PRAS Matrix.
Views should be passed in for all arrays.
"""
function assign_to_gen_stor_matrices!(
    formulation::HydroEnergyReservoirPRAS,
    g_s::PSY.Device,
    s2p_meta::S2P_metadata,
    charge_cap_array,
    discharge_cap_array,
    inflow_array,
    energy_cap_array,
    gridinj_cap_array,
)
    if (PSY.has_time_series(g_s))
        if (
            get_inflow(formulation) in
            PSY.get_name.(PSY.get_time_series_multiple(g_s, s2p_meta.filter_func))
        )
            charge_cap_array .= get_pras_array_from_timseries(
                g_s,
                s2p_meta.filter_func,
                get_inflow(formulation),
                PSY.get_inflow(g_s),
            )
            discharge_cap_array .= charge_cap_array
            inflow_array .= charge_cap_array
        else
            fill!(charge_cap_array, floor(Int, PSY.get_inflow(g_s)))
            fill!(discharge_cap_array, floor(Int, PSY.get_inflow(g_s)))
            fill!(inflow_array, floor(Int, PSY.get_inflow(g_s)))
        end
        if (
            get_storage_capacity(storage_capacity) in
            PSY.get_name.(PSY.get_time_series_multiple(g_s, s2p_meta.filter_func))
        )
            energy_cap_array .= get_pras_array_from_timseries(
                g_s,
                s2p_meta.filter_func,
                get_storage_capacity(formulation),
                PSY.get_storage_capacity(g_s),
            )
        else
            fill!(energy_cap_array, floor(Int, PSY.get_storage_capacity(g_s)))
        end
        if (
            get_max_active_power(formulation) in
            PSY.get_name.(PSY.get_time_series_multiple(g_s, s2p_meta.filter_func))
        )
            gridinj_cap_array .= get_pras_array_from_timseries(
                g_s,
                s2p_meta.filter_func,
                get_max_active_power(formulation),
                PSY.get_max_active_power(g_s),
            )
        else
            fill!(gridinj_cap_array, floor(Int, PSY.get_max_active_power(g_s)))
        end
    else
        fill!(charge_cap_array, floor(Int, PSY.get_inflow(g_s)))
        fill!(discharge_cap_array, floor(Int, PSY.get_inflow(g_s)))
        fill!(energy_cap_array, floor(Int, PSY.get_storage_capacity(g_s)))
        fill!(inflow_array, floor(Int, PSY.get_inflow(g_s)))
        fill!(gridinj_cap_array, floor(Int, PSY.get_max_active_power(g_s)))
    end
end

"""
    $(TYPEDSIGNATURES)

Apply PRASGeneratorStorage to create PRAS matrices for generator storage
"""
function process_genstorage(
    gen_stor::Array{PSY.Device},
    s2p_meta::S2P_metadata,
    component_to_formulation::Dict{PSY.Device, PRASGeneratorStorage},
)
    if (length(gen_stor) == 0)
        gen_stor_names = String[]
    else
        gen_stor_names = PSY.get_name.(gen_stor)
    end

    n_genstors = length(gen_stor_names)

    gen_stor_charge_cap_array = Matrix{Int64}(undef, n_genstors, s2p_meta.N)
    gen_stor_discharge_cap_array = Matrix{Int64}(undef, n_genstors, s2p_meta.N)
    gen_stor_enrgy_cap_array = Matrix{Int64}(undef, n_genstors, s2p_meta.N)
    gen_stor_inflow_array = Matrix{Int64}(undef, n_genstors, s2p_meta.N)
    gen_stor_gridinj_cap_array = Matrix{Int64}(undef, n_genstors, s2p_meta.N)

    λ_genstors = Matrix{Float64}(undef, n_genstors, s2p_meta.N)
    μ_genstors = Matrix{Float64}(undef, n_genstors, s2p_meta.N)

    for (idx, g_s) in enumerate(gen_stor)
        assign_to_gen_stor_matrices!(
            component_to_formulation[g_s],
            g_s,
            s2p_meta,
            view(gen_stor_charge_cap_array, idx, :),
            view(gen_stor_discharge_cap_array, idx, :),
            view(gen_stor_inflow_array, idx, :),
            view(gen_stor_enrgy_cap_array, idx, :),
            view(gen_stor_gridinj_cap_array, idx, :),
        )

        λ_genstors[idx, :], μ_genstors[idx, :] = get_outage_time_series_data(g_s, s2p_meta)
    end

    gen_stor_gridwdr_cap_array = zeros(Int64, n_genstors, s2p_meta.N) # Not currently available/ defined in PowerSystems
    gen_stor_charge_eff = ones(n_genstors, s2p_meta.N)                # Not currently available/ defined in PowerSystems
    gen_stor_discharge_eff = ones(n_genstors, s2p_meta.N)             # Not currently available/ defined in PowerSystems
    gen_stor_cryovr_eff = ones(n_genstors, s2p_meta.N)                # Not currently available/ defined in PowerSystems

    return PRASCore.GeneratorStorages{
        s2p_meta.N,
        s2p_meta.pras_timestep,
        s2p_meta.pras_resolution,
        PRASCore.MW,
        PRASCore.MWh,
    }(
        gen_stor_names,
        get_generator_category.(gen_stor),
        gen_stor_charge_cap_array,
        gen_stor_discharge_cap_array,
        gen_stor_enrgy_cap_array,
        gen_stor_charge_eff,
        gen_stor_discharge_eff,
        gen_stor_cryovr_eff,
        gen_stor_inflow_array,
        gen_stor_gridwdr_cap_array,
        gen_stor_gridinj_cap_array,
        λ_genstors,
        μ_genstors,
    )
end

"""
    $(TYPEDSIGNATURES)

Use a RATemplate to create a PRAS system from a Sienna system.

# Arguments

- `sys::PSY.System`: Sienna PowerSystems System
- `template::RATemplate`: RATemplate
- `export_location::Union{Nothing, String}`: Export location for PRAS SystemModel

# Returns

- `PRASCore.SystemModel`: PRAS SystemModel

# Examples

```julia
generate_pras_system(sys, template)
```

Note that the original system will only be set to NATURAL_UNITS.
"""
function generate_pras_system(
    sys::PSY.System,
    template::RATemplate,
    export_location::Union{Nothing, String}=nothing,
)::PRASCore.SystemModel

    # PRAS needs Sienna\Data PowerSystems.jl System to be in NATURAL_UNITS
    PSY.set_units_base_system!(sys, PSY.UnitSystem.NATURAL_UNITS)

    # Check if any GeometricDistributionForcedOutage objects exist in the System
    outages = PSY.get_supplemental_attributes(PSY.GeometricDistributionForcedOutage, sys)

    # If no GeometricDistributionForcedOutage objects exist, add them to relevant components in the System
    if (outages.length == 0)
        add_default_data!(sys)
    end
    #######################################################
    # PRAS timestamps
    # Need this to select timeseries values of interest
    # TODO: Is it okay to assume each System will have a
    # SingleTimeSeries?
    #######################################################
    #
    s2p_meta = S2P_metadata()

    ts_counts = PSY.get_time_series_counts(sys)
    if (!iszero(ts_counts.static_time_series_count))
        s2p_meta.has_static_timeseries = true
        filter_func = x -> (typeof(x) <: PSY.StaticTimeSeries)
        s2p_meta.filter_func = filter_func
        all_ts = PSY.get_time_series_multiple(sys, filter_func)
        start_datetime = PSY.IS.get_initial_timestamp(first(all_ts))
        s2p_meta.first_timestamp = start_datetime
        s2p_meta.first_timeseries = first(all_ts)
    else
        if (!iszero(ts_counts.forecast_count))
            s2p_meta.has_forecasts = true
            filter_func = x -> (typeof(x) <: PSY.Forecast)
            s2p_meta.filter_func = filter_func
            all_ts = PSY.get_time_series_multiple(sys, filter_func)
            start_datetime = PSY.IS.get_initial_timestamp(first(all_ts))
            s2p_meta.first_timestamp = start_datetime
            s2p_meta.first_timeseries = first(all_ts)
        end
    end

    # Ensure Sienna/Data System has either static time series or forecasts
    if !(s2p_meta.has_static_timeseries) && !(s2p_meta.has_forecasts)
        error(
            "Sienna/Data PowerSystems System has no time series data (static time series data and/or forecasts)",
        )
    end
    # add N to S2P_metadata object
    add_N!(s2p_meta)

    # TODO: Is it okay to just get the first elemnt of vector returned by PSY.get_time_series_resolutions?
    sys_res =
        round(Dates.Millisecond(first(PSY.get_time_series_resolutions(sys))), Dates.Hour)
    if iszero(sys_res.value)
        sys_res = round(
            Dates.Millisecond(first(PSY.get_time_series_resolutions(sys))),
            Dates.Minute,
        )
        s2p_meta.pras_resolution = Dates.Minute
    end
    s2p_meta.pras_timestep = sys_res.value
    start_datetime_tz = TimeZones.ZonedDateTime(s2p_meta.first_timestamp, TimeZones.tz"UTC")
    finish_datetime_tz =
        start_datetime_tz + s2p_meta.pras_resolution((s2p_meta.N - 1) * sys_res)
    my_timestamps =
        StepRange(start_datetime_tz, s2p_meta.pras_resolution(sys_res), finish_datetime_tz)

    @info "The first timestamp of PRAS System being built is : $(start_datetime_tz) and last timestamp is : $(finish_datetime_tz) "
    #######################################################
    # Ensure no double counting of HybridSystem subcomponents
    # TODO: Not sure if we need this anymore.
    #######################################################
    dup_uuids = Base.UUID[]
    h_s_comps = PSY.get_available_components(PSY.HybridSystem, sys)
    for h_s in h_s_comps
        push!(dup_uuids, PSY.IS.get_uuid.(PSY._get_components(h_s))...)
    end
    # Add HybridSystem sub component UUIDs to s2p_meta
    if !(isempty(dup_uuids))
        s2p_meta.hs_uuids = dup_uuids
    end
    #######################################################
    # PRAS Regions - Areas in PowerSystems.jl
    #######################################################
    @info "Processing $(template.aggregation) objects in Sienna/Data PowerSystems System... "
    regions = collect(PSY.get_components(template.aggregation, sys))
    if !(length(regions) == 0)
        @info "The Sienna/Data PowerSystems System has $(length(regions)) regions based on PSY AggregationTopology : $(template.aggregation)."
    else
        error(
            "No regions in the Sienna/Data PowerSystems System. Cannot proceed with the process of making a PRAS SystemModel.",
        )
    end

    region_load = get_region_loads(sys, s2p_meta, regions)
    new_regions =
        PRASCore.Regions{s2p_meta.N, PRASCore.MW}(PSY.get_name.(regions), region_load)

    @info "Processing Generators in PSY System... "
    gens_to_formula =
        build_component_to_formulation(PRASGenerator, sys, template.device_models)
    gens, region_gen_idxs =
        get_generator_region_indices(sys, s2p_meta, regions, gens_to_formula)
    new_generators = process_generators(gens, s2p_meta, gens_to_formula)

    # **TODO Future : time series for storage devices
    @info "Processing Storages in PSY System... "
    stors_to_formula =
        build_component_to_formulation(StoragePRAS, sys, template.device_models)
    stors, region_stor_idxs =
        get_storage_region_indices(sys, s2p_meta, regions, stors_to_formula)
    new_storage = process_storage(stors, s2p_meta, stors_to_formula)

    # **TODO Consider all combinations of HybridSystem (Currently only works for DER+ESS)
    @info "Processing GeneratorStorages in PSY System... "
    gen_stors_to_formula =
        build_component_to_formulation(PRASGeneratorStorage, sys, template.device_models)
    gen_stors, region_genstor_idxs =
        get_gen_storage_region_indices(sys, regions, gen_stors_to_formula)
    new_gen_stors = process_genstorage(gen_stors, s2p_meta, gen_stors_to_formula)

    #######################################################
    # Network
    #######################################################
    if (length(regions) > 1)
        #######################################################
        # PRAS Lines
        #######################################################
        @info "Collecting all inter regional lines in Sienna/Data PowerSystems System..."

        lines = collect(
            PSY.get_components(
                x -> (
                    PSY.get_available(x) &&
                    typeof(x) ∉ TransformerTypes &&
                    typeof(x) != PSY.AreaInterchange
                ),
                PSY.Branch,
                sys,
            ),
        )
        #######################################################
        # Inter-Regional Line Processing
        #######################################################
        inter_regional_lines =
            filter(x -> !(x.arc.from.area.name == x.arc.to.area.name), lines)
        sorted_lines, interface_reg_idxs, interface_line_idxs =
            get_sorted_lines(inter_regional_lines, PSY.get_name.(regions))
        new_lines, new_interfaces = make_pras_interfaces(
            sorted_lines,
            interface_reg_idxs,
            interface_line_idxs,
            PSY.get_name.(regions),
            sys,
            s2p_meta,
        )

        pras_system = PRASCore.SystemModel(
            new_regions,
            new_interfaces,
            new_generators,
            region_gen_idxs,
            new_storage,
            region_stor_idxs,
            new_gen_stors,
            region_genstor_idxs,
            new_lines,
            interface_line_idxs,
            my_timestamps,
        )

        @info "Successfully built a PRAS SystemModel of type $(typeof(pras_system))."

        export_pras_system(pras_system, export_location::Union{Nothing, String})

        return pras_system

    else
        load_vector = vec(sum(region_load, dims=1))
        pras_system = PRASCore.SystemModel(
            new_generators,
            new_storage,
            new_gen_stors,
            my_timestamps,
            load_vector,
        )

        @info "Successfully built a PRAS SystemModel of type $(typeof(pras_system))."

        export_pras_system(pras_system, export_location::Union{Nothing, String})

        return pras_system
    end
end

const DEFAULT_DEVICE_MODELS = [
    DeviceRAModel(PSY.ThermalGen, PRASGenerator),
    DeviceRAModel(PSY.RenewableGen, PRASGenerator),
    DeviceRAModel(PSY.HydroDispatch, PRASGenerator),
    DeviceRAModel(PSY.EnergyReservoirStorage, EnergyReservoirLossless),
    DeviceRAModel(PSY.HybridSystem, HybridSystemPRAS),
    DeviceRAModel(PSY.HydroEnergyReservoir, HydroEnergyReservoirPRAS),
]

const _LUMPED_RENEWABLE_DEVICE_MODELS = [
    DeviceRAModel(PSY.ThermalGen, PRASGenerator),
    DeviceRAModel(PSY.RenewableGen, PRASGenerator, lump_renewable_generation=true),
    DeviceRAModel(PSY.HydroDispatch, PRASGenerator),
    DeviceRAModel(PSY.EnergyReservoirStorage, EnergyReservoirLossless),
    DeviceRAModel(PSY.HybridSystem, HybridSystemPRAS),
    DeviceRAModel(PSY.HydroEnergyReservoir, HydroEnergyReservoirPRAS),
]

const DEFAULT_TEMPLATE = RATemplate(PSY.Area, DEFAULT_DEVICE_MODELS)

"""
    $(TYPEDSIGNATURES)

Sienna/Data PowerSystems.jl System is the input and an object of PRAS SystemModel is returned.
...

# Arguments

  - `sys::PSY.System`: Sienna/Data PowerSystems.jl System
  - `aggregation<:PSY.AggregationTopology`: "PSY.Area" (or) "PSY.LoadZone" {Optional}
  - `lump_region_renewable_gens::Bool`: Whether to lumps PV and Wind generators in a region because usually these generators don't have FOR data {Optional}
  - `export_location::String`: Export location of the .pras file
    ...

# Returns

    - `PRASCore.SystemModel`: PRAS SystemModel object

# Examples

```julia-repl
julia> generate_pras_system(psy_sys, PSY.Area)
PRAS SystemModel
```
"""
function generate_pras_system(
    sys::PSY.System,
    aggregation::Type{AT},
    lump_region_renewable_gens::Bool=false,
    export_location::Union{Nothing, String}=nothing,
)::PRASCore.SystemModel where {AT <: PSY.AggregationTopology}
    if lump_region_renewable_gens
        template = RATemplate(aggregation, _LUMPED_RENEWABLE_DEVICE_MODELS)
    else
        template = RATemplate(aggregation, DEFAULT_DEVICE_MODELS)
    end
    generate_pras_system(sys, template, export_location)
end

"""
    generate_pras_system(sys_location::String, aggregation; kwargs...)

Generate a PRAS SystemModel from a Sienna/Data PowerSystems System JSON file.

# Arguments

  - `sys_location::String`: Location of the Sienna/Data PowerSystems System JSON file
  - `aggregation::Type{AT}`: Aggregation topology type
  - `lump_region_renewable_gens::Bool`: Lumping of region renewable generators
  - `export_location::Union{Nothing, String}`: Export location of the .pras file

# Returns

  - `PRASCore.SystemModel`: PRAS SystemModel
"""
function generate_pras_system(
    sys_location::String,
    aggregation::Type{AT},
    lump_region_renewable_gens=false,
    export_location::Union{Nothing, String}=nothing,
) where {AT <: PSY.AggregationTopology}
    @info "Running checks on the Sienna/Data PowerSystems System location provided ..."
    runchecks(sys_location)

    @info "The Sienna/Data PowerSystems System is being de-serialized from the System JSON ..."
    sys = try
        PSY.System(sys_location; time_series_read_only=true, runchecks=false)
    catch
        error(
            "Sienna/Data PowerSystems System could not be de-serialized using the location of JSON provided. Please check the location and make sure you have permission to access time_series_storage.h5",
        )
    end

    generate_pras_system(sys, aggregation, lump_region_renewable_gens, export_location)
end
