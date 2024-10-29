#######################################################
# Main function to build a instance of PRAS SystemModel
# from Sienna/Data PowerSytems.jl System
#######################################################
function make_pras_system(sys::PSY.System;
                          aggregation::Union{Nothing, Type{AT}} = nothing, availability=true, lump_region_renewable_gens=false,
                          export_location::Union{Nothing, String} = nothing) where {AT<:PSY.AggregationTopology}
    """
    make_pras_system(psy_sys,system_model)

    Sienna/Data PowerSystems.jl System is the input and an object of PRAS SystemModel is returned.
    ...
    # Arguments
    - `sys::PSY.System`: Sienna/Data PowerSystems.jl System
    - `aggregation<:PSY.AggregationTopology`: "PSY.Area" (or) "PSY.LoadZone" {Optional} 
    - `availability::Bool`: Takes into account avaialability of StaticInjection components when building the PRAS System {Optional} 
    - `lump_region_renewable_gens::Bool`: Whether to lumps PV and Wind generators in a region because usually these generators don't have FOR data {Optional} 
    - `export_location::String`: Export location of the .pras file  
    ...

    # Examples
    ```julia-repl
    julia> make_pras_system(psy_sys)
    PRAS SystemModel
    ```
    """
    # PRAS needs Sienna\Data PowerSystems.jl System to be in NATURAL_UNITS
    PSY.set_units_base_system!(sys, PSY.UnitSystem.NATURAL_UNITS); 

    # Check if any GeometricDistributionForcedOutage objects exist in the System
    outages = PSY.get_supplemental_attributes(PSY.GeometricDistributionForcedOutage, sys)
    
    # If no GeometricDistributionForcedOutage objects exist, add them to relevant components in the System
    if (outages.length == 0)
        @warn "No forced outage data available in the Sienna/Data System. Using generic outage data ..."
        df_outage = DataFrames.DataFrame(CSV.File(OUTAGE_INFO_FILE, types=Dict(:tech => String, :PrimeMovers => String, :ThermalFuels => String), missingstring = "NA"));
        
        outage_values =outage_data[]
        for row in eachrow(df_outage)
            if (ismissing(row.ThermalFuels))
                push!(outage_values, outage_data(PSY.PrimeMovers(row.PrimeMovers),row.ThermalFuels,row.NameplateLimit_MW,(row.FOR/100),row.MTTR))
            else
                push!(outage_values, outage_data(PSY.PrimeMovers(row.PrimeMovers),PSY.ThermalFuels(row.ThermalFuels),row.NameplateLimit_MW,(row.FOR/100),row.MTTR))
            end
        end
        
        # Add min capacity fields to outage_data objects
        add_min_capacity!(outage_values)
        # Adding generic data to components in the System
        for outage_val in outage_values
            λ, μ = rate_to_probability(outage_val.FOR, outage_val.MTTR)
            transition_data = PSY.GeometricDistributionForcedOutage(;mean_time_to_recovery=outage_val.MTTR,outage_transition_probability=λ)

            comps = 
            if ismissing(outage_val.fuel)
                PSY.get_components(x -> (PSY.get_prime_mover_type(x) ==  outage_val.prime_mover && outage_val.min_capacity <= PSY.get_max_active_power(x) < outage_val.max_capacity), PSY.Generator, sys)
            else
                PSY.get_components(x -> (PSY.get_prime_mover_type(x) ==  outage_val.prime_mover && PSY.get_fuel(x) ==  outage_val.fuel && outage_val.min_capacity <= PSY.get_max_active_power(x) < outage_val.max_capacity) , PSY.ThermalGen, sys)
            end

            for comp in comps
                PSY.add_supplemental_attribute!(sys, comp, transition_data)
            end
        end
    end
    #######################################################
    # PRAS timestamps
    # Need this to select timeseries values of interest
    # TODO: Is it okay to assume each System will have a 
    # SingleTimeSeries? 
    #######################################################
    # 
    s2p_meta = S2P_metadata()
    
    ts_counts = PSY.get_time_series_counts(sys);
    if (~iszero(ts_counts.static_time_series_count))
        s2p_meta.has_static_timeseries = true
        filter_func = x -> (typeof(x) <: PSY.StaticTimeSeries)
        s2p_meta.filter_func = filter_func
        all_ts = PSY.get_time_series_multiple(sys, filter_func);
        start_datetime = PSY.IS.get_initial_timestamp(first(all_ts));
        s2p_meta.first_timestamp = start_datetime
        s2p_meta.first_timeseries = first(all_ts)
    else
        if (~iszero(ts_counts.forecast_count))
            s2p_meta.has_forecasts = true
            filter_func = x -> (typeof(x) <: PSY.Forecast)
            s2p_meta.filter_func = filter_func
            all_ts = PSY.get_time_series_multiple(sys, filter_func);
            start_datetime = PSY.IS.get_initial_timestamp(first(all_ts));
            s2p_meta.first_timestamp = start_datetime
            s2p_meta.first_timeseries = first(all_ts)
        end
    end

    # Ensure Sienna/Data System has either static time series or forecasts
    if ~(s2p_meta.has_static_timeseries) && ~(s2p_meta.has_forecasts)
        error("Sienna/Data PowerSystems System has no time series data (static time series data and/or forecasts)")
    end
    # add N to S2P_metadata object
    add_N!(s2p_meta)

   # TODO: Is it okay to just get the first elemnt of vector returned by PSY.get_time_series_resolutions?
    sys_res_in_hour = round(Dates.Millisecond(first(PSY.get_time_series_resolutions(sys))), Dates.Hour);
    start_datetime_tz = TimeZones.ZonedDateTime(s2p_meta.first_timestamp,TimeZones.tz"UTC");
    finish_datetime_tz = start_datetime_tz +  Dates.Hour((s2p_meta.N-1)*sys_res_in_hour);
    my_timestamps = StepRange(start_datetime_tz, Dates.Hour(sys_res_in_hour), finish_datetime_tz);

    @info "The first timestamp of PRAS System being built is : $(start_datetime_tz) and last timestamp is : $(finish_datetime_tz) "
    #######################################################
    # Ensure no double counting of HybridSystem subcomponents
    # TODO: Not sure if we need this anymore.
    #######################################################
    dup_uuids =Base.UUID[];
    h_s_comps = availability ? PSY.get_components(PSY.get_available, PSY.HybridSystem, sys) : PSY.get_components(PSY.HybridSystem, sys)
    for h_s in h_s_comps
        push!(dup_uuids,PSY.IS.get_uuid.(PSY._get_components(h_s))...)
    end
    # Add HybridSystem sub component UUIDs to s2p_meta 
    if ~(isempty(dup_uuids))
        s2p_meta.hs_uuids = dup_uuids
    end

    # TODO: Do we still need to do this? From now, PSS/e parser
    # will return PSY.StandardLoad objects
    if (length(PSY.get_components(PSY.PowerLoad,sys)) > 0)
        s2p_meta.load_type = PSY.PowerLoad
    end
    #######################################################
     # PRAS Regions - Areas in SIIP
    #######################################################
    @info "Processing $(aggregation) objects in Sienna/Data PowerSystems System... "
    regions = collect(PSY.get_components(aggregation, sys));
    if ~(length(regions)==0)
        @info "The Sienna/Data PowerSystems System has $(length(regions)) regions based on PSY AggregationTopology : $(aggregation)."
    else
        error("No regions in the Sienna/Data PowerSystems System. Cannot proceed with the process of making a PRAS SystemModel.")
    end 

    region_names = PSY.get_name.(regions);
    num_regions = length(region_names);

    region_load = Array{Int64,2}(undef,num_regions,s2p_meta.N);
   
    for (idx,region) in enumerate(regions)
        reg_load_comps = availability ? get_available_components_in_aggregation_topology(s2p_meta.load_type, sys, region) :
                                             PSY.get_components_in_aggregation_topology(s2p_meta.load_type, sys, region)
        if (length(reg_load_comps) > 0)

            region_load[idx,:]=floor.(Int,sum(get_ts_values.(get_first_ts.(PSY.get_time_series_multiple.(reg_load_comps,s2p_meta.filter_func, name = "max_active_power")))
                               .*PSY.get_max_active_power.(reg_load_comps))); # Any issues with using the first of time_series_multiple?
        else
            region_load[idx,:] = zeros(Int64,s2p_meta.N)
        end
    end

    new_regions = PRAS.Regions{s2p_meta.N,PRAS.MW}(region_names, region_load);
    #######################################################
    # Generator Region Indices
    #######################################################
    gens=Array{PSY.Generator}[];
    start_id = Array{Int64}(undef,num_regions); 
    region_gen_idxs = Array{UnitRange{Int64},1}(undef,num_regions); 
    reg_wind_gens = []
    reg_pv_gens = []

    if (lump_region_renewable_gens)
        for (idx,region) in enumerate(regions)
            reg_ren_comps = availability ? get_available_components_in_aggregation_topology(PSY.RenewableGen, sys, region) :
                                                 PSY.get_components_in_aggregation_topology(PSY.RenewableGen, sys, region)
            wind_gs = filter(x -> (PSY.get_prime_mover_type(x) == PSY.PrimeMovers.WT), reg_ren_comps)
            pv_gs =  filter(x -> (PSY.get_prime_mover_type(x) == PSY.PrimeMovers.PVe), reg_ren_comps)
            reg_gen_comps = availability ? get_available_components_in_aggregation_topology(PSY.Generator, sys, region) :
                                                PSY.get_components_in_aggregation_topology(PSY.Generator, sys, region)
            gs = filter(x -> (~(typeof(x) == PSY.HydroEnergyReservoir) && ~(iszero(PSY.get_max_active_power(x))) &&  PSY.IS.get_uuid(x) ∉ union(s2p_meta.hs_uuids,PSY.IS.get_uuid.(wind_gs),PSY.IS.get_uuid.(pv_gs))),reg_gen_comps)
            push!(gens,gs)
            push!(reg_wind_gens,wind_gs)
            push!(reg_pv_gens,pv_gs)

            if (idx==1)
                start_id[idx] = 1
            else 
                if (length(reg_wind_gens[idx-1]) > 0 && length(reg_pv_gens[idx-1]) > 0)
                    start_id[idx] =start_id[idx-1]+length(gens[idx-1])+2
                elseif (length(reg_wind_gens[idx-1]) > 0 || length(reg_pv_gens[idx-1]) > 0)
                    start_id[idx] =start_id[idx-1]+length(gens[idx-1])+1
                else
                    start_id[idx] =start_id[idx-1]+length(gens[idx-1])
                end
            end

            if (length(reg_wind_gens[idx]) > 0 && length(reg_pv_gens[idx]) > 0)
                region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx])+2)
            elseif (length(reg_wind_gens[idx]) > 0 || length(reg_pv_gens[idx]) > 0)
                region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx])+1)
            else
                region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx]))
            end
        end
    else
        for (idx,region) in enumerate(regions)
            reg_gen_comps = availability_flag ? get_available_components_in_aggregation_topology(PSY.Generator, sys, region) :
                                                PSY.get_components_in_aggregation_topology(PSY.Generator, sys, region)
            gs= [g for g in reg_gen_comps if (typeof(g) != PSY.HydroEnergyReservoir && PSY.get_max_active_power(g)!=0 && PSY.IS.get_uuid(g) ∉ dup_uuids)]
            gs = filter(x -> (~(typeof(x) == PSY.HydroEnergyReservoir) && ~(iszero(PSY.get_max_active_power(x))) &&  PSY.IS.get_uuid(x) ∉ s2p_meta.hs_uuids),reg_gen_comps)
            push!(gens,gs)
            idx==1 ? start_id[idx] = 1 : start_id[idx] =start_id[idx-1]+length(gens[idx-1])
            region_gen_idxs[idx] = range(start_id[idx], length=length(gens[idx]))
        end
    end
    #######################################################
    # Storages Region Indices
    #######################################################
    stors=[];
    start_id = Array{Int64}(undef,num_regions);
    region_stor_idxs = Array{UnitRange{Int64},1}(undef,num_regions);

    for (idx,region) in enumerate(regions)
        #push!(stors,[s for s in PSY.get_components_in_aggregation_topology(PSY.Storage, sys, region)])
        reg_stor_comps = availability ? get_available_components_in_aggregation_topology(PSY.Storage, sys, region) :
                                             PSY.get_components_in_aggregation_topology(PSY.Storage, sys, region)
        push!(stors,filter(x -> (PSY.IS.get_uuid(x) ∉ dup_uuids), reg_stor_comps))
        idx==1 ? start_id[idx] = 1 : start_id[idx] =start_id[idx-1]+length(stors[idx-1])
        region_stor_idxs[idx] = range(start_id[idx], length=length(stors[idx]))
    end
    #######################################################
    # GeneratorStorages Region Indices
    #######################################################
    gen_stors=[];
    start_id = Array{Int64}(undef,num_regions);
    region_genstor_idxs = Array{UnitRange{Int64},1}(undef,num_regions);

    for (idx,region) in enumerate(regions)
        reg_gen_stor_comps = availability ? get_available_components_in_aggregation_topology(PSY.Generator, sys, region) :
                                                 PSY.get_components_in_aggregation_topology(PSY.Generator, sys, region)
        gs= filter(x -> (typeof(x) == PSY.HydroEnergyReservoir || typeof(x)==PSY.HybridSystem), reg_gen_stor_comps)
        push!(gen_stors,gs)
        idx==1 ? start_id[idx] = 1 : start_id[idx] =start_id[idx-1]+length(gen_stors[idx-1])
        region_genstor_idxs[idx] = range(start_id[idx], length=length(gen_stors[idx]))
    end
    #######################################################
    # PRAS Generators
    #######################################################
    @info "Processing Generators in PSY System... "
    
    # Lumping Wind and PV Generators per Region
    if (lump_region_renewable_gens)
        for i in 1: num_regions
            if (length(reg_wind_gens[i])>0)
                # Wind
                temp_lumped_wind_gen = PSY.RenewableDispatch(nothing)
                PSY.set_name!(temp_lumped_wind_gen,"Lumped_Wind_"*region_names[i])
                PSY.set_prime_mover_type!(temp_lumped_wind_gen,PSY.PrimeMovers.WT)
                ext = PSY.get_ext(temp_lumped_wind_gen)
                ext["region_gens"] = reg_wind_gens[i]
                push!(gens[i],temp_lumped_wind_gen)
            end
            if (length(reg_pv_gens[i])>0)
                # PV
                temp_lumped_pv_gen = PSY.RenewableDispatch(nothing)
                PSY.set_name!(temp_lumped_pv_gen,"Lumped_PV_"*region_names[i])
                PSY.set_prime_mover_type!(temp_lumped_pv_gen,PSY.PrimeMovers.PVe)
                ext = PSY.get_ext(temp_lumped_pv_gen)
                ext["region_gens"] = reg_pv_gens[i]
                push!(gens[i],temp_lumped_pv_gen)
            end
        end
    end

    gen=[];
    for i in 1: num_regions
        if (length(gens[i]) != 0)
            push!(gen,gens[i]...)
        end
    end
    
    if(length(gen) ==0)
        gen_names = String[];
    else
        gen_names = PSY.get_name.(gen);
    end

    gen_categories = get_generator_category.(gen);
    n_gen = length(gen_names);

    gen_cap_array = Matrix{Int64}(undef, n_gen, s2p_meta.N);
    λ_gen = Matrix{Float64}(undef, n_gen, s2p_meta.N);
    μ_gen = Matrix{Float64}(undef, n_gen, s2p_meta.N);

    for (idx,g) in enumerate(gen) 
        if (lump_region_renewable_gens && (PSY.get_prime_mover_type(g) == PSY.PrimeMovers.WT || PSY.get_prime_mover_type(g) == PSY.PrimeMovers.PVe))
            reg_gens_DA = PSY.get_ext(g)["region_gens"];
            gen_cap_array[idx,:] = round.(Int,sum(get_ts_values.(get_first_ts.(PSY.get_time_series_multiple.(reg_gens_DA,s2p_meta.filter_func, name = "max_active_power")))
                                   .*PSY.get_max_active_power.(reg_gens_DA)));
        else
            if (PSY.has_time_series(g) && ("max_active_power" in PSY.get_name.(PSY.get_time_series_multiple(g,s2p_meta.filter_func))))
                gen_cap_array[idx,:] = floor.(Int,get_ts_values(get_first_ts(PSY.get_time_series_multiple(g,s2p_meta.filter_func, name = "max_active_power")))
                                       *PSY.get_max_active_power(g));
                if ~(all(gen_cap_array[idx,:] .>=0))
                    @warn "There are negative values in max active time series data for $(PSY.get_name(g)) of type $(gen_categories[idx]) is negative. Using zeros for time series data." 
                    gen_cap_array[idx,:] = zeros(Int,s2p_meta.N);
                end
            else
                if (PSY.get_max_active_power(g) > 0)
                    gen_cap_array[idx,:] = fill.(floor.(Int,PSY.get_max_active_power(g)),1,s2p_meta.N);
                else
                    @warn "Max active power for $(PSY.get_name(g)) of type $(gen_categories[idx]) is negative. Using zeros for time series data." 
                    gen_cap_array[idx,:] = zeros(Int,s2p_meta.N); # to handle components with negative active power (usually UNAVAIALABLE)
                end
            end
        end

        # Get GeometricForcedOutage SupplementalAttribute of the generator g
        outage_sup_attrs = PSY.get_supplemental_attributes(PSY.GeometricDistributionForcedOutage, g)
        
        if (length(outage_sup_attrs) > 0)
            transition_data = first(outage_sup_attrs)
            λ = PSY.get_outage_transition_probability(transition_data)
            μ = 1 / PSY.get_mean_time_to_recovery(transition_data) 

            if (PSY.has_time_series(transition_data, PSY.SingleTimeSeries))
                λ_gen[idx,:] = PSY.get_time_series_values(PSY.SingleTimeSeries, transition_data, "outage_probability")
                μ_gen[idx,:] = PSY.get_time_series_values(PSY.SingleTimeSeries, transition_data, "recovery_probability")
            else
                λ_gen[idx,:] = fill.(λ,1,s2p_meta.N); 
                μ_gen[idx,:] = fill.(μ,1,s2p_meta.N); 
            end
        else
            λ_gen[idx,:] = zeros(Float64,1,s2p_meta.N); 
            μ_gen[idx,:] = ones(Float64,1,s2p_meta.N); 
        end
       
    end

    new_generators = PRAS.Generators{s2p_meta.N,1,PRAS.Hour,PRAS.MW}(gen_names, get_generator_category.(gen), gen_cap_array , λ_gen ,μ_gen);
        
    #######################################################
    # PRAS Storages
    # **TODO Future : time series for storage devices
    #######################################################
    @info "Processing Storages in PSY System... "

    stor=[];
    for i in 1: num_regions
        if (length(stors[i]) != 0)
            append!(stor,stors[i])
        end
    end

    if(length(stor) ==0)
        stor_names=String[];
    else
        stor_names = PSY.get_name.(stor);
    end

    stor_categories = get_generator_category.(stor);

    n_stor = length(stor_names);

    stor_charge_cap_array = Matrix{Int64}(undef, n_stor, N);
    stor_discharge_cap_array = Matrix{Int64}(undef, n_stor, N);
    stor_energy_cap_array = Matrix{Int64}(undef, n_stor, N);
    stor_chrg_eff_array = Matrix{Float64}(undef, n_stor, N);
    stor_dischrg_eff_array  = Matrix{Float64}(undef, n_stor, N);
    λ_stor = Matrix{Float64}(undef, n_stor, N);   
    μ_stor = Matrix{Float64}(undef, n_stor, N);

    for (idx,s) in enumerate(stor)
        stor_charge_cap_array[idx,:] = fill(floor(Int,getfield(PSY.get_input_active_power_limits(s), :max)),1,N);
        stor_discharge_cap_array[idx,:] = fill(floor(Int,getfield(PSY.get_output_active_power_limits(s), :max)),1,N);
        stor_energy_cap_array[idx,:] = fill(floor(Int,getfield(PSY.get_state_of_charge_limits(s),:max)),1,N);
        stor_chrg_eff_array[idx,:] = fill(getfield(PSY.get_efficiency(s), :in),1,N);
        stor_dischrg_eff_array[idx,:]  = fill.(getfield(PSY.get_efficiency(s), :out),1,N);

        if (~outage_flag)
            @warn "No outage information is available for $(PSY.get_name(s)) of type $(stor_categories[idx]). Using nominal outage and recovery probabilities for this generator."
            λ = 0.0;
            μ = 1.0;
        else
            ext = PSY.get_ext(s)
            if (!(haskey(ext,"outage_probability") && haskey(ext,"recovery_probability")))
                @warn "No outage information is available in ext field of $(PSY.get_name(s)) of type $(stor_categories[idx]). Using nominal outage and recovery probabilities for this generator."
                λ = 0.0;
                μ = 1.0;
            else
                λ = ext["outage_probability"];
                μ = ext["recovery_probability"];
            end
        end
        
        λ_stor[idx,:] = fill.(λ,1,N); 
        μ_stor[idx,:] = fill.(μ,1,N); 
    end
    
    stor_cryovr_eff = ones(n_stor,N);   # Not currently available/ defined in PowerSystems
    
    new_storage = PRAS.Storages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(stor_names,get_generator_category.(stor),
                                            stor_charge_cap_array,stor_discharge_cap_array,stor_energy_cap_array,
                                            stor_chrg_eff_array,stor_dischrg_eff_array, stor_cryovr_eff,
                                            λ_stor,μ_stor);

    #######################################################
    # PRAS Generator Storages
    # **TODO Consider all combinations of HybridSystem (Currently only works for DER+ESS)
    #######################################################
    @info "Processing GeneratorStorages in PSY System... "

    gen_stor=[];
    for i in 1: num_regions
        if (length(gen_stors[i]) != 0)
            append!(gen_stor,gen_stors[i])
        end
    end
    
    if(length(gen_stor) == 0)
        gen_stor_names=String[];
    else
        gen_stor_names = PSY.get_name.(gen_stor);
    end

    gen_stor_categories = string.(typeof.(gen_stor)); 
    
    n_genstors = length(gen_stor_names);

    gen_stor_charge_cap_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_discharge_cap_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_enrgy_cap_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_inflow_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_gridinj_cap_array = Matrix{Int64}(undef, n_genstors, N);

    λ_genstors = Matrix{Float64}(undef, n_genstors, N);   
    μ_genstors = Matrix{Float64}(undef, n_genstors, N);  

    for (idx,g_s) in enumerate(gen_stor)
        if(typeof(g_s) ==PSY.HydroEnergyReservoir)
            if (PSY.has_time_series(g_s))
                if ("inflow" in PSY.get_name.(PSY.get_time_series_multiple(g_s)))
                    gen_stor_charge_cap_array[idx,:] = floor.(Int,get_forecast_values(get_first_ts(PSY.get_time_series_multiple(g_s, name = "inflow")))
                                                       *PSY.get_inflow(g_s));
                    gen_stor_discharge_cap_array[idx,:] = floor.(Int,get_forecast_values(get_first_ts(PSY.get_time_series_multiple(g_s, name = "inflow")))
                                                          *PSY.get_inflow(g_s));
                    gen_stor_inflow_array[idx,:] = floor.(Int,get_forecast_values(get_first_ts(PSY.get_time_series_multiple(g_s, name = "inflow")))
                                                   *PSY.get_inflow(g_s));
                else
                    gen_stor_charge_cap_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                    gen_stor_discharge_cap_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                    gen_stor_inflow_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                end
                if ("storage_capacity" in PSY.get_name.(PSY.get_time_series_multiple(g_s)))
                    gen_stor_enrgy_cap_array[idx,:] = floor.(Int,get_forecast_values(get_first_ts(PSY.get_time_series_multiple(g_s, name = "storage_capacity")))
                                                      *PSY.get_storage_capacity(g_s));
                else
                    gen_stor_enrgy_cap_array[idx,:] = fill.(floor.(Int,PSY.get_storage_capacity(g_s)),1,N);
                end
                if ("max_active_power" in PSY.get_name.(PSY.get_time_series_multiple(g_s)))
                    gen_stor_gridinj_cap_array[idx,:] = floor.(Int,get_forecast_values(get_first_ts(PSY.get_time_series_multiple(g_s, name = "max_active_power")))
                                                        *PSY.get_max_active_power(g_s));
                else
                    gen_stor_gridinj_cap_array[idx,:] = fill.(floor.(Int,PSY.get_max_active_power(g_s)),1,N);
                end
            else
                gen_stor_charge_cap_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                gen_stor_discharge_cap_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                gen_stor_enrgy_cap_array[idx,:] = fill.(floor.(Int,PSY.get_storage_capacity(g_s)),1,N);
                gen_stor_inflow_array[idx,:] = fill.(floor.(Int,PSY.get_inflow(g_s)),1,N);
                gen_stor_gridinj_cap_array[idx,:] = fill.(floor.(Int,PSY.get_max_active_power(g_s)),1,N);
            end  
        else
            gen_stor_charge_cap_array[idx,:] = fill.(floor.(Int,getfield(PSY.get_input_active_power_limits(PSY.get_storage(g_s)), :max)),1,N);
            gen_stor_discharge_cap_array[idx,:] = fill.(floor.(Int,getfield(PSY.get_output_active_power_limits(PSY.get_storage(g_s)), :max)),1,N);
            gen_stor_enrgy_cap_array[idx,:] = fill.(floor.(Int,getfield(PSY.get_state_of_charge_limits(PSY.get_storage(g_s)), :max)),1,N); 
            gen_stor_gridinj_cap_array[idx,:] = fill.(floor.(Int,PSY.getfield(PSY.get_output_active_power_limits(g_s), :max)),1,N);
            
            if (PSY.has_time_series(PSY.get_renewable_unit(g_s)) && ("max_active_power" in PSY.get_name.(PSY.get_time_series_multiple(PSY.get_renewable_unit(g_s)))))
                gen_stor_inflow_array[idx,:] = floor.(Int,get_forecast_values(get_first_ts(PSY.get_time_series_multiple(PSY.get_renewable_unit(g_s), name = "max_active_power")))
                                               *PSY.get_max_active_power(PSY.get_renewable_unit(g_s))); 
            else
                gen_stor_inflow_array[idx,:] = fill.(floor.(Int,PSY.get_max_active_power(PSY.get_renewable_unit(g_s))),1,N); 
            end
        end
        
        if (~outage_flag)
            if (typeof(g_s) ==PSY.HydroEnergyReservoir)
                p_m = string(PSY.get_prime_mover_type(g_s))
                p_m_idx = findall(x -> x == p_m, getfield.(outage_values,:prime_mover))

                temp_cap = floor(Int,PSY.get_max_active_power(g_s))
                
                if (length(p_m_idx)>1)
                    for (x,y) in zip(p_m_idx,getfield.(outage_values[p_m_idx],:capacity))
                        temp=0
                        if (temp<temp_cap<y)
                            gen_idx = x

                            f_or = getfield(outage_values[gen_idx],:FOR)
                            mttr_hr = getfield(outage_values[gen_idx],:MTTR)

                            (λ,μ) = outage_to_rate(f_or,mttr_hr)
                            break
                        else
                            temp = y
                        end
                    end
                end
            else
                @warn "No outage information is available for $(PSY.get_name(g_s)) of type $(gen_stor_categories[idx]). Using nominal outage and recovery probabilities for this generator."
                λ = 0.0;
                μ = 1.0;
            end

        else
            ext = PSY.get_ext(g_s)
            if (!(haskey(ext,"outage_probability") && haskey(ext,"recovery_probability")))
                @warn "No outage information is available in ext field of $(PSY.get_name(g_s)) of type $(gen_stor_categories[idx]). Using nominal outage and recovery probabilities for this generator."
                λ = 0.0;
                μ = 1.0;
            else
                λ = ext["outage_probability"];
                μ = ext["recovery_probability"];
            end
        end
        
        λ_genstors[idx,:] = fill.(λ,1,N); 
        μ_genstors[idx,:] = fill.(μ,1,N);
    end
    
    gen_stor_gridwdr_cap_array = zeros(Int64,n_genstors, N); # Not currently available/ defined in PowerSystems
    gen_stor_charge_eff = ones(n_genstors,N);                # Not currently available/ defined in PowerSystems
    gen_stor_discharge_eff = ones(n_genstors,N);             # Not currently available/ defined in PowerSystems
    gen_stor_cryovr_eff = ones(n_genstors,N);                # Not currently available/ defined in PowerSystems

    
    new_gen_stors = PRAS.GeneratorStorages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(gen_stor_names,get_generator_category.(gen_stor),
                                                    gen_stor_charge_cap_array, gen_stor_discharge_cap_array, gen_stor_enrgy_cap_array,
                                                    gen_stor_charge_eff, gen_stor_discharge_eff, gen_stor_cryovr_eff,
                                                    gen_stor_inflow_array, gen_stor_gridwdr_cap_array, gen_stor_gridinj_cap_array,
                                                    λ_genstors, μ_genstors);

    #######################################################
    # PRAS SystemModel
    #######################################################
    if (system_model=="Zonal")
        #######################################################
        # PRAS Lines 
        #######################################################
        @info "Collecting all inter regional lines in PSY System..."

        lines = availability_flag ? 
        collect(PSY.get_components(x -> (typeof(x) ∉ TransformerTypes && PSY.get_available(x)), PSY.Branch, sys)) :
        collect(PSY.get_components(x -> (typeof(x) ∉ TransformerTypes), PSY.Branch, sys));

        #######################################################
        # Inter-Regional Line Processing
        #######################################################
        regional_lines = filter(x -> ~(x.arc.from.area.name == x.arc.to.area.name),lines);
        sorted_lines, interface_reg_idxs, interface_line_idxs = get_sorted_lines(regional_lines, region_names);
        new_lines, new_interfaces = make_pras_interfaces(sorted_lines, interface_reg_idxs, interface_line_idxs,N);
    
        pras_system = PRAS.SystemModel(new_regions, new_interfaces, new_generators, region_gen_idxs, new_storage, region_stor_idxs, new_gen_stors,
                          region_genstor_idxs, new_lines,interface_line_idxs,my_timestamps);

        @info "Successfully built a PRAS $(system_model) system of type $(typeof(pras_system))."

        if (pras_sys_exp_loc !== nothing)
            if ~(isprasfile(pras_sys_exp_loc))
                error("PRAS System export location should be a .pras file. $(pras_sys_exp_loc) is not a valid location.")
            else
                PRAS.savemodel(pras_system,pras_sys_exp_loc, string_length =100, verbose = true, compression_level = 9)
                @info "PRAS System exported can be found here : $(pras_sys_exp_loc)"
            end
        end
    
        return pras_system
    
    elseif (system_model =="Copper Plate")
        load_vector = vec(sum(region_load,dims=1));
        pras_system = PRAS.SystemModel(new_generators, new_storage, new_gen_stors, my_timestamps, load_vector);

        @info "Successfully built a PRAS $(system_model) system of type $(typeof(pras_system))."

        if (pras_sys_exp_loc !== nothing)
            if ~(isprasfile(pras_sys_exp_loc))
                error("PRAS System export location should be a .pras file. $(pras_sys_exp_loc) is not a valid location.")
            else
                PRAS.savemodel(pras_system,pras_sys_exp_loc, string_length =100, verbose = true, compression_level = 9)
                @info "PRAS System exported can be found here : $(pras_sys_exp_loc)"
            end
        end
    
        return pras_system
    else
        error("Unrecognized SystemModel; Please specify correctly if SystemModel is Single-Node or Zonal.")
    end
end

#######################################################
# Main Function to make the PRAS System
#######################################################
function make_pras_system(sys_location::String;
                          system_model::Union{Nothing, String} = nothing,aggregation::Union{Nothing, String} = nothing,
                          period_of_interest::Union{Nothing, UnitRange} = nothing,outage_flag=true,lump_pv_wind_gens=false,availability_flag=true, 
                          outage_csv_location::Union{Nothing, String} = nothing, pras_sys_exp_loc::Union{Nothing, String} = nothing)

    @info "Running checks on the System location provided ..."
    runchecks(sys_location)
    
    @info "The PowerSystems System is being de-serialized from the System JSON ..."
    sys = 
    try
        PSY.System(sys_location;time_series_read_only = true,runchecks = false);
    catch
        error("The PSY System could not be de-serialized using the location of JSON provided. Please check the location and make sure you have permission to access time_series_storage.h5")
    end

    make_pras_system(sys,system_model = system_model,aggregation = aggregation,period_of_interest = period_of_interest,
                     outage_flag = outage_flag,lump_pv_wind_gens = lump_pv_wind_gens,availability_flag = availability_flag, 
                     outage_csv_location = outage_csv_location, pras_sys_exp_loc = pras_sys_exp_loc) 
end

