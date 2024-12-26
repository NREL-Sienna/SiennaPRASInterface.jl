"""
    get_sorted_region_tuples(lines::Vector{PSY.Branch}, region_names::Vector{String})

Get sorted (reg_from, reg_to) tuples of inter-regional lines.
"""
function get_sorted_region_tuples(lines::Vector{PSY.Branch}, region_names::Vector{String})
    region_idxs = Dict(name => idx for (idx, name) in enumerate(region_names))

    line_from_to_reg_idxs = similar(lines, Tuple{Int, Int})

    for (l, line) in enumerate(lines)
        from_name = PSY.get_name(PSY.get_area(PSY.get_from_bus(line)))
        to_name = PSY.get_name(PSY.get_area(PSY.get_to_bus(line)))

        from_idx = region_idxs[from_name]
        to_idx = region_idxs[to_name]

        line_from_to_reg_idxs[l] =
            from_idx < to_idx ? (from_idx, to_idx) : (to_idx, from_idx)
    end

    return line_from_to_reg_idxs
end

"""
    get_sorted_lines(lines::Vector{PSY.Branch}, region_names::Vector{String})

Get sorted lines, interface region indices, and interface line indices.

# Arguments

  - `lines::Vector{PSY.Branch}`: Lines
  - `region_names::Vector{String}`: Region names

# Returns

  - `sorted_lines::Vector{PSY.Branch}`: Sorted lines
  - `interface_reg_idxs::Vector{Tuple{Int, Int}}`: Interface region indices
  - `interface_line_idxs::Vector{UnitRange{Int}}`: Interface line indices
"""
function get_sorted_lines(lines::Vector{PSY.Branch}, region_names::Vector{String})
    line_from_to_reg_idxs = get_sorted_region_tuples(lines, region_names)
    line_ordering = sortperm(line_from_to_reg_idxs)

    sorted_lines = lines[line_ordering]
    sorted_from_to_reg_idxs = line_from_to_reg_idxs[line_ordering]
    interface_reg_idxs = unique(sorted_from_to_reg_idxs)

    # Ref tells Julia to use interfaces as Vector, only broadcasting over
    # lines_sorted
    interface_line_idxs = searchsorted.(Ref(sorted_from_to_reg_idxs), interface_reg_idxs)

    return sorted_lines, interface_reg_idxs, interface_line_idxs
end

"""
    make_pras_interfaces(
        sorted_lines::Vector{PSY.Branch},
        interface_reg_idxs::Vector{Tuple{Int64, Int64}},
        interface_line_idxs::Vector{UnitRange{Int64}},
        s2p_meta::S2P_metadata,
    )

Converts PSY branches and interaces indices into PRAS Lines and Interfaces.

# Returns

  - `new_lines::PRASCore.Lines`: PRAS Lines
  - `new_interfaces::PRASCore.Interfaces`: PRAS Interfaces
"""
function make_pras_interfaces(
    sorted_lines::Vector{PSY.Branch},
    interface_reg_idxs::Vector{Tuple{Int64, Int64}},
    interface_line_idxs::Vector{UnitRange{Int64}},
    s2p_meta::S2P_metadata,
)
    num_interfaces = length(interface_reg_idxs)
    interface_regions_from = first.(interface_reg_idxs)
    interface_regions_to = last.(interface_reg_idxs)
    num_lines = length(sorted_lines)

    # Lines
    line_names = PSY.get_name.(sorted_lines)
    line_cats = string.(typeof.(sorted_lines))

    line_forward_cap = Matrix{Int64}(undef, num_lines, s2p_meta.N)
    line_backward_cap = Matrix{Int64}(undef, num_lines, s2p_meta.N)
    line_λ = Matrix{Float64}(undef, num_lines, s2p_meta.N) # Not currently available/ defined in PowerSystems
    line_μ = Matrix{Float64}(undef, num_lines, s2p_meta.N) # Not currently available/ defined in PowerSystems

    for i in 1:num_lines
        line_forward_cap[i, :] =
            fill.(
                floor.(Int, getfield(line_rating(sorted_lines[i]), :forward_capacity)),
                1,
                s2p_meta.N,
            )
        line_backward_cap[i, :] =
            fill.(
                floor.(Int, getfield(line_rating(sorted_lines[i]), :backward_capacity)),
                1,
                s2p_meta.N,
            )

        line_λ[i, :], line_μ[i, :] = get_outage_time_series_data(sorted_lines[i], s2p_meta)
    end

    new_lines = PRASCore.Lines{
        s2p_meta.N,
        s2p_meta.pras_timestep,
        s2p_meta.pras_resolution,
        PRASCore.MW,
    }(
        line_names,
        line_cats,
        line_forward_cap,
        line_backward_cap,
        line_λ,
        line_μ,
    )
    interface_forward_capacity_array = Matrix{Int64}(undef, num_interfaces, s2p_meta.N)
    interface_backward_capacity_array = Matrix{Int64}(undef, num_interfaces, s2p_meta.N)

    for i in 1:num_interfaces
        interface_forward_capacity_array[i, :] =
            sum(line_forward_cap[interface_line_idxs[i], :], dims=1)
        interface_backward_capacity_array[i, :] =
            sum(line_backward_cap[interface_line_idxs[i], :], dims=1)
    end

    new_interfaces = PRASCore.Interfaces{s2p_meta.N, PRASCore.MW}(
        interface_regions_from,
        interface_regions_to,
        interface_forward_capacity_array,
        interface_backward_capacity_array,
    )

    return new_lines, new_interfaces
end
