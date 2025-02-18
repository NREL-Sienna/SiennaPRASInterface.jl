function test_names_equal(
    x::Vector{String},
    y::Vector{String};
    left="Left-side",
    right="Right-side",
)
    x = sort(x)
    y = sort(y)
    idx = 1
    while idx <= length(x) && idx <= length(y)
        if x[idx] == y[idx]
            idx += 1
        elseif x[idx] != y[idx]
            prev_idx = max(idx - 2, 1)
            next_idx = min(idx + 2, length(x))
            @error "$left $(x[prev_idx:next_idx]) != $right $(y[prev_idx:next_idx])"
            return false
        end
    end
    if idx > length(x) && idx > length(y)
        return true
    end
    prev_idx = max(idx - 2, 1)
    next_idx = min(idx + 2, length(x))
    comparison = idx <= length(x) ? "more" : "less"
    @error "$left has $comparison elements $(x[prev_idx:next_idx]) than $right $(y[prev_idx:next_idx])"
    return false
end

function array_all_equal(x::AbstractVector{T}, v::T) where {T}
    for (i, xi) in enumerate(x)
        if !isapprox(xi, v)
            @error "array[$i]: $xi != $v"
            return false
        end
    end
    return true
end
