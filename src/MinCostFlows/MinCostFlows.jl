module MinCostFlows

# Efficiently solves the min-cost network flow problem using the Relaxation dual ascent method of Bertsekas (1985).

export FlowProblem,
    solveflows!,
    updateinjection!,
    updateflowlimit!,
    updateflowcost!,
    flows,
    costs,
    limits,
    injections,
    prices

include("lists.jl")
include("FlowProblem.jl")
include("solveflows.jl")
include("solveflows_singlenode.jl")
include("solveflows_multinode.jl")
include("update.jl")
include("utils.jl")

end
