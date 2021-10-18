# A struct that contain all the data required for converging ep
struct EPModel{T<:Real}
    scalefact::T
    updatealg!::Function
    epfields::EPFields{T}
    epmat::AbstractEPMat
    alpha::T
    beta_vec::SparseVector{T, Int}
    stat::Dict # Just data about the ep run time
end

function EPModel(
        S::AbstractArray{T,2}, b::AbstractArray{T}, lb::AbstractArray{T}, ub::AbstractArray{T};
        alpha::Real=Inf,                                        # inverse temperature
        beta_vec::AbstractVector{T} = spzeros(T, size(S, 2)),   # maxent inverse temperature vector
        solution::Union{EPOut{T}, Nothing} = nothing,           # start from a solution
        expval = nothing                                        # fix posterior probability experimental values for std and mean
    ) where {T<:Real}

    # Some checks
    M, N = size(S)
    M > N && @warn("M = $M ≥ N = $N")
    any(lb .> ub) && error("lower bound fluxes > upper bound fluxes. Consider swapping lower and upper bounds")

    # The scalefactor is just the maximum absolute bound (lb or ub).
    scalefact = get_scalefactor(lb, ub)

    # Create EPFields. If a solution is not given, the EPfields will be fresh
    epfields::EPFields = isnothing(solution) ? EPFields(N, expval, eltype(S)) : 
        deepcopy(solution.sol) # preserve the original solution!

    # making a local copy to rescale
    lb, ub, b = copy.([lb, ub, b]) 

    #=
    Scale down μ, s, av, va of epfields and ub, lb and Y using the 
    previous computed scalefactor.
    If epfields is fresh, it only will have effect on ub, lb and Y
    =#
    scaleepfield!(inv(scalefact), epfields, ub, lb, b) # scaling fields in [0,1]

    epmat = (alpha < Inf) ? EPMat(S, b, lb, ub, alpha) : EPMatT0(S, b, lb, ub)

    # One iteration of EP
    updatealg! = alpha == Inf ? eponesweepT0! : eponesweep!

    beta_vec = prepare_beta_vec(epmat, beta_vec)

    return EPModel{T}(scalefact, updatealg!, epfields, epmat, alpha, beta_vec, Dict())

end

EPModel(metnet::MetNet; kwargs...) = EPModel(metnet.S, metnet.b, metnet.lb, metnet.ub; kwargs...)