function produce_epout!(epfields, epmat, alpha, scalefact, returnstatus, iter; drop_epfields)

    #= Scale back μ, s, av, va of epfields and lub, llb and Y =#
    scaleepfield!(scalefact, epfields)
    μ, σ = epfields.μ, epfields.s
    av, va = epfields.av, epfields.va
    if isinf(alpha)
        sperm = sortperm(epmat.usperm)
        μ, σ, av, va = μ[sperm], σ[sperm], av[sperm], va[sperm]
    end
    
    sol = drop_epfields ? nothing : epfields
    return  EPOut(μ, σ, av, va, sol, returnstatus, iter)
end

produce_epout!(epmodel, returnstatus = UNSET_STATUS, iter = -1; drop_epfields = false) = 
    produce_epout!(epmodel.epfields, epmodel.epmat, epmodel.alpha, epmodel.scalefact, returnstatus, iter; drop_epfields)
produce_epout(epmodel, returnstatus = UNSET_STATUS, iter = -1; drop_epfields = false) = 
    produce_epout!(deepcopy(epmodel.epfields), epmodel.epmat, epmodel.alpha, epmodel.scalefact, returnstatus, iter; drop_epfields)