let

    net = MetNets.toy_model()
    M, N = size(net)
    net = MetEP.fva_preprocess(net)
    objider = MetNets.TOY_MODEL_BIOMASS_IDER
    objidx = MetNets.rxnindex(net, objider)

    
    # maxent ep
    epout = MetEP.maxent_ep(net; verbose = false)
    beta_vec = zeros(N)
    Δbeta = 1.05
    beta_vec[objidx] = 10
    for it in 1:500
        epout_ = try
            beta_vec[objidx] *= Δbeta
            MetEP.maxent_ep(net; 
                beta_vec, verbose = false, 
                solution = epout
            )
        catch err;
            @warn("Error", err)
            break
        end
        !MetEP.is_converged(epout_) && break
        epout = epout_
    end

    # fbaout
    fbaout = MetLP.fba!(net, objidx)

    ep_objval = MetNets.av(net, epout, objidx)
    fba_objval = MetNets.av(net, fbaout, objidx)

    @show ep_objval
    @show fba_objval

    ep_av = MetNets.av(epout)
    fba_av = MetNets.av(fbaout)

    atol = maximum(abs, fba_av) * 0.1
    @test all(isapprox.(ep_av, fba_av; atol))

end
