function run_maxent_ep_test(net, objider; checkall = false)
    
    # setup
    net = MetEP.fva_preprocess(net)
    net = MetNets.force_dims(net)
    M, N = size(net)
    objidx = MetNets.rxnindex(net, objider)

    # fbaout
    fbaout = MetLP.fba(net, objidx)
    fba_objval = MetNets.av(net, fbaout, objidx)

    # maxent ep
    epout = MetEP.maxent_ep(net; verbose = false)
    beta_vec = zeros(N)
    Δbeta = 1.1
    beta_vec[objidx] = 10
    ep_objval = MetNets.av(net, epout, objidx)
    
    obj_th = 1e-4
    prog = ProgressThresh(obj_th)
    for it in 1:5000
        epout_ = try
            beta_vec[objidx] *= Δbeta
            MetEP.maxent_ep(net; 
                beta_vec, verbose = false, 
                epsconv = 1e-4,
                solution = epout
            )
        catch err;
            beta = beta_vec[objidx]
            @warn("Error", it, beta, err)
            break
        end
        # !MetEP.is_converged(epout_) && break
        epout = epout_
        
        ep_objval = MetNets.av(net, epout, objidx)
        update!(prog, abs(ep_objval - fba_objval) / fba_objval)
        isapprox(ep_objval, fba_objval; rtol = obj_th) && break
    end
    finish!(prog)
    
    @show ep_objval
    @show fba_objval
    @test isapprox(ep_objval, fba_objval; rtol = obj_th)
    
    if checkall
        ep_av = MetNets.av(epout)
        fba_av = MetNets.av(fbaout)
        @show maximum(abs, ep_av - fba_av)

        rtol = 1e-2
        @test all(isapprox.(ep_av, fba_av; rtol = rtol))
    end

end

let
    @info("maxent_ep Toy Model")
    net = MetNets.toy_model()
    objider = MetNets.TOY_MODEL_BIOMASS_IDER
    checkall = true
    run_maxent_ep_test(net, objider; checkall)
    println()

    @info("maxent_ep EColi Core Model")
    net = MetNets.ecoli_core_model()
    objider = MetNets.ECOLI_MODEL_BIOMASS_IDER
    checkall = false # fba (max biomass) is degenerated
    run_maxent_ep_test(net, objider; checkall)
end