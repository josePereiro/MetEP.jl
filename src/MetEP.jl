
#= The EP implementation is derived from:
Alfredo Braunstein, Anna Paola Muntoni, and Andrea Pagnani, “An Analytic Approximation of the Feasible 
Space of Metabolic Networks,” Nature Communications 8 (April 6, 2017): https://doi.org/10.1038/ncomms14915.
Package metabolicEP (https://github.com/anna-pa-m/Metabolic-EP) 
It was modified just for implement MaxEnt as described in the Appendix of:
Jorge Fernandez-de-Cossio-Diaz and Roberto Mulet, “Maximum Entropy and Population Heterogeneity in 
Continuous Cell Cultures.” PLOS Computational Biology 15, no. 2 (February 27, 2019): e1006823. 
https://doi.org/10.1371/journal.pcbi.1006823.
=#

module MetEP

    using SparseArrays
    using ProgressMeter

    import MetNets
    import MetNets: MetNet, AbstractMetState
    import ExtractMacro: @extract
    using LinearAlgebra
    import Printf: @printf
    import SpecialFunctions: erf
    import Distributions: Truncated, Normal, mean, var

    import MetLP
    import MetLP: fva_preprocess

    # TODO: Add Chemostat HR code
    
    # Types
    include("Types/AbstractEPMat.jl")
    include("Types/EPOut.jl")
    include("Types/EPAlgs.jl")
    include("Types/EPFields.jl")
    include("Types/EPMat.jl")
    include("Types/EPMatT0.jl")
    include("Types/EPModel.jl")
    
    # EP
    include("EP/compute_mom5d.jl")
    include("EP/converge_ep!.jl")
    include("EP/echelonize.jl")
    include("EP/epconverge.jl")
    include("EP/eponesweep.jl")
    include("EP/eponesweepT0.jl")
    include("EP/fast_maxent_ep.jl")
    include("EP/get_scalefactor.jl")
    include("EP/matchmom.jl")
    include("EP/maxent_ep.jl")
    include("EP/newav.jl")
    include("EP/newμs.jl")
    include("EP/prepare_beta_vec.jl")
    include("EP/prepareinput.jl")
    include("EP/produce_epout.jl")
    include("EP/Q_sigma.jl")
    include("EP/scaleepfield.jl")
    include("EP/status.jl")
    include("EP/toy_maxent_ep.jl")
    include("EP/utils.jl")

    # Utils
    include("Utils/summary.jl")
    include("Utils/getters.jl")

end