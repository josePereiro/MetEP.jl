function get_join(epmodel::EPModel)
    @extract epmodel: epmat scalefact

    if epmat isa EPMat
        Σ = deepcopy(epmat.Σ)
        v̄ = deepcopy(epmat.v)
        rmul!(v̄, scalefact)
        rmul!(Σ, scalefact^2)
        return MvNormal(v̄, Σ)
    else epmat isa EPMatT0
        # TODO rebuild the full join
        # Now v̄ are just the independent and are unsorted
        Σ = deepcopy(epmat.Σi)
        v̄ = deepcopy(epmat.vi)
        rmul!(v̄, scalefact)
        rmul!(Σ, scalefact^2)
        return MvNormal(v̄, Σ)
    end
end