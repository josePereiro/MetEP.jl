const UNCONVERGED_STATUS = :unconverged
const CONVERGED_STATUS = :converged

is_converged(out::EPOut) = (out.status == CONVERGED_STATUS)
