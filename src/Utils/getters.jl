import MetNets: av, va

# MetState interface
av(out::EPOut) = out.av
va(out::EPOut) = out.va

# EPOut
μ(out::EPOut) = out.μ
σ(out::EPOut) = out.σ