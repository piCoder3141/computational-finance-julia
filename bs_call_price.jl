using Random, Distributions

function black_scholes_price(S, K, T, σ, rf)
    d1 = (log(S/K) + (rf + σ^2/2)*T) / (σ * sqrt(T))
    d2 = d1 - σ*sqrt(T)
    N = Normal()
    (cdf.(N, d1) * S) - (cdf.(N, d2) * K * exp(-rf * T))
end

println(black_scholes_price(300, 250, 2, 0.15, 0.03))
