import Distributions

function black_scholes_price(S, K, T, σ, rf; option_type="call")
    d1 = (log(S/K) + (rf + σ^2/2)*T) / (σ * sqrt(T))
    d2 = d1 - σ*sqrt(T)
    N = Distributions.Normal()
    if option_type == "call"
        (Distributions.cdf.(N, d1) * S) - (Distributions.cdf.(N, d2) * K * exp(-rf * T))
    else
        (Distributions.cdf.(N, -d2) * K * exp(-rf * T)) - (Distributions.cdf.(N, -d1) * S)
    end
end

call_price = black_scholes_price(300, 250, 2, 0.15, 0.03; option_type = "call")
put_price  = black_scholes_price(300, 250, 2, 0.15, 0.03; option_type = "put")
println(call_price, " ", put_price, " ", abs(call_price - put_price - (300 - 250 * exp(-0.06))) < 1e-6)
