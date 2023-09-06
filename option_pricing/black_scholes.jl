using Random, Distributions

function black_scholes_price(S, K, T, σ, rf, option_type="call")
    d1 = (log(S/K) + (rf + σ^2/2)*T) / (σ * sqrt(T))
    d2 = d1 - σ*sqrt(T)
    N = Normal()
    if option_type == "call"
        (cdf.(N, d1) * S) - (cdf.(N, d2) * K * exp(-rf * T))
    else
        (cdf.(N, -d2) * K * exp(-rf * T)) - (cdf.(N, -d1) * S)
    end
end

call_price = black_scholes_price(300, 250, 2, 0.15, 0.03, "call")
put_price  = black_scholes_price(300, 250, 2, 0.15, 0.03, "put")
println(call_price, " ", put_price, " ", call_price - put_price)
