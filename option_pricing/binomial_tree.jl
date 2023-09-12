# Is this a comment? Seems like it.
# What are the inpus to the Black-Scholes formula?
# S; K; rf; σ; T
# Steps:
# 1. Choose δ_t; n_t = T / δ_t
# 2. Set u = exp(σ * sqrt(δ_t)); d = 1/u
# 3. Set p_u = \frac{exp(r * δ_t) - d}{u - d}; p_d = 1 - p_u;
# 4. Create (n_t + 1) x (n_t + 1) matrix.
# 5. Update evolution of stock price.
# 6. Set payoff for t = T and backpropagate.

function option_payoff(stock_prices, strike_price; option_type="call")
    if option_type == "call"
        max.(stock_prices .- strike_price, 0.0)
    elseif option_type == "put"
        max.(strike_price .- stock_prices, 0.0)
    end
end

function binomial_tree_price(S, K, T, σ, rf; option_type="call", δ_t=0.001)
    n_t = convert(Int64, T / δ_t)
    u = exp(σ * sqrt(δ_t))
    d = exp(-σ * sqrt(δ_t))
    p_u = (exp(rf * δ_t) - d) / (u - d)
    p_d = 1 - p_u

    stock_prices = S * broadcast(^, u, [ (i > j) ? NaN : (j-i) - (i-1) for i=1:(n_t+1), j=1:(n_t+1) ])

    option_prices = copy(stock_prices)
    option_prices[:, n_t+1] = option_payoff(option_prices[:, n_t+1], K; option_type=option_type)
    for t = n_t:-1:1
        option_prices[1:t, t] = exp(-rf * δ_t) * (p_u * option_prices[1:t, t+1] + p_d * option_prices[2:t+1, t+1])
    end

    option_prices[1, 1]
end

call_price = binomial_tree_price(300, 250, 2, 0.15, 0.03; option_type = "call")
put_price  = binomial_tree_price(300, 250, 2, 0.15, 0.03; option_type = "put")
println(call_price, " ", put_price, " ", abs(call_price - put_price - (300 - 250 * exp(-0.06))) < 1e-6)
