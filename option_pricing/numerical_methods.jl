# Solve the Black-Scholes PDE using numerical schemes.
import LinearAlgebra

function option_payoff(S_array, K; option_type="call")
    if option_type == "call"
        max.(S_array .- K, 0.0)
    elseif option_type == "put"
        max.(K .- S_array, 0.0)
    end
end

function boundary_conditions(S, K, T, rf, t_i; boundary="lower", option_type="call")
    if option_type == "call"
        if boundary == "lower"
            return 0
        elseif boundary == "upper"
            return S .- K * exp.(-rf .* (T .- t_i)) 
        end
    elseif option_type == "put"
        if boundary == "lower"
            return K * exp.(-rf .* (T .- t_i)) 
        elseif boundary == "upper"
            return 0
        end
    end
    return NaN
end

function explicit_scheme()
end

function implicit_scheme(K, T, σ, rf; option_type="call", n_time=50, n_space=500)
    S = K * 4
    ds = S / n_space
    dt = T / n_time

    s_j = collect(range(0, S, length=n_space+1))
    t_i = collect(range(0, T, length=n_time+1))

    s_idx = collect(range(0, n_space))

    a_j = (0.5 * rf * dt .* s_idx) - (0.5 * dt * (σ^2) .* (s_idx .^ 2)) 
    b_j = (1 + dt * rf) .+ (dt * (σ^2) .* (s_idx .^ 2))
    c_j = (-0.5 * rf * dt .* s_idx) - (0.5 * dt * (σ^2) .* (s_idx .^ 2)) 
    D = LinearAlgebra.Tridiagonal(a_j[3: n_space], b_j[2:n_space], c_j[2:n_space-1])

    F = zeros((n_time+1, n_space+1))
    F[n_time+1, :] = option_payoff(s_j, K; option_type=option_type)
    F[:, 1] .= boundary_conditions(S, K, T, rf, t_i; boundary="lower", option_type=option_type)
    F[:, n_space+1] .= boundary_conditions(S, K, T, rf, t_i; boundary="upper", option_type=option_type)

    for t = n_time:-1:1
        offset = zeros(n_space-1)
        offset[1] = c_j[1] * F[t+1, 1]
        offset[n_space-1] = a_j[n_space+1] * F[t+1, n_space+1]

        F[t, 2:n_space] = D \ (F[t+1, 2:n_space] - offset)
    end

    return F
end
