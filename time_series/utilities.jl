module Utilities

function lag(x, n::Int)
    padded_vals = [NaN * ones(n); x[1: length(x) - n]]
    return padded_vals
end

end
