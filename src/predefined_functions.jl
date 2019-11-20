function lognormal(x, μ, σ)
  (x-μ)^2/2/σ^2
end

function logpoisson(n, x...)
  λ = sum([x...])
  λ - n*log(λ) + n*log(n) - n
end
