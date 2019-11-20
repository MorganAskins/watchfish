function lognormal(x, μ, σ)
  if σ == Inf
    return 0
  end
  if σ ≈ 0
    return Inf
  end
  (x-μ)^2/2/σ^2
end

function logpoisson(n, x...)
  λ = sum([x...])
  if λ <= 0
    return Inf
  end
  λ - n*log(λ) + n*log(n) - n
end
