"""
lognormal(x, μ, σ)

```math
\\frac{(x-μ)^2}{2σ^2}
```
"""
function lognormal(x, μ, σ)
  if σ == Inf
    return 0
  end
  if σ ≈ 0
    return Inf
  end
  (x-μ)^2/2/σ^2
end

"""
logpoisson(n, x...)

  λ > 0
  λ = sum([x...])
  λ - n*log(λ) + n*log(n) - n
"""
function logpoisson(n, x...)
  λ = sum([x...])
  if λ <= 0
    return Inf
  end
  λ - n*log(λ) + n*log(n) - n
end
