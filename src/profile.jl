function profile!(name, results; stop=5, step=1.0, prior=nothing)
  idx = 0
  for (i, p) in enumerate(results.model.params)
    if p.name == name
      idx = i
    end
  end
  if idx == 0
    return 0, 0
  end
  #component = results.model.component_dict[name]
  param = results.model.params[idx]
  ## Go left
  x_left = []
  nll_left = []
  nlow = copy(results.lower_bounds)
  nhigh = copy(results.upper_bounds)
  p0 = copy(results.min_parameters)
  val = 0
  xeval = p0[idx]
  count = 0
  while (val-results.min_objective) < stop
    nlow[idx] = xeval
    nhigh[idx] = xeval
    results.opt.lower_bounds = nlow
    results.opt.upper_bounds = nhigh
    p0[idx] = xeval
    init_step = p0./10.0
    init_step[idx] = xeval * 1e-6
    results.opt.initial_step=init_step
    val, minx, ret = optimize!(results.opt, p0)
    likelihood_value = exp(-(val-results.min_objective))
    if (likelihood_value == Inf) || (likelihood_value < 0)
      break
    end
    push!(nll_left, likelihood_value )
    push!(x_left, xeval)
    xeval -= step
    count += 1
    if count > 1000
      println("trials > 1000 -- ", val)
      @show nll_left
      break
    end
  end
  ## Go right
  x_right = []
  nll_right = []
  nlow = copy(results.lower_bounds)
  nhigh = copy(results.upper_bounds)
  p0 = copy(results.min_parameters)
  val = 0
  xeval = p0[idx]
  count = 0
  while (val-results.min_objective) < stop
    nlow[idx] = xeval
    nhigh[idx] = xeval
    results.opt.lower_bounds = nlow
    results.opt.upper_bounds = nhigh
    p0[idx] = xeval
    init_step = p0./10.0
    init_step[idx] = xeval * 1e-6
    results.opt.initial_step=init_step
    val, minx, ret = optimize!(results.opt, p0)
    likelihood_value = exp(-(val-results.min_objective))
    if (likelihood_value == Inf) || (likelihood_value < 0)
      break
    end
    push!(nll_right, likelihood_value )
    push!(x_right, xeval)
    xeval += step
    count += 1
    if count > 1000
      println("trials > 1000 -- ", val)
      @show nll_right
      break
    end
  end
  x = append!(reverse(x_left), x_right)
  nll = append!(reverse(nll_left), nll_right)
  ## update the component with the likelihood results
  if prior != nothing
    nll = nll .* prior.(x)
  end
  param.likelihood_x = x
  param.likelihood_y = nll
  return x, nll
end

function compute_profiled_uncertainties!(results 
                                         ; CL=nothing, σ=nothing, mode="FC")
  if CL != nothing
    α = CL
  elseif σ != nothing
    α = erf(σ/sqrt(2))
  else
    α = 0.90
  end
  for param in results.model.params
    profile!(param.name, results)
    uncertainty!(param.name, results; CL=α, mode=mode)
  end
end

function uncertainty!(name, results; 
                      CL=nothing, σ=nothing, mode="FC", prior=nothing)
  if CL != nothing
    α = CL
  elseif σ != nothing
    α = erf(σ/sqrt(2))
  else
    α = erf(1/sqrt(2)) # 1 Sigma default
  end
  idx = 0
  for (i, p) in enumerate(results.model.params)
    if p.name == name
      idx = i
    end
  end
  param = results.model.params[idx]
  x, y = param.likelihood_x, param.likelihood_y
  # Apply the arbitrary prior function f(x)
  if prior != nothing
    y = y .* prior.(x)
  end
  # Mode can be FC (Feldman-Cousins), Mode-centered, mean-centered, left, right, etc.
  # Return cumulative, x_low, x_high
  y = y / sum(y)
  cumulative = cumsum(y)
  β = 1 - α
  if mode == "FC"
    cy = sortperm(y; rev=true)
    cum_y = cumsum(y[cy])
    cum_x = x[cy]
    region = cum_x[ cum_y .<= α ]
    left = minimum(region)
    right = maximum(region)
  end
  param.low, param.high = left, right
  left, right
end
