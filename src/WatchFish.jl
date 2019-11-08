module WatchFish

using NLopt
using DataFrames
using DataStructures
using SpecialFunctions

# Note, should this really be mutable?
mutable struct Component
  name::String
  μ::Float64
  σ::Float64
  fit::Float64
  low::Float64
  high::Float64
  likelihood_x::Array{Float64}
  likelihood_y::Array{Float64}
  function Component(name::String, μ::Float64, σ::Float64)
    new(name, μ, σ, 0.0, 0.0, 0.0,
        Array{Float64}(undef, 0),
        Array{Float64}(undef, 0),
       )
  end
end

mutable struct Model
  component_dict::OrderedDict{String, Component}
  dims::Int64
  function Model()
    new( 
        OrderedDict{String, Component}(), 
        0, 
       ) 
  end
end

mutable struct Results
  opt
  model::Model
  min_objective::Float64
  min_parameters::Array{Float64}
  iterations::Int64
  lower_bounds::Array{Float64}
  upper_bounds::Array{Float64}
end

function add_component!(m::Model, name::String, μ; σ=Inf)
  get!(m.component_dict, name, Component(name, μ, σ))
  m.dims = length(m.component_dict)
end

function constraint(x, μ, σ)
  if σ == Inf
    return 0
  end
  if σ ≈ 0
    return Inf
  end
  (x-μ)^2/2/σ^2
end

function run_fish!(m::Model, data)
  N  = data
  β  = Array{Float64}(undef, 0)
  σ  = Array{Float64}(undef, 0)
  p0 = Array{Float64}(undef, 0)
  name = Array{String}(undef, 0)
  for (k, v) in m.component_dict
    push!(β, v.μ)
    push!(p0, v.μ)
    push!(σ, v.σ)
    push!(name, k)
  end

  function objective(x::Vector, grad::Vector)
    if length(grad)>0
      grad = 2x
    end
    λ = sum(x)
    if λ <= 0
      return Inf
    end
    ξ = sum( constraint.(x, β, σ) )
    λ - N*log(λ) + N*log(N) - N + ξ
  end
  opt = Opt(:LN_SBPLX, m.dims)
  opt.ftol_rel = 1e-4
  lb = [0.0 for a in 1:length(β)]
  opt.lower_bounds = lb
  opt.min_objective = objective
  (minf, minx, ret) = optimize!(opt, p0)
  numevals = opt.numevals
  for (idx, (k, v)) in enumerate(m.component_dict)
    v.fit = minx[idx]
  end
  Results( opt, m, copy(minf), copy(minx), numevals,
          copy(opt.lower_bounds), copy(opt.upper_bounds) )
end

## Show the results
function Base.show(io::IO, m::Results)
  compact = get(io, :compact, false)
  if !compact
    println("Fit converged after $(m.iterations) iterations")
    println("Best fit at: -LnL = $(m.min_objective)")
    println("Values of $(m.min_parameters)")
  end
end

function pretty_results(r::Results)
  df = DataFrame(
                 Name = String[],
                 Estimate = Float64[],
                 Constraint = Float64[],
                 Fit = Float64[],
                 Interval_Low = Float64[],
                 Interval_High = Float64[],
                )
  for (k,v) in r.model.component_dict
    push!( df,
          (
           v.name,
           v.μ,
           v.σ,
           v.fit,
           v.low,
           v.high,
          )
         )
  end
  return df
end

function profile!(name, results; stop=5, step=1.0, prior=nothing)
  idx = 0
  for (i, (k,v)) in enumerate(results.model.component_dict)
    if k == name
      idx = i
    end
  end
  if idx == 0
    return 0, 0
  end
  component = results.model.component_dict[name]
  # Go left
  x_left = []
  nll_left = []
  nlow = copy(results.lower_bounds)
  nhigh = copy(results.upper_bounds)
  p0 = copy(results.min_parameters)
  val = 0
  xeval = p0[idx]
  while (val-results.min_objective) < stop
    nlow[idx] = xeval
    nhigh[idx] = xeval
    results.opt.lower_bounds = nlow
    results.opt.upper_bounds = nhigh
    p0[idx] = xeval
    val, minx, ret = optimize!(results.opt, p0)
    push!(nll_left, exp(-(val-results.min_objective)) )
    push!(x_left, xeval)
    xeval -= step
  end
  # Go right
  x_right = []
  nll_right = []
  nlow = copy(results.lower_bounds)
  nhigh = copy(results.upper_bounds)
  p0 = copy(results.min_parameters)
  val = 0
  xeval = p0[idx]
  while (val-results.min_objective) < stop
    nlow[idx] = xeval
    nhigh[idx] = xeval
    results.opt.lower_bounds = nlow
    results.opt.upper_bounds = nhigh
    p0[idx] = xeval
    val, minx, ret = optimize!(results.opt, p0)
    push!(nll_right, exp(-(val-results.min_objective)) )
    push!(x_right, xeval)
    xeval += step
  end
  x = append!(reverse(x_left), x_right)
  nll = append!(reverse(nll_left), nll_right)
  # update the component with the likelihood results
  if prior != nothing
    nll = nll .* prior.(x)
  end
  component.likelihood_x = x
  component.likelihood_y = nll
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

  for (k, v) in results.model.component_dict
    profile!(k, results)
    uncertainty!(k, results; CL=α, mode=mode)
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
  comp = results.model.component_dict[name]
  x, y = comp.likelihood_x, comp.likelihood_y
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
  comp.low, comp.high = left, right
  left, right
end

## Plotting tools
include("plotter.jl")

# WatchFish exports based on a blacklist; where functions
# which begin with "_" are not exported.
const _EXCLUDE_SYMBOLS=[Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__, all=true)
  sym_string = string(sym)
  if sym in _EXCLUDE_SYMBOLS || startswith(sym_string, "_")
    continue
  end
  if !(Base.isidentifier(sym) || (startswith(sym_string, "@") &&
                                  Base.isidentifier(sym_string[2:end])))
    continue
  end
  @eval export $sym
end

end # End module
