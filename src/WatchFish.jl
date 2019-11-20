module WatchFish

using NLopt
using DataFrames
using DataStructures
using SpecialFunctions

abstract type Model end

# Primary definitions
include("predefined_functions.jl")
include("variable.jl")
include("likelihood.jl")
include("model.jl")
include("plotter.jl")
# Kit definitions

# Move me (below)
function optimize_model!(m::Model, n::NLogLikelihood; 
                         lower_bounds=nothing, upper_bounds=nothing)
  function objective(x::Vector, grad::Vector)
    if length(grad)>0
      grad = 2x
    end
    n.objective(x)
  end
  opt = Opt(:LN_SBPLX, n.numparams)
  opt.ftol_rel = 1e-4
  if lower_bounds != nothing
    opt.lower_bounds = lower_bounds
  end
  if upper_bounds != nothing
    opt.upper_bounds = upper_bounds
  end
  opt.min_objective = objective
  p0 = [v.init for v in n.parameters]
  (minf, minx, ret) = optimize!(opt, p0)
  numevals = opt.numevals
  for (i,p) in enumerate(m.params)
    p.fit = minx[i]
  end
  Results( opt, m, copy(minf), copy(minx), numevals,
           copy(opt.lower_bounds), copy(opt.upper_bounds) )
end

mutable struct CountingExperiment <: Model
  params
  pdflist
  nll::NLogLikelihood
  counts::Constant
  function CountingExperiment()
    new( [], [], NLogLikelihood(),
        Constant("counts", 0.0) )
  end
end

function add_component!( m::CountingExperiment, name::String, μ; σ=Inf )
  p = Parameter( name; init=μ )
  push!(m.params, p)
  if σ != Inf
    a = Constant( name*"_mu", μ )
    b = Constant( name*"_sig", σ )
    nlp = NLogPDF( "lognormal", p, a, b )
    push!(m.pdflist, nlp)
  end
end

function set_counts!( m::CountingExperiment, counts::Number )
  m.counts.init = Float64(counts)
end

function minimize!( m::CountingExperiment )
  poisson_pdf = NLogPDF("logpoisson", m.counts, (m.params)...)
  likelihood = NLogLikelihood([m.pdflist..., poisson_pdf])
  add_likelihood!( m, likelihood )
  optimize_model!( m, likelihood; 
                   lower_bounds=[0.0 for a in 1:length(m.params)] )
end

function add_likelihood!( m::Model, nll::NLogLikelihood )
  m.nll = nll
  m.params = nll.parameters
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


function constraint(x, μ, σ)
  if σ == Inf
    return 0
  end
  if σ ≈ 0
    return Inf
  end
  (x-μ)^2/2/σ^2
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
#                 Estimate = Float64[],
#                 Constraint = Float64[],
                 Fit = Float64[],
                 Interval_Low = Float64[],
                 Interval_High = Float64[],
                )
  for v in r.model.params
    push!( df,
          (
           v.name,
#           v.μ,
#           v.σ,
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
  ## Go right
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

function add_dataset(name, df::DataFrame)
  eval(:( $name = $df ) )
end


# todo, move macros
macro addfunction(f)
  :($f)
end

macro dataset( name, df )
  :( $(esc(name)) = $df )
end

## Plotting tools

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
