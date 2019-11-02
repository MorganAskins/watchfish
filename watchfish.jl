module watchfish

export Component, Model, Results, add_component!, run_fish, pretty_results
using NLopt
using DataFrames
using DataStructures

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

function run_fish(m::Model, data)
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

end # End module
