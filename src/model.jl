abstract type Model end

mutable struct Results
  opt
  model::Model
  min_objective::Float64
  min_parameters::Array{Float64}
  iterations::Int64
  lower_bounds::Array{Float64}
  upper_bounds::Array{Float64}
end

## Generic model functions
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

function getparam(m::Model, name::String)
  for (i, p) in enumerate(m.params)
    if p.name == name
      return p
    end
  end
  return nothing
end

function add_likelihood!( m::Model, nll::NLogLikelihood )
  m.nll = nll
  m.params = nll.parameters
end
