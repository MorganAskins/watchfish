"""
abstract type Model

Models hold the parameters of interest and the associated likelihood/cost
functions used in optimization. Individual models should be derived from
this base object.
"""
abstract type Model end

"""
Results

Created upon optimizing a model, contains the model and all of
its parameters as well as the optimized results found by
minimizing the cost function.
"""
mutable struct Results
  opt
  model::Model
  min_objective::Float64
  min_parameters::Array{Float64}
  iterations::Int64
  lower_bounds::Array{Float64}
  upper_bounds::Array{Float64}
end

"""
    optimize_model!(m::Model, nll::NLogLikelihood;
                    lower_bounds=nothing, upper_bounds=nothing)

Currently NLopt is used for minimizing the objective function.
optimize_model! produces the objective function from the given
NLogLikelihood which is then passed to NLopt to be minimized.
optimize_model! returns a set of Results.

TODO: Move lower_bounds and upper_bounds into options dict
"""
function optimize_model!(m::Model, nll::NLogLikelihood; 
                         lower_bounds=nothing, upper_bounds=nothing,
                         options=Dict() )
  function objective(x::Vector, grad::Vector)
    if length(grad)>0
      grad = 2x
    end
    nll.objective(x)
  end
  opt = Opt(:LN_SBPLX, nll.numparams)
  # A note here on tolerence. With nlopt one can specify four
  # types of tolerence for a stopping condition. Either on the
  # value of the objective function, or on the values of the parameters.
  # For both of these it can be relative to the scale of each, or
  # absolute. For simplicity I have forced in tolerance on the objective
  # function at an absolute scale. This is because I am assuming that
  # the function is a log-likelihood, and thus has statistical information
  # in its derivative, but may be shifting by a constant NOT a scale factor.
  p0 = [v.init for v in nll.parameters]
  opt.ftol_abs     = get(options, "ftol_abs", 1e-6)
  opt.ftol_rel     = get(options, "ftol_rel", 0)
  opt.xtol_rel     = get(options, "xtol_rel", 0)
  opt.xtol_abs     = get(options, "xtol_abs", 0)
  opt.maxeval      = get(options, "maxeval", 0)
  opt.maxtime      = get(options, "maxtime", 0)
  initial_step     = get(options, "initial_step", p0./10.0)
  initial_step[(initial_step .== 0)] .= 0.1
  opt.initial_step = initial_step
  opt.stopval      = get(options, "stopval", -Inf)
  opt.lower_bounds = nll.lower_bounds
  opt.upper_bounds = nll.upper_bounds

  opt.min_objective = objective
  try
    nll.objective(p0)
  catch e
    @error "Problems contructing objective function"
    println(e)
  end
  # Test before optimizing
  (minf, minx, ret) = optimize!(opt, p0)
  ## If we fit twice we can go coarse the first time and then re-evaluate
  initial_step = minx ./ 10.0
  initial_step[(initial_step .== 0)] .= 0.1
  opt.initial_step = initial_step
  (minf, minx, ret) = optimize!(opt, minx)

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
