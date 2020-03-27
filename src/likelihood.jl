import StatsBase
import Interpolations; itp=Interpolations
import KernelDensity

"""
    NLogPDF(f::String, p::Variable...)

Negative log PDF used when contructing the complete likelihood
function. The NLogPDF should have a unique name `f::String`.
"""
mutable struct NLogPDF 
  func::String
  params::Array{Variable}
  function NLogPDF(f::String, p::Variable...)
    new(f, [p...])
  end
end

"""
    NLogLikelihood(pf::Array{NLogPDF})

Given an array of `NLogPDF`, a complete negative log-likelihood
function is contructed which can then be passed to a minimizer.
"""
mutable struct NLogLikelihood
  objective::Function
  numparams::Int64
  variableList
  parameters
  lower_bounds
  upper_bounds
  function NLogLikelihood(pf::Array{NLogPDF})
    var_list = []
    for p in pf
      for v in p.params
        push!(var_list, v)
      end
    end
    vv = unique(var_list)
    params = [v for v in vv if v.constant == false]
    lower_bounds = [v.bounds[1] for v in vv if v.constant == false]
    upper_bounds = [v.bounds[2] for v in vv if v.constant == false]
    nparams = length(params)
    npdf = length(pf)
    matches = []
    for (idx, p) in enumerate(pf)
      match = []
      for v in p.params
        loc = findlast(x->x==v, params)
        push!(match, loc)
      end
      push!(matches, match)
    end
    seval = ""
    for (idx, p) in enumerate(pf)
      a = matches[idx]
      name = pf[idx].func
      seval *= "$name("
      for (i,v) in zip(a, p.params)
        if v.constant
          val_use = v.init
          seval *= "$val_use,"
        else
          seval *= "x[$i],"
        end
      end
      seval *= ") +"
    end
    seval = seval[1:end-1]
    @debug "Likelihood function" seval
    ff = eval(Meta.parse("x->"*seval))
    new(ff, nparams, vv, params, lower_bounds, upper_bounds)
  end
  function NLogLikelihood()
    new(x->x, 0, [], [])
  end
end

"""
    HistogramPDF(df::DataFrame, axis::Symbol; kwargs...)

Given a set of data `df`, along an axis `axis`, build a histogram
to be used as a PDF in a fit.

# Arguments
- `df::DataFrame`: Data to build the PDF
- `axis::Symbol`: Axis along which to slice the data
- `ndim::Integer=1`: Number of dimensions
- `bins=100`: Number of bins, or collection of bins.
- `interpolate=Interpolations.Constant()`: Constant or Linear interpolation
- `extrapolate=0`: Flat(), Linear(), or 0 extrapolation.
"""
function HistogramPDF(df::DataFrame, axis::Symbol; kwargs...)
  kwargs = Dict(kwargs)
  ndim = get(kwargs, :dimensions, 1)
  bins = get(kwargs, :bins, 100)
  itp_model = get(kwargs, :interpolate, itp.Constant())
  etp_model = get(kwargs, :extrapolate, 0)
  h = StatsBase.fit( StatsBase.Histogram, df[!, axis], bins)
  hx, hy = (h.edges[1][2:end] + h.edges[1][1:end-1])/2.0, h.weights
  #hx, hy = h.edges[1][2:end], h.weights
  #hy = hy ./ sum(hy) ./ (hx[2]-hx[1])
  hy = hy ./ sum(hy)
  itp.extrapolate( 
      itp.interpolate((hx,), hy, itp.Gridded(itp_model)), etp_model
  )
end

"""
    KdePDF(df::DataFrame, axis::Symbol; kwargs...)

Given a set of data `df`, along an axis `axis`, produce a kernel
density estimate to use as a PDF in a fit.

# Arguments
- `df::DataFrame`: Data to build the PDF
- `axis::Symbol`: Axis along which to slice the data
- `ndim::Integer=1`: Number of dimensions
"""
function KdePDF(df::DataFrame, axis::Symbol; kwargs...)
  kwargs = Dict(kwargs)
  # :dimensions not used currently
  ndim = get(kwargs, :dimensions, 1)

  k = KernelDensity.kde(df[!, axis])
  imd = KernelDensity.InterpKDE(k)
  p(x) = KernelDensity.pdf(imd, x)
  return p
end
