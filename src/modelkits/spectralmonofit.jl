import Distributions: Poisson

struct SingleSpectra
  f
  p::Parameter
  o::Observable
  df::DataFrame
end

"""
    SpectralMonofit()

Multidimensional fit to one-dimensional distributions (neglecting
correlations between observables).

# Example
```julia-repl
m = SpectralMonofit()
```
"""
mutable struct SpectralMonofit <: Model
  params
  pdflist::Array{NLogPDF}
  nll::NLogLikelihood
  # Additional parameters
  spectra::Array{SingleSpectra}
  datasets::Array{DataFrame}
  observables::Array{Observable}
  function SpectralMonofit(observables::Symbol...)
    new( [], [], NLogLikelihood(), 
        Array{SingleSpectra}(undef, 0), Array{DataFrame}(undef, 0), Array{Observable}(undef, 0) )
  end
end

struct SpectralPDF
  p::Parameter
  f::String
end

function generate_mock_dataset(m::SpectralMonofit; kwargs...)
  kwargs = Dict(kwargs)
  # Clean out the old data
  for d in m.datasets
    deleterows!(d, 1:size(d, 1))
  end
  for p in m.params
    p.fit = rand(Poisson(p.init))
  end
  for sp in m.spectra
    append!(sp.df[!, sp.o.sname], 
            randomFromSpectrum(sp.f; count=sp.p.fit, min=sp.o.min, max=sp.o.max))
  end
  # Should global the data
  for d in m.datasets
    add_dataset(Symbol(d), d)
  end
end

function constructPDF!(m::SpectralMonofit, p::Parameter, shapes, observables, ds::Symbol)
  df = eval(ds)
  if !(df in m.datasets)
    push!(m.datasets, df)
  end
  seval = ""
  for (shp, obs) in zip(shapes, observables)
    fSymbol = Symbol(p.name*"_"*String(obs)*"_"* String(ds))
    add_function(fSymbol, shp)
    #push!(m.spectra, SingleSpectra(p, shp, obs, df))
    push!(m.spectra, SingleSpectra(shp, p, observableFromSymbol(m,obs), df))
    seval *= "$fSymbol($ds[!, :$obs]).*"
  end
  seval = seval[1:end-2]
  @debug seval
  unique_name = "spectralPDF_"*p.name*"_"*String(ds)
  add_function(unique_name, eval(Meta.parse("()->"*seval)))
  SpectralPDF(p, unique_name)
end

function combinePDFs!(m::SpectralMonofit, spectralpdfs::Array{SpectralPDF}, 
                      ds::Symbol; kwargs... )
  kwargs = Dict(kwargs)
  livetime = get(kwargs, :livetime, nothing)
  eff = get(kwargs, :efficiency, Array{Nothing}(undef, size(spectralpdfs)))

  seval = "-sum(batlog.("
  fcall = "("
  extended = "+sum("
  plist = []
  flist = []
  if livetime != nothing
    push!(plist, livetime)
    tname = livetime.name
    fcall *= "$tname,"
  end
  for (spdf,e) in zip(spectralpdfs, eff)
    p = spdf.p
    f = spdf.f
    push!(plist, p)
    push!(flist, f)
    vname = p.name
    if e != nothing
      ename = e.name
      seval *= "$ename*"
      fcall *= "$ename,"
      extended *= "$ename*"
      push!(plist, e)
    end
    if livetime != nothing
      tname = livetime.name
      seval *= "$tname*"
      extended *= "$tname*"
    end
    seval *= "$vname*$f().+"
    fcall *= "$vname,"
    extended *= "$vname+"
  end
  seval=seval[1:end-2]
  seval*="))"
  fcall=fcall[1:end-1]
  fcall*=")->"
  extended=extended[1:end-1]
  n = size(eval(ds), 1)
  #extended*=")"*"+$n*log($n)-$n"
  extended*=")"*"+size($ds,1)*log(size($ds,1))-size($ds,1)"
  unique_name = "fullSpectralPDF_"*String(ds)
  @debug fcall*seval*extended
  add_function(unique_name, eval(Meta.parse(fcall*seval*extended)))
  push!(m.pdflist, NLogPDF(unique_name, plist...))
end

function add_parameter!(m::SpectralMonofit, name::String, init; σ=Inf )
  p = getparam(m, name)
  if p != nothing
    @warn name*" already exists in model."
    return p
  end
  p = Parameter(name; init=init)
  push!(m.params, p)
  if σ != Inf
    a = Constant(name*"_mu", init)
    b = Constant(name*"_sig", σ)
    nlp = NLogPDF("lognormal", p, a, b)
    push!(m.pdflist, nlp)
  end
  return p
end

function add_constant!(m::SpectralMonofit, name::String, init)
  k = getparam(m, name)
  if k != nothing
    @warn name*" already exists in model."
    return k
  end
  Constant(name, init)
end

function add_observable!(m::SpectralMonofit, name::Symbol, min::Float64, max::Float64)
  push!(m.observables, Observable(String(name), name, min, max))
end

function observableFromSymbol(m::SpectralMonofit, name::Symbol)
  for obs in m.observables
    if obs.sname == name
      return obs
    end
  end
  return nothing
end

function minimize!( m::SpectralMonofit; options=Dict() )
  likelihood = NLogLikelihood( m.pdflist )
  add_likelihood!( m, likelihood )
  optimize_model!( m, likelihood; options=options )
end

"""
    randomFromSpectrum(f; min=0.0, max=1.0, top=1.0, count=1)

Return an `Array{Float64}` set of random elements from a 1D pdf given by `f`
between the values of `min` and `max` with a maximum probability of `top`.
"""
function randomFromSpectrum(f; min=0.0, max=1.0, top=1.0, count=1)
  xArr = Array{Float64}(undef, 0)
  for c in collect(1:count)
    while true
      x = rand()*(max - min) + min
      y = rand()*top
      if y < f(x)
        push!(xArr, x)
        break
      end
    end
  end
  return xArr
end
