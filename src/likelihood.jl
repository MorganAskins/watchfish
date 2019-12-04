mutable struct NLogPDF 
  func::String
  params::Array{Variable}
  function NLogPDF(f::String, p::Variable...)
    new(f, [p...])
  end
end

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
    #@show seval
    ff = eval(Meta.parse("x->"*seval))
    new(ff, nparams, vv, params, lower_bounds, upper_bounds)
  end
  function NLogLikelihood()
    new(x->x, 0, [], [])
  end
end
