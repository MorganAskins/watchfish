"""
WatchFish

Julia library designed to build statistical models with associated systematic
uncertainties and provide inference into particular parameters of interest.
Both Bayesian and Frequentist analysis are supported.
"""
module WatchFish

using NLopt
using DataFrames
using DataStructures
using SpecialFunctions
using PyPlot
using Printf

# Primary definitions
include("predefined_functions.jl")
include("variable.jl")
include("likelihood.jl")
include("model.jl")
include("profile.jl")
include("plotter.jl")
include("show.jl")
# Kit definitions
include("modelkits/base.jl")
# Finally, macros
include("macros.jl")

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
