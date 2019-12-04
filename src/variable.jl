abstract type Variable end

mutable struct Parameter <: Variable
  name::String
  info::String
  min::Float64
  max::Float64
  init::Float64
  units::String
  constant::Bool
  bounds::Tuple{Float64, Float64}
  # Post fit information
  fit::Float64
  low::Float64
  high::Float64
  likelihood_x::Array{Float64}
  likelihood_y::Array{Float64}
  function Parameter( name::String; info::String="", min::Float64=-Inf,
                      max::Float64=Inf, init::Float64=0.0, units::String="", 
                      constant=false, bounds::Tuple{Float64,Float64}=(-Inf, Inf) )
    new( name, info, min, max, init, units, constant, bounds,
         0.0, 0.0, 0.0,
         Array{Float64}(undef, 0),
         Array{Float64}(undef, 0),
       )
  end
end

mutable struct Constant <: Variable
  name::String
  info::String
  init::Float64
  units::String
  constant::Bool
  function Constant( name::String, init::Float64; 
                     info::String="", units::String="" ) 
    new( name, info, init, units, true )
  end
end

struct Dataset <: Variable
  name::String
  init::String
  constant::Bool
  #init::Symbol
  function Dataset( name::String ; info::String="" )
    new( name, name, true )
  end
end
