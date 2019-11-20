function getparam(m::Model, name::String)
  for (i, p) in enumerate(m.params)
    if p.name == name
      return p
    end
  end
  return nothing
end
