"""
@addfunction f

Using `@addfunction` will write in a new function into the
current namespace at read-time.

# Examples
```julia-repl
addfunction lognormalconstraint(x, μ, σ) = (x-μ)^2/2/σ
```
"""
macro addfunction(f)
  :($f)
end

## Is this obsolete?
macro dataset( name, df )
  :( $(esc(name)) = $df )
end

"""
@add_dataset(name, df::DataFrame)

Introduce a new dataset into the namespace at read-time.
"""
function add_dataset(name, df::DataFrame)
  eval(:( $name = $df ) )
end
