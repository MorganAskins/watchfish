macro addfunction(f)
  :($f)
end

macro dataset( name, df )
  :( $(esc(name)) = $df )
end

## Do functions which eval expressions count?
function add_dataset(name, df::DataFrame)
  eval(:( $name = $df ) )
end
