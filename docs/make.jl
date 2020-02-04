push!(LOAD_PATH, "../src/")
using Documenter, Batman 

makedocs(
  sitename="Batman.jl",
  modules=[Batman],
  pages = Any[
    "Home" => "index.md"
  ]
)

deploydocs(
  repo = "github.com/MorganAskins/Batman.jl.git",
)
