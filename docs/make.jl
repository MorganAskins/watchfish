push!(LOAD_PATH, "../src/")
using Documenter, Batman 

makedocs(
  sitename="BATMAN Documentation"
)

deploydocs(
  repo = "github.com/MorganAskins/Batman.jl.git",
)
