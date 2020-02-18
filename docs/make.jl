push!(LOAD_PATH, "../src/")
using Documenter, Batman 

makedocs(
  sitename="Batman.jl",
  modules=[Batman],
  pages = Any[
    "Home" => "index.md",
    "Library" => map(s -> "lib/$(s)", sort(readdir(joinpath(@__DIR__, "src/lib"))))
  ]
)

deploydocs(
  repo = "github.com/MorganAskins/Batman.jl.git",
  devbranch = "dev",
  devurl = "dev",
  versions = ["stable"=>"v^", "v#.#", "dev" => "dev", "master"=>"master"]
)
