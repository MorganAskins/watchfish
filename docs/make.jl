push!(LOAD_PATH, "../src/")
using Documenter, WatchFish

makedocs(
  sitename="WatchFish Documentation"
)

deploydocs(
  repo = "github.com/MorganAskins/watchfish.git",
)
