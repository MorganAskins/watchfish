using Random

mutable struct MCMCSample
  pdf
end

function _transition_model(x)
  randn()*0.1 + x
end

function _acceptance(v, vnew, x, y)
  if (x<-1) || (x>1) || (y<-1) || (y>1)
    return false
  end
  if vnew < v
    return true
  else
    accept = rand()
    return (accept < exp(v - vnew))
  end
end

function run_mcmc!( mcmc::MCMCSample, seed ; iterations=1_000, walkers=1)
  accepted = []
  rejected = []
