using Batman, Test

m = CountingExperiment()
add_component!(m, "Signal", 20.0)
add_component!(m, "Bkg 1", 30.0; 0.1)
add_component!(m, "Bkg 2", 30.0; 45.0)
add_component!(m, "Bkg 3", 12.0; 12.0)

set_counts!(m, 100)

results = minimize!(m)

compute_profiled_uncertainties!(results; σ=1)

signal = getparam(m, "Signal")

@test signal.fit ≈ 18 rtol=0.5
