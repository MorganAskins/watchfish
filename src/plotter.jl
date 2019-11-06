# Plotting functions
using PyPlot

function interval_plot(results, name)
  comp = results.model.component_dict[name]
  # Check if these are already computed, if not compute
  if length(comp.likelihood_x) == 0
    compute_profiled_uncertainties!(results, 0.68)
  end
  x, y = comp.likelihood_x, comp.likelihood_y
  low, high = comp.low, comp.high
  mid = comp.fit
  ns, ps = abs(mid - low), abs(mid - high)
  labeluse = "$name: " *
    " \$\\mu = $(floor(mid*100)/100)" *
    "^{+$(floor(ps*100)/100)}" *
    "_{-$(ns)}\$"
  plot(x, y, label=labeluse)
  vlines(low, 0, y[x .>= low][1], color="red")
  @show y[x .>= low][1]
  vlines(high, 0, y[x .>= high][1], color="red")
  ylim(0, 1.2)
  xlim(minimum(x), maximum(x))
  legend(loc=2, frameon=false)
end
