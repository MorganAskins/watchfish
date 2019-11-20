# Watchfish

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/morganaskins/watchfish/master)
[Presentation.io](https://morganaskins.github.io/watchfish)

Given $M$ components, each with an estimated rate $\vec{\beta}$ determined by a
normal distribution with uncertainty $\vec{\sigma}$, calculate the confidence
itervals and perform a hypothesis tests for each parameter $b$.

Nominally each event corresponds to a set of observables $\vec{x}$ of $N$
measurements, for any given measurement, the probability for that particular
measurement to come from a particular components is given by

$$ P_i(\vec{x}) $$

The prior probability is then formed through a combination of these components
such that the total probability is 

$$ \mathbf{P} = \sum_i^M P_i(\vec{x}) $$

The likelihood for a full data set of $N$ measurements is the product of each
event total probability

$$\mathcal{L}(\vec{x}) = \prod_j^N \left( \sum_i^M b_iP_i(\vec{x}) \right) / \sum_i^Mb_i $$

We can extend the likelihood by proclaiming that each components as well as the
sum of components are simply a stochastic process, produces the extended
likelihood:

$$\mathcal{L}(\vec{x}) = \frac{\text{e}^{-\sum_i^Mb_i}}{N!} \prod_j^N \left( \sum_i^M b_iP_i(\vec{x}) \right) $$

Finally, we can claim that we have _a priori_ knowledge of the parameters,
whether it be through side-band analysis or external constraints, by including
those constraints via some prior probability. Given no specific knowledge of
the shape of that prior, we will consider the information we receive on the
variables to be normally distributed and multiply the likelihood by those
constraints

$$\mathcal{L}(\vec{x}) = \frac{\text{e}^{-\sum_i^Mb_i}}{N!} \prod_j^N \left( \sum_i^M b_iP_i(\vec{x}) \right) \frac{1}{\sqrt{2\pi \sigma_j^2}}\text{exp}\left({\frac{-(\beta_i-b_i)^2}{2\sigma_i^2}}\right)$$

A few definitions to simplify things:
$$ \lambda := \sum_i^Mb_i $$

Then then our objective function $\mathcal{O} = -\text{Ln}\mathcal{L}$

$$\mathcal{O} = \lambda + \text{Ln}N! -\sum_j^N\text{Ln}\left( \sum_i^M b_iP_i(\vec{x}) \right) + \sum_i^M \left( \frac{(\beta_i-b_i)^2}{2\sigma_i^2} + \text{Ln}\sqrt{2\pi \sigma_i} \right)$$

Finally, for a counting analysis we assume that an optimal set of cuts has been
applied which optimizes the sensitivity to a particular parameter, which
simplifies the likelihood such that
$$ P_i(\vec{x}) := 1 $$
Also, because the shape of the likelihood space is independent of constant
parameters, we can drop the $\text{Ln}\sqrt{2\pi \sigma_i}$ terms. We could
also remove the $\text{Ln}N!$ term as well, but for numerical stability we will
keep it around, but use Sterling's approximation: $\text{Ln}N! \approx
N\text{Ln}N - N$. The remaining objective function we will thus use is:

$$\mathcal{O} = \lambda - N\text{Ln}\lambda + N\text{Ln}N - N + \sum_i^M \left( \frac{(\beta_i-b_i)^2}{2\sigma_i^2} \right)$$

_Note: If the different values of $\beta$ differ by orders of magnitude, it
might be worth forming an affine invariant form of the likelihood, otherwise
the $\text{Ln}\sqrt{2\pi \sigma_i}$ term should not matter_
