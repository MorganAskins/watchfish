# Watchfish

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/morganaskins/watchfish/master)
[![Build Status](https://travis-ci.com/MorganAskins/watchfish.svg?branch=master)](https://travis-ci.com/MorganAskins/watchfish)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://morganaskins.github.io/watchfish)

Given <img src="svgs/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode" align=middle width=17.73973739999999pt height=22.465723500000017pt/> components, each with an estimated rate <img src="svgs/c90119f20c10a72dd5dccbfc89cd0785.svg?invert_in_darkmode" align=middle width=11.826559799999991pt height=32.16441360000002pt/> determined by a
normal distribution with uncertainty <img src="svgs/3f9f71491501df368c7b0ef70db38d54.svg?invert_in_darkmode" align=middle width=10.747741949999991pt height=23.488575000000026pt/>, calculate the confidence
itervals and perform a hypothesis tests for each parameter <img src="svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?invert_in_darkmode" align=middle width=7.054796099999991pt height=22.831056599999986pt/>.

Nominally each event corresponds to a set of observables <img src="svgs/19e3f7018228f8a8c6559d0ea5500aa2.svg?invert_in_darkmode" align=middle width=10.747741949999991pt height=23.488575000000026pt/> of <img src="svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998994999999pt height=22.465723500000017pt/>
measurements, for any given measurement, the probability for that particular
measurement to come from a particular components is given by

<p align="center"><img src="svgs/21122488bdced611e7ed8bffcd42b543.svg?invert_in_darkmode" align=middle width=38.20684395pt height=16.438356pt/></p>

The prior probability is then formed through a combination of these components
such that the total probability is 

<p align="center"><img src="svgs/33f2760534db4fe806e49280c548fc68.svg?invert_in_darkmode?sanitize=true" align=middle width=99.53078024999999pt height=47.806078649999996pt/></p>

The likelihood for a full data set of <img src="svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998994999999pt height=22.465723500000017pt/> measurements is the product of each
event total probability

<p align="center"><img src="svgs/a364c165d0dd83eefa02e86048508140.svg?invert_in_darkmode" align=middle width=234.3143352pt height=50.399845649999996pt/></p>

We can extend the likelihood by proclaiming that each components as well as the
sum of components are simply a stochastic process, produces the extended
likelihood:

<p align="center"><img src="svgs/1624ea1a01ed7e584701d845dd89b4d7.svg?invert_in_darkmode" align=middle width=249.07579454999998pt height=50.399845649999996pt/></p>

Finally, we can claim that we have _a priori_ knowledge of the parameters,
whether it be through side-band analysis or external constraints, by including
those constraints via some prior probability. Given no specific knowledge of
the shape of that prior, we will consider the information we receive on the
variables to be normally distributed and multiply the likelihood by those
constraints

<p align="center"><img src="svgs/3a59d3e79684a56191cc0718fda97fe6.svg?invert_in_darkmode" align=middle width=444.07050929999997pt height=57.205834949999996pt/></p>

A few definitions to simplify things:
<p align="center"><img src="svgs/5dbbdf06403eb766a9eaf09d34a85148.svg?invert_in_darkmode" align=middle width=74.26263239999999pt height=47.806078649999996pt/></p>

Then then our objective function <img src="svgs/cffa95e5b84679edc10428c782fe8e2f.svg?invert_in_darkmode" align=middle width=78.99089384999999pt height=22.465723500000017pt/>

<p align="center"><img src="svgs/4d56d8302518412d75f39e43357f5a75.svg?invert_in_darkmode" align=middle width=504.59342834999995pt height=50.399845649999996pt/></p>

Finally, for a counting analysis we assume that an optimal set of cuts has been
applied which optimizes the sensitivity to a particular parameter, which
simplifies the likelihood such that
<p align="center"><img src="svgs/a2654de0d1534690d918217913385c8e.svg?invert_in_darkmode" align=middle width=72.90990795pt height=16.438356pt/></p>
Also, because the shape of the likelihood space is independent of constant
parameters, we can drop the <img src="svgs/96a94d4478eee20a423dcae01dfa99ba.svg?invert_in_darkmode" align=middle width=66.15028695pt height=26.045612999999992pt/> terms. We could
also remove the <img src="svgs/be53185ed16aa9438e74a8f287753791.svg?invert_in_darkmode" align=middle width=38.97263864999999pt height=22.831056599999986pt/> term as well, but for numerical stability we will
keep it around, but use Sterling's approximation: <img src="svgs/83b15a890ca32ebd305a780c615eec74.svg?invert_in_darkmode" align=middle width=145.38783435pt height=22.831056599999986pt/>. The remaining objective function we will thus use is:

<p align="center"><img src="svgs/fadac3b98225fbb22296866ea07908f8.svg?invert_in_darkmode" align=middle width=355.9963737pt height=47.806078649999996pt/></p>

_Note: If the different values of <img src="svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode" align=middle width=10.16555099999999pt height=22.831056599999986pt/> differ by orders of magnitude, it
might be worth forming an affine invariant form of the likelihood, otherwise
the <img src="svgs/96a94d4478eee20a423dcae01dfa99ba.svg?invert_in_darkmode" align=middle width=66.15028695pt height=26.045612999999992pt/> term should not matter_
