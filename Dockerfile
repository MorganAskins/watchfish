FROM julia:1.3
LABEL maintainer="Morgan Askins <maskins@berkeley.edu>"

RUN julia -e 'using Pkg;\
              Pkg.REPLMode.pkgstr("add https://github.com/MorganAskins/Batman.jl.git ;precompile");\
              using Batman'
