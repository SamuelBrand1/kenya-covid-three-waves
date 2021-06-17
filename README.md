# kenya-covid-three-waves
 Simulation and inference code for SARS-CoV-2 transmission within and between socio-economic groups in Kenya over the first three waves of COVID-19 in Kenya.

This code repository contains the underlying model code for the paper *COVID-19 Transmission Dynamics Underlying Epidemic Waves in Kenya*.
## Prerequisites and recommended background knowledge:
* Access and basic familiarity with the [Julia programming language](https://julialang.org/).
* The underlying dynamical system representing the unobserved infection process is a modification to the basic SEIR model that is described in [this paper](https://journals.sagepub.com/doi/full/10.1177/0962280217747054), albeit we implement a continuous time rather than discrete version of the model.
* Solutions of the infection process are generated using the performant, and well documented, package [SciML/DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). Familiarisation with this package is desirable.
* Hamiltonian MCMC (HMC) is implemented using the [dynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl) package. The log-likelihood function for parameters is directly defined in KenyaSerology, and log-likelihood gradients (necessary for HMC) are calculated using forward-mode automatic differentiation. The combination of ODE solutions and log-likelihood function gradients in code was inspired by [DiffEqBayes.jl](https://github.com/SciML/DiffEqBayes.jl). A good conceptual introduction to HMC can be found [here](https://arxiv.org/abs/1701.02434).
* MCMC posterior draws for parameters are stored as a `Chains` struct from the `MCMCChains` [package](https://github.com/TuringLang/MCMCChains.jl). The MCMC draws for each county are available from the `JLD2` objects stored in `\modelfits`, e.g.
```julia
using JLD2,MCMCChains
import KenyaCoVSD
@load("modelfits/Nairobi_model.jld2") #Loads an object called model into Main scope which contains information about Nairobi county
model.MCMC_results.chain #Summary information of MCMC posterior draws
```

## Data sets in this repository

`/opendatacsvs`folder contains .csv datafiles (see *supplementary information* for the main manuscript to read data file captions). This data is also present on the repository in the form of .jld2 datafiles, which are directly used in the tutorial notebooks (see below).

## Supplementary Information

Please find the supplementary information containing further details on the data analysis, and mathematical/statistical reasoning and assumptions behind the model in `/supplementary_info` as a Word file.

## Getting started

The source code for the module `KenyaCoVSD` is located in `\src`. `KenyaCoVSD` is not part of the Julia package ecosystem general directory, therefore, to use this module:

1. Clone this project using e.g. git clone.
2. Activate the Julia REPL.
3. cd to the directory containing the `Project.toml` and `Manifest.toml` files for the cloned repository.
4. Open the Package manager using `]`, activate and instantiate your copy of `KenyaCoVSD` by

```julia
(v1.6) pkg> activate .
(KenyaCoVSD) pkg> instantiate
```
This will download and precompile the dependencies for `KenyaCoVSD`.
