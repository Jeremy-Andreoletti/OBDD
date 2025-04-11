# Set working directory
include(homedir()*"/Nextcloud/Recherche/1_Methods/INSANE/Source_INSANE.jl");

using Distributions
using DataFrames
using DelimitedFiles
using StatsBase
using Random: seed!
using Distributed
using Plots

## Simulate datasets to assess model calibration
# Parameters
tor     = 1.0      # time of origin
ψ_epoch = Float64[tor/2]

# Tree number
i = isempty(ARGS) ? 0 : parse(Int64, ARGS[1])

## Load trees and predict parameters
# Inference parameters
veryShortMCMC = false
shortMCMC = false

λa_prior = (1.5, 1.0)
μa_prior = (1.5, 1.0)
α_prior  = (0.0, 0.5)
σλ_prior = (5.0, 0.5)
σμ_prior = (5.0, 0.5)
ψ_prior  = (2.0, 1.0)
f_epoch  = Int64[0]
niter    = veryShortMCMC ? 100 : (shortMCMC ? 50_000 : 2_000_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? 50 : 400)
nburn    = veryShortMCMC ? 0 : (shortMCMC ? 5_000 : 100_000)
nflushθ  = nthin
nflushΞ  = Int64(ceil(niter/100))
tune_int = 100 
ϵi       = 0.2
λi       = NaN
μi       = NaN
ψi       = NaN
λtni     = 1.0
μtni     = 1.0
obj_ar   = 0.4
αi       = 0.0
σλi      = 0.1
σμi      = 0.1
λtni     = 0.1
obj_ar   = 0.234
pupdp    = (0.0, 0.05, 0.01, 0.0, 0.1, 0.2)
δt       = 1e-3
survival = true
prints   = 5
mxthf    = Inf
tρ       = Dict("" => 1.0)

# MCMC with fossils
dat_path = "datasets_calibration_fixedAlpha/dataset$i/"
@show dat_path

# With fossils
out_file = "outputs/FBDD_SimulatedOBDD_fixedAlpha_dataset$(i)_$(niter)iter"
isdir("outputs/") || mkdir("outputs/")

# reconstructed tree
tree_rec = read_newick(dat_path*"tree_reconstructed.tre", true)
@show tree_rec

seed!(i); insane_gbmfbd(tree_rec::sTf_label,
              λa_prior = λa_prior,
              μa_prior = μa_prior,
              α_prior  = α_prior,
              σλ_prior = σλ_prior,
              σμ_prior = σμ_prior,
              ψ_prior  = ψ_prior,
              ψ_epoch  = ψ_epoch,
              f_epoch  = f_epoch,
              niter    = niter,
              nthin    = nthin,
              nburn    = nburn,
              nflushθ  = nflushθ,
              nflushΞ  = nflushΞ,
              ofile    = out_file,
              tune_int = tune_int,
              ϵi       = ϵi,
              λi       = λi,
              μi       = μi,
              ψi       = ψi,
              αi       = αi,
              σλi      = σλi,
              σμi      = σμi,
              λtni     = λtni,
              obj_ar   = obj_ar,
              pupdp    = pupdp,
              δt       = δt,
              survival = survival,
              mxthf    = mxthf,
              prints   = prints,
              tρ       = tρ)





