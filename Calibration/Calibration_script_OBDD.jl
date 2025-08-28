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
tor      = 1.0      # time of origin
ψω_epoch = Float64[tor/2]

# Tree number
i = isempty(ARGS) ? 0 : parse(Int64, ARGS[1])

## Load trees and predict parameters
# Inference parameters
veryShortMCMC = false
shortMCMC = false

λa_prior = (1.5, 1.0)
μa_prior = (1.5, 1.0)
αλ_prior = (0.0, 0.1)
αμ_prior = (0.0, 0.1)
σλ_prior = (5.0, 0.5)
σμ_prior = (5.0, 0.5)
ψ_prior  = (2.0, 1.0)
ω_prior  = (2.0, 0.2)
f_epoch  = Int64[0]
niter    = veryShortMCMC ? 100 : (shortMCMC ? 100_000 : 2_000_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? niter+1 : niter+1)
nburn    = veryShortMCMC ? 0 : (shortMCMC ? 5_000 : 50_000)
nflushθ  = Int64(ceil(niter/20_000))
nflushΞ  = Int64(ceil(niter/100))
tune_int = 100
ϵi       = 0.2
λi       = NaN
μi       = NaN
ψi       = NaN
ωi       = NaN
αλi      = 0.0
αμi      = 0.0
σλi      = 0.1
σμi      = 0.1
pupdp    = (0.02, 0.02, 0.1, 0.01, 0.01, 0.02, 0.1, 0.2)
δt       = 1e-3
survival = true
prints   = 5
stnλ     = 0.5
stnμ     = 0.5
mxthf    = Inf
tρ       = Dict("" => 1.0)

## MCMC
dat_path = "datasets_calibration_OBDD/dataset$i/"
@show dat_path

out_file = "outputs/OBDD_SimulatedOBDD_dataset$(i)_$(niter)iter"
isdir("outputs/") || mkdir("outputs/")

# reconstructed tree
tree_rec = read_newick(dat_path*"tree_reconstructed.tre", true)
ωtimes = filesize(dat_path*"occurrences.csv")>0 ? readdlm(dat_path*"occurrences.csv")[:] : Float64[]
@show tree_rec
@show ωtimes

seed!(i)
insane_gbmobd(tree_rec::sTf_label,
              ωtimes,
              λa_prior = λa_prior,
              μa_prior = μa_prior,
              αλ_prior = αλ_prior,
              αμ_prior = αμ_prior,
              σλ_prior = σλ_prior,
              σμ_prior = σμ_prior,
              ψ_prior  = ψ_prior,
              ω_prior  = ω_prior,
              ψω_epoch = ψω_epoch,
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
              ωi       = ωi,
              αλi      = αλi,
              αμi      = αμi,
              σλi      = σλi,
              σμi      = σμi,
              pupdp    = pupdp,
              δt       = δt,
              survival = survival,
              mxthf    = mxthf,
              prints   = prints,
              stnλ     = stnλ,
              stnμ     = stnμ,
              tρ       = tρ)



