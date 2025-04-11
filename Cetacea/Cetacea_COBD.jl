# Set working directory

include(homedir()*"/Nextcloud/Recherche/1_Methods/INSANE/Source_INSANE.jl");

using Distributions
using DataFrames
using DelimitedFiles
using StatsBase
using Random: seed!
using Distributed
using Plots

# using Profile
# using PProf

## Read in data
trees_rec = read_newick("Data/Cetacea_Safe_posterior.trees", true)
ω_ranges  = readdlm("Data/Cetacea_occurrences_min_max_age_species_corrected.csv", ',')

## Sample a tree in the MCMC trace
seed_nb = isempty(ARGS) ? 0 : parse(Int64, ARGS[1])
seed!(seed_nb)
tree = rand(trees_rec)

## Sample occurrence times uniformly within their age uncertainty range
ωtimes = map(row -> row[1]+rand()*(row[2]-row[1]), eachrow(ω_ranges))

### COBD inference ###

veryShortMCMC = false
shortMCMC = false

λ_prior  = (1.0, 1.0)
μ_prior  = (1.0, 1.0)
ψ_prior  = (1.0, 10.0)
ω_prior  = (1.0, 5.0)
ψω_epoch = Float64[33.9, 28.1, 23.03, 20.44, 7.246, 5.333]
# ψω_epoch = Float64[33.9, 28.1]
f_epoch  = Int64[0]
niter    = veryShortMCMC ? 10 : (shortMCMC ? 10_000 : 2_000_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? niter+1 : niter+1)
nburn    = veryShortMCMC ? 0 : (shortMCMC ? 1_000 : 100_000)
nflushθ  = Int64(ceil(niter/2_000))
nflushΞ  = Int64(ceil(niter/100))
ofile    = "COBD_Cetaceans_$(niter)iter_seed$(seed_nb)" ; isdir("outputs/") || mkdir("outputs/")
ϵi       = 0.4
λi       = NaN
μi       = NaN
ψi       = NaN
ωi       = NaN
# pupdp    = (0.02, 0.02, 0.01, 0.01, 0.1)
pupdp    = (0.02, 0.02, 0.01, 0.01, 0.5)
survival = true
prints   = 5
mxthf    = 0.05
tρ       = Dict("" => 1.0)

# Profile.clear()
# seed!(seed_nb); out = @profile insane_cobd(tree::sTf_label, 
seed!(seed_nb); out = insane_cobd(tree::sTf_label, 
    		                          ωtimes    ,
    		                          λ_prior   = λ_prior,
    		                          μ_prior   = μ_prior,
    		                          ψ_prior   = ψ_prior,
    		                          ω_prior   = ω_prior,
    		                          ψω_epoch  = ψω_epoch,
    		                          f_epoch   = f_epoch,
    		                          niter     = niter,
    		                          nthin     = nthin,
    		                          nburn     = nburn,
                                  nflushθ   = nflushθ,
                                  nflushΞ   = nflushΞ,
    		                          ofile     = "outputs/"*ofile,
    		                          ϵi        = ϵi,
    		                          λi        = λi,
    		                          μi        = μi,
    		                          ψi        = ψi,
    		                          ωi        = ωi,
    		                          pupdp     = pupdp,
    		                          survival  = survival,
    		                          prints    = prints,
    		                          mxthf     = mxthf,
    		                          tρ        = tρ)
# pprof()

out_trees = out[2]
out_trees = iread("outputs/$(ofile).txt")

ENV["GKSwstype"] = "nul"
isdir("Animations/") || mkdir("Animations/") ; 
isdir("Images/") || mkdir("Images/") ; 
gr(dpi=400, size=(500,300))
anim_tree = @animate for tree_i in out_trees
  plot(tree_i, shownodes=(true, true, true), showda=true, shsizes=[1.0, 1.0, 1.0])
end
mp4(anim_tree, "Animations/$(ofile)_anim_tree.mp4", fps=5)

plotω(ltt.(out_trees), ωtimes, 0.05) ; plot!(ltt(tree), scale=:identity)
png("Images/$(ofile)_LTT.png")
