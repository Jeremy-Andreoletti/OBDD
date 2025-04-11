# Set working directory

include(homedir()*"/Nextcloud/Recherche/1_Methods/INSANE/Source_INSANE.jl");

using Distributions
using DataFrames
using DelimitedFiles
using StatsBase
using Random: seed!
using Distributed
using Plots

#using Profile
#using PProf

## Read in data
trees_rec = read_newick("Data/Cetacea_Safe_posterior.trees", true)

## Sample a tree in the MCMC trace
seed_nb = isempty(ARGS) ? 0 : parse(Int64, ARGS[1])
seed!(seed_nb)
tree = rand(trees_rec)

### CFBD inference ###

veryShortMCMC = false
shortMCMC = false

λ_prior  = (1.0, 1.0)
μ_prior  = (1.0, 1.0)
ψ_prior  = (1.0, 10.0)
ψ_epoch  = Float64[33.9, 28.1, 23.03, 20.44, 7.246, 5.333]
# ψ_epoch = Float64[33.9, 28.1]
f_epoch  = Int64[0]
niter    = veryShortMCMC ? 10 : (shortMCMC ? 10_000 : 2_000_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? 100 : 1_000)
nburn    = veryShortMCMC ? 0 : (shortMCMC ? 1_000 : 100_000)
nflush   = nthin
ofile    = "CFBD_Cetaceans_$(niter)iter_seed$(seed_nb)" ; isdir("outputs/") || mkdir("outputs/")
ϵi       = 0.4
λi       = NaN
μi       = NaN
ψi       = NaN
# pupdp    = (0.02, 0.02, 0.01, 0.1)
pupdp    = (0.02, 0.02, 0.01, 0.5)
survival = true
prints   = 5
mxthf    = 0.05
tρ       = Dict("" => 1.0)

# Profile.clear()
# seed!(seed_nb); out = @profile insane_cfbd(tree::sTf_label, 
seed!(seed_nb); out = insane_cfbd(tree::sTf_label, 
    		                          λ_prior   = λ_prior,
    		                          μ_prior   = μ_prior,
    		                          ψ_prior   = ψ_prior,
    		                          ψ_epoch   = ψ_epoch,
    		                          f_epoch   = f_epoch,
    		                          niter     = niter,
    		                          nthin     = nthin,
    		                          nburn     = nburn,
    		                          nflush    = nflush,
    		                          ofile     = "outputs/"*ofile,
    		                          ϵi        = ϵi,
    		                          λi        = λi,
    		                          μi        = μi,
    		                          ψi        = ψi,
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
anim_tree = @animate for tree_i in out_trees[ Integer.(round.(collect(range(1,length(out_trees), length=100))))]
  plot(tree_i, shownodes=(true, true, true), showda=true, shsizes=[1.0, 1.0, 1.0])
end
mp4(anim_tree, "Animations/$(ofile)_anim_tree.mp4", fps=5)

plot(ltt.(last(out_trees, 500)), 0.05) ; plot!(ltt(tree), scale=:identity)
png("Images/$(ofile)_LTT.png")
