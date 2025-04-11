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

## Sample a tree in the MCMC trace
seed_nb = isempty(ARGS) ? 0 : parse(Int64, ARGS[1])
seed!(seed_nb)
tree = rand(trees_rec)

## Remove fossils in the tree
tree = sT_label(remove_fossils(tree))

## Remove the stem branch
rm_stem!(tree);

### CBD inference ###
veryShortMCMC = false
shortMCMC = false

λ_prior  = (1.0, 1.0)
μ_prior  = (1.0, 1.0)
niter    = veryShortMCMC ? 10 : (shortMCMC ? 100_000 : 2_000_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? 100 : 1_000)
nburn    = veryShortMCMC ? 0 : (shortMCMC ? 10_000 : 50_000)
nflush   = nthin
ofile    = "CBD_Cetaceans_$(niter)iter_seed$(seed_nb)" ; isdir("outputs/") || mkdir("outputs/")
ϵi       = 0.4
λi       = NaN
μi       = NaN
pupdp    = (0.2,0.2,0.2)
survival = true
mxthf    = 0.05
tρ       = Dict("" => 1.0)

# Profile.clear()
# seed!(0); out = @profile insane_cbd(tree::sTf_label, 
seed!(seed_nb); out = insane_cbd(tree::sT_label, 
    		                         λ_prior   = λ_prior,
    		                         μ_prior   = μ_prior,
    		                         niter     = niter,
    		                         nthin     = nthin,
    		                         nburn     = nburn,
    		                         nflush    = nflush,
    		                         ofile     = "outputs/"*ofile,
    		                         ϵi        = ϵi,
    		                         λi        = λi,
    		                         μi        = μi,
    		                         pupdp     = pupdp,
    		                         survival  = survival,
    		                         mxthf     = mxthf,
    		                         tρ        = tρ)
# pprof()

#out_trees = out[2]
# out_trees = iread("outputs/$(ofile).txt")

# ENV["GKSwstype"] = "nul"
# isdir("Animations/") || mkdir("Animations/") ; 
# isdir("Images/") || mkdir("Images/") ; 
# gr(dpi=400, size=(500,300))
# anim_tree = @animate for tree_i in out_trees
#   plot(tree_i, shownodes=(true, true, true), showda=true, shsizes=[1.0, 1.0, 1.0])
# end
# mp4(anim_tree, "Animations/$(ofile)_anim_tree.mp4", fps=5)

# plot(ltt.(out_trees), 0.05) ; plot!(ltt(tree), scale=:identity)
# png("Images/$(ofile)_LTT.png")


