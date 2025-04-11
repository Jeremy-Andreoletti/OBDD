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

## Sample a tree in the MCMC trace
seed_nb = isempty(ARGS) ? 0 : parse(Int64, ARGS[1])
seed!(seed_nb)
tree = rand(trees_rec)

## Remove fossils in the tree
tree = sT_label(remove_fossils(tree))

## Remove the stem branch
rm_stem!(tree);


### BDD inference ###

veryShortMCMC = false
shortMCMC = false

λa_prior = (1.5, 1.0)
μa_prior = (1.5, 1.0)
α_prior  = (0.0, 10.0)
σλ_prior = (5.0, 0.5)
σμ_prior = (5.0, 0.5)
niter    = veryShortMCMC ? 10 : (shortMCMC ? 50_000 : 10_000_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? niter + 1 : niter + 1)
nburn    = veryShortMCMC ? 0 : (shortMCMC ? 1_000 : 50_000)
nflushθ  = Int64(ceil(niter/20_000))
nflushΞ  = Int64(ceil(niter/100))
ofile    = "BDD_Cetaceans_$(niter)iter_seed$(seed_nb)"
ϵi       = 0.2
λi       = NaN
μi       = NaN
αi       = 0.0
σλi      = 0.1
σμi      = 0.1
pupdp    = (0.02, 0.1, 0.01, 0.1, 0.5)
δt       = 1e-3
survival = true
mxthf    = 0.05
prints   = 5
stnλ     = 0.5
stnμ     = 0.5
tρ       = Dict("" => 1.0)


seed!(seed_nb); insane_gbmbd(tree::sT_label,
                             λa_prior = λa_prior,
                             μa_prior = μa_prior,
                             α_prior  = α_prior,
                             σλ_prior = σλ_prior,
                             σμ_prior = σμ_prior,
                             niter    = niter,
                             nthin    = nthin,
                             nburn    = nburn,
                             nflushθ  = nflushθ,
                             nflushΞ  = nflushΞ,
                             ofile    = "outputs/"*ofile,
                             ϵi       = ϵi,
                             λi       = λi,
                             μi       = μi,
                             αi       = αi,
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

# out_trees = iread("outputs/$(ofile).txt")
# #tree_nb = findall(x->x==treelength(remove_extinct(out_trees[1])), treelength.(sT_label.(remove_fossils.(trees_rec))))[1]
# #tree = sT_label(remove_fossils(trees_rec[tree_nb]))

# ENV["GKSwstype"] = "nul"
# gr(dpi=400, size=(500,300))
# anim_tree = @animate for tree_i in out_trees
#   plot(tree_i, shownodes=(true, true, true), showda=true, shsizes=[1.0, 1.0, 1.0])
# end
# mp4(anim_tree, "Animations/$(ofile)_anim_tree.mp4", fps=5)

# gr(dpi=400, size=(500,300))
# anim_tree = @animate for tree_i in out_trees
#   plot(tree_i, b, shownodes=(false, false, false))
# end
# mp4(anim_tree, "Animations/$(ofile)_anim_tree_λ.mp4", fps=5)

# plot(out_trees[lastindex(out_trees)], shownodes=(false, false, true), showda=true, shsizes=[1.0, 1.0, 1.0])
# png("Images/$(ofile)_tree.png")

# median_tree = iquantile(remove_unsampled.(out_trees), 0.5)
# plot(median_tree, b, shownodes=(false, false, false), clim=(0,1.5))
# png("Images/$(ofile)_tree_λ.png")
# plot(median_tree, d, shownodes=(false, false, false), clim=(0,1.5))
# png("Images/$(ofile)_tree_μ.png")
# plot(median_tree, nd, shownodes=(false, false, false), clim=(-1.5,1.5))
# png("Images/$(ofile)_tree_nd.png")

# plot(out_trees, b, δt, ylab="λ(t)")
# png("Images/$(ofile)_λTT.png")

# plot(out_trees, d, δt, ylab="μ(t)")
# png("Images/$(ofile)_μTT.png")

# plot(ltt.(out_trees), 0.05) ; plot!(ltt(tree), scale=:identity)
# png("Images/$(ofile)_LTT.png")



