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

### FBDD inference ###

veryShortMCMC = false
shortMCMC = false

λa_prior = (1.5, 1.0)
μa_prior = (1.5, 1.0)
αλ_prior = (0.0, 0.1)
αμ_prior = (0.0, 0.1)
σλ_prior = (5.0, 0.5)
σμ_prior = (5.0, 0.5)
ψ_prior  = (1.0, 5.0)
ψ_epoch  = Float64[33.9, 28.1, 23.03, 20.44, 7.246, 5.333]
f_epoch  = Int64[0]
niter    = veryShortMCMC ? 10 : (shortMCMC ? 50_000 : 10_000_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? niter+1 : niter+1)
nburn    = veryShortMCMC ? 0 : (shortMCMC ? 1_000 : 50_000)
nflushθ  = Int64(ceil(niter/20_000))
nflushΞ  = Int64(ceil(niter/100))
ofile    = "FBDD_Cetaceans_rateαprior0.1_$(niter)iter_seed$(seed_nb)"
ϵi       = 0.2
λi       = NaN
μi       = NaN
ψi       = NaN
αλi      = 0.0
αμi      = 0.0
σλi      = 0.1
σμi      = 0.1
pupdp    = (0.02, 0.02, 0.1, 0.01, 0.01, 0.1, 0.5)
δt       = 1e-3
survival = true
mxthf    = 0.05
prints   = 5
stnλ     = 0.5
stnμ     = 0.5
tρ       = Dict("" => 1.0)


seed!(seed_nb); insane_gbmfbd(tree::sTf_label,
                              λa_prior = λa_prior,
                              μa_prior = μa_prior,
                              αλ_prior = αλ_prior,
                              αμ_prior = αμ_prior,
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
                              ofile    = "outputs/"*ofile,
                              ϵi       = ϵi,
                              λi       = λi,
                              μi       = μi,
                              ψi       = ψi,
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

out_trees = iread("outputs/$(ofile).txt")
#tree = trees_rec[argmin(abs.(treelength.(trees_rec) .- treelength(remove_extinct(out_trees[1]))))]
out_trees = reduce(vcat, [iread("outputs/FBDD_Cetaceans_rateαprior0.1_$(niter)iter_seed$(seed_nb).txt") for seed_nb in 0:4])

ENV["GKSwstype"] = "nul"
gr(dpi=400, size=(500,300))
anim_tree = @animate for tree_i in out_trees
  plot(tree_i, shownodes=(true, true, true), showda=true, shsizes=[1.0, 1.0, 1.0])
end
mp4(anim_tree, "Animations/$(ofile)_anim_tree.mp4", fps=5)

λmax = ceil(maximum(b.(out_trees))[1]*2)
anim_tree = @animate for tree_i in out_trees
  plot(tree_i, b, shownodes=(false, false, false), clim=(0,λmax))
end
mp4(anim_tree, "Animations/$(ofile)_anim_tree_λ.mp4", fps=5)

plot(out_trees[lastindex(out_trees)], shownodes=(false, false, true), showda=true, shsizes=[1.0, 1.0, 1.0])
png("Images/$(ofile)_tree.png")

median_tree = iquantile(remove_unsampled.(out_trees), 0.5)
plot(median_tree, b, shownodes=(false, false, false), clim=(0,1.5))
png("Images/$(ofile)_tree_λ.png")
plot(median_tree, d, shownodes=(false, false, false), clim=(0,1.5))
png("Images/$(ofile)_tree_μ.png")
plot(median_tree, nd, shownodes=(false, false, false), clim=(-1.5,1.5))
png("Images/$(ofile)_tree_nd.png")

plot(out_trees, b, δt, ylab="λ(t)")
png("Images/$(ofile)_λTT.png")

plot(out_trees, d, δt, ylab="μ(t)")
png("Images/$(ofile)_μTT.png")

plot(out_trees, b, δt, ylab="λ(t) and μ(t)")
plot!(out_trees, d, δt)
png("Images/$(ofile)_λTT_μTT.png")

plot(ltt.(out_trees), 0.05) ; plot!(ltt(tree), scale=:identity)
png("Images/$(ofile)_LTT.png")



