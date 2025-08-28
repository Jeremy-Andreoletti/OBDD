include(homedir()*"/Nextcloud/Recherche/1_Methods/INSANE/Source_INSANE.jl");

using Distributions
using DataFrames
using DelimitedFiles
using StatsBase
using Random: seed!
using Distributed
using Plots

ntrees = 1000  # number of trees

## Simulate datasets to assess model calibration
# Parameters
tor     = 2.0      # time of origin
divrate = false    # wether to reparametrize μ into a diversification rate
const TIMEOUT = 10 # maxumum simulation time, in seconds

# Prior distributions
Gamma11 = Gamma(1.0,1.0)
Gamma12 = Gamma(1.0,2.0)
Gamma15 = Gamma(1.0,5.0)
Exp001 = Exponential(0.01)

# Simulate trees
trees = sTfbd[]
rtrees = sTfbd[]
Ωtimes = Vector{Float64}[]
pars = DataFrame(λ=Float64[], μ=Float64[], ψ1=Float64[], ψ2=Float64[], ω1=Float64[], ω2=Float64[], 
				 tree_length=Float64[], n_extinct=Int64[])
ψω_epoch = Float64[tor/2]

i = 1
s = 0
while i <= ntrees
	# Draw parameters from prior distributions
	seed!(s) ; s+=1
	λ = rand(Gamma11)
	if divrate
		μ = -1.0
		while μ<0
			λminusμ = rand(Exp001)
			μ = λ - λminusμ
		end
	else
		μ = rand(Gamma11)
	end
	ψ1 = rand(Gamma11)
	ψ2 = rand(Gamma11)
	ω1 = rand(Gamma12)
	ω2 = rand(Gamma12)
	
	parms = round.((λ,μ,ψ1,ψ2,ω1,ω2); sigdigits=3)
	@show parms

	tree, ωtimes = sim_cobd(tor, λ, μ, [ψ1, ψ2], [ω1, ω2], ψω_epoch)

	# Condition on the survival
	if ntipsalive(tree)>0
		if ntips(tree)>200
			@warn "More than 200 tips"
		else 
			@show tree
			push!(trees, tree)
			push!(rtrees, reconstructed(deepcopy(tree)))
			push!(pars, (λ, μ, ψ1, ψ2, ω1, ω2, treelength(tree), ntipsextinct(tree)))
			push!(Ωtimes, ωtimes)
			i += 1
		end
	end
end

# Save trees and true parameters
for i in Base.OneTo(ntrees)
	@show i
	path = "datasets_calibration/dataset$i/"
	isdir(path) || mkdir(path)

	write_newick(trees[i], path*"tree_full")
	write_newick(rtrees[i], path*"tree_reconstructed")

	open(path*"occurrences.csv", "w") do io
		writedlm(io, Ωtimes[i])
	end

	open(path*"params.csv", "w") do io
		writedlm(io, [names(pars)])
		writedlm(io, [pars[i,:]])
	end
end


## Load trees and predict parameters
# Inference parameters
veryShortMCMC = false
shortMCMC = false

λ_prior  = (1.0, 0.5)
μ_prior  = (1.0, 1.0)
ψ_prior  = (1.0, 1.0)
ω_prior  = (1.0, 0.2)
f_epoch  = Int64[0]
niter    = veryShortMCMC ? 10 : (shortMCMC ? 10_000 : 500_000)
nthin    = veryShortMCMC ? 1 : (shortMCMC ? 10 : 100)
nburn    = veryShortMCMC ? 0 : (shortMCMC ? 1_000 : 50_000)
nflush   = nthin
ϵi       = 0.4
λi       = NaN
μi       = NaN
ψi       = NaN
ωi       = NaN
pupdp    = (0.02, 0.02, 0.01, 0.01, 0.1)
# pupdp    = (0.02, 0.02, 0.01, 0.01, 0.2)
survival = true
prints   = 5
mxthf    = Inf
tρ       = Dict("" => 1.0)

# Load trees
trees = sTf_label[]; rtrees = sTf_label[] ; Ωtimes = Vector{Float64}[]
pars = DataFrame(λ=Float64[], μ=Float64[], ψ1=Float64[], ψ2=Float64[], ω1=Float64[], ω2=Float64[], 
				 tree_length=Float64[], n_extinct=Int64[])
for i in 1:ntrees
	dat_path = "datasets_calibration_COBD/dataset$i/"
	# complete tree
	tree = read_newick(dat_path*"tree_full.tre", false)
	extinguishpasttips!(tree, accerr)
	
	# reconstructed tree
	tree_rec = read_newick(dat_path*"tree_reconstructed.tre", true)
	
	push!(pars, readdlm(dat_path*"params.csv")[2,:])
	push!(Ωtimes, filesize(dat_path*"occurrences.csv")>0 ? readdlm(dat_path*"occurrences.csv")[:] : Float64[])
	push!(trees, tree)
	push!(rtrees, tree_rec)
end

# MCMC with fossils
for i in 1:ntrees
	dat_path = "datasets_calibration_COBD/dataset$i/"
	@show dat_path

	# With fossils
	seed_nb = i
	out_file = "outputs/COBD_SimulatedCOBD_dataset$(i)_$(niter)iter_seed$(seed_nb)"
	isdir("outputs/") || mkdir("outputs/")
	
	# reconstructed tree
	tree_rec = read_newick(dat_path*"tree_reconstructed.tre", true)
	ωtimes = filesize(dat_path*"occurrences.csv")>0 ? readdlm(dat_path*"occurrences.csv")[:] : Float64[]
	@show tree_rec
	@show ωtimes
	
	seed!(seed_nb)
	out = insane_cobd(tree_rec::sTf_label, 
                      ωtimes    ::Vector{Float64},
                      λ_prior   = λ_prior,
                      μ_prior   = μ_prior,
                      ψ_prior   = ψ_prior,
                      ω_prior   = ω_prior,
                      ψω_epoch  = ψω_epoch,
                      f_epoch   = f_epoch,
                      niter     = niter,
                      nthin     = nthin,
                      nburn     = nburn,
                      nflush    = nflush,
                      ofile     = out_file,
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
	
end







