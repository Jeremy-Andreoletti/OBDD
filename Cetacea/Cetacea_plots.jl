# Set working directory

include(homedir()*"/Nextcloud/Recherche/1_Methods/INSANE/Source_INSANE.jl");

using Distributions
using DataFrames
using DelimitedFiles
using StatsBase
using Random: seed!
using Distributed
using Plots
using Measures

## Read in data
trees_rec = read_newick("Data/Cetacea_Safe_posterior.trees", true)

# Automatically find all output files in the "outputs" directory with a specific pattern
output_files = readdir("outputs"; join=true)  # Get full paths
output_files = filter(x -> occursin(r".txt", x), output_files)
output_files = filter(x -> !occursin(r"\.\_", x), output_files)

# Group files by their common model name, ignoring the seed part
model_groups = Dict{String, Vector{String}}()
for file in output_files
    if occursin(r"_seed", file)
        # Extract the part before "_seed" as the model name
        model_name = match(r"^(.*)000", basename(file)).match
        if haskey(model_groups, model_name)
            push!(model_groups[model_name], file)
        else
            model_groups[model_name] = [file]
        end
    end
end
model_groups

# Iterate over each model, merge seeds, and generate plots/animations
ENV["GKSwstype"] = "nul"
gr(dpi=400, size=(500,300))
# for (model_name, files) in Dict("BDD_Cetaceans_20000000" => model_groups["BDD_Cetaceans_20000000"])
# for (model_name, files) in Dict("COBD_rmFossils_Cetaceans_50000000" => model_groups["COBD_rmFossils_Cetaceans_50000000"])
for (model_name, files) in model_groups
    println("Processing model: $model_name")
    @show occursin(r"DD", model_name)

    files = files[filesize.(files) .!= 0]
    
    if isempty(files)
        print("Empty files\n\n")
        continue
    end

    # Set `ofile` based on the current model name
    ofile = model_name*"iter"

    # MCMC traces
    trace_files = replace.(files, ".txt" => ".log")
    traces = readdlm.(trace_files, '\t', header=true)
    combined_traces = DataFrame(vcat.(traces...)[1], vec(traces[1][2]))

    if occursin(r"^BDD|CBD", model_name)
        trees_rec_mrca = sTf_label.(trees_rec)
        trees_rec_mrca = remove_fossils.(trees_rec_mrca)
        rm_stem!.(trees_rec_mrca)
        tor = maximum(treeheight.(trees_rec_mrca))
    else
        tor = maximum(treeheight.(trees_rec))
    end
    
    # Animations and plots for this specific model
    if occursin(r"DD", model_name)
        δt        = 1e-3
        out_trees = iread.(files, ix=1:100)
        out_trees = out_trees[length.(out_trees) .> 1]

        isempty(out_trees) && continue

        out_trees_seed0 = out_trees[1]
        @show length(out_trees_seed0)
        λmax = ceil(maximum(b.(out_trees_seed0))[1]*2)
        anim_tree = @animate for tree_i in out_trees_seed0
          plot(tree_i, b, shownodes=(false, false, false), clim=(0,λmax))
        end
        mp4(anim_tree, "Animations/$(ofile)_anim_tree_λ.mp4", fps=5)

        tree_type = typeof(out_trees_seed0[lastindex(out_trees_seed0)])
        last_tree_seed0 = tree_type(out_trees_seed0[lastindex(out_trees_seed0)])
        reorder!(last_tree_seed0)
        plot(last_tree_seed0, xlab = "Time before present", tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm, 
             shownodes=(false, false, true), showda=true, shsizes=[1.0, 1.0, 1.0], size=(500,500))
        savefig("Images/$(ofile)_tree.pdf")

        median_tree_seed0 = iquantile(remove_unsampled.(out_trees_seed0), 0.5)
        reorder!(median_tree_seed0)
        plot(median_tree_seed0, b, shownodes=(false, false, false), xlab = "Time before present", 
             tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm, size=(500,500))
        savefig("Images/$(ofile)_seed0_tree_λ.pdf")
        plot(median_tree_seed0, d, shownodes=(false, false, false), xlab = "Time before present", 
             tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm, size=(500,500))
        savefig("Images/$(ofile)_seed0_tree_μ.pdf")
        plot(median_tree_seed0, nd, shownodes=(false, false, false), xlab = "Time before present", 
             tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm, size=(500,500))
        savefig("Images/$(ofile)_seed0_tree_nd.pdf")

        out_trees = reduce(vcat, out_trees)
        @show length(out_trees)
        plot(out_trees, b, δt, fillcolor="#005AB5", linecolor=:darkblue, ylim=(0,1.5))
        plot!(out_trees, d, δt, xlab = "Time before present", ylab="Diversification rates", fillcolor="#DC3220", linecolor=:darkred, 
              tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm)
        savefig("Images/$(ofile)_λTT_μTT.pdf")
    else
        if occursin(r"COBD", model_name)
            out_trees = reduce(vcat, iread.(files, ix=1:100))
        else
            out_trees = reduce(vcat, iread.(files, ix=1:100:10000))
        end
        @show length(out_trees)

        plot(combined_traces.lambda, tor, fillcolor="#005AB5", linecolor=:darkblue, ylim=(0,1.5))
        plot!(combined_traces.mu, tor, fillcolor="#DC3220", linecolor=:darkred, ylim=(0,1.5), 
              xlab = "Time before present", ylab = "Diversification rates", 
              tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm)
        savefig("Images/$(ofile)_λTT_μTT.pdf")
    end

    ψω_epoch = Float64[33.9, 28.1, 23.03, 20.44, 7.246, 5.333]
    if occursin(r"OBD", model_name)
        omega_traces = [combined_traces[!, col] for col in names(combined_traces, r"omega")]
        plot(omega_traces, tor, ψω_epoch, fillcolor=:chocolate1, linecolor=:chocolate3, ylim=(0,1.1), xlab = "Time before present", 
             ylab = "Sampling rate", tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm)
        if occursin(r"rmFossils", model_name)
            savefig("Images/$(ofile)_ωTT.pdf")
        else
            psi_traces = [combined_traces[!, col] for col in names(combined_traces, r"psi")]
            plot!(psi_traces, tor, ψω_epoch, fillcolor=:darkorange4, linecolor=:darkorange4, ylim=(0,1.1), xlab = "Time before present", 
                  ylab = "Sampling rates", tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm)
            savefig("Images/$(ofile)_ψTT_ωTT.pdf")
        end
    elseif occursin(r"FBD", model_name)
        psi_traces = [combined_traces[!, col] for col in names(combined_traces, r"psi")]
        plot(psi_traces, tor, ψω_epoch, fillcolor=:darkorange4, linecolor=:darkorange4, ylim=(0,1.1), xlab = "Time before present", 
             ylab = "Sampling rate", tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm)
        savefig("Images/$(ofile)_ψTT.pdf")
    end
        
    anim_tree = @animate for tree_i in out_trees
        plot(tree_i, shownodes=(true, true, true), showda=true, shsizes=[1.0, 1.0, 1.0])
    end
    mp4(anim_tree, "Animations/$(ofile)_anim_tree.mp4", fps=5)

    plot(ltt.(out_trees), 0.05, ylim=(1,400), scale=:identity, xlab = "Time before present", ylab="Number of lineages", 
         tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm)
    for seed_nb in 0:4
        seed!(seed_nb)
        tree = rand(trees_rec)
        if occursin(r"^BDD|CBD|rmFossils", model_name)
            tree = remove_fossils(tree)
        end
        plot!(ltt(tree), scale=:identity, xlab = "Time before present", ylab="Number of lineages", 
              tickfontsize = 13, guidefontsize = 17, left_margin = 2mm, bottom_margin = 3mm, top_margin = 3mm)
    end
    savefig("Images/$(ofile)_LTT.pdf")
    
    print("\n")
end

## Plot the number of occurrences over time
# Load age uncertainty ranges
ω_ranges = readdlm("Data/Cetacea_occurrences_min_max_age_species_corrected.csv", ',')

# Initialize a container to store sampled times for each seed
sampled_times_per_seed = []

# Loop over seeds and sample uniformly within the uncertainty range
for seed_nb in 0:4
    seed!(seed_nb)
    sampled = map(row -> row[1] + rand() * (row[2] - row[1]), eachrow(ω_ranges))
    push!(sampled_times_per_seed, sampled)
end

# Create a histogram with reversed x-axis and overlay histograms per seed
plot(
    xlabel="Time (Ma)", ylabel="Number of Occurrences", 
    legend=:topleft, xflip=true
)

for (i, sampled_times) in enumerate(sampled_times_per_seed)
    histogram!(
        sampled_times, bins=60, alpha=0.5, label="Seed $(i-1)"
    )
end

# Save the figure
savefig("Images/OccurrencesTT.pdf")

