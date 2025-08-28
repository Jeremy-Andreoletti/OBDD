# Overview
This project introduces the **Occurrence Birth-Death Diffusion (OBDD)** model and applies it to study the diversification dynamics of cetaceans.  
It builds on the [`Tapestree.jl`](https://github.com/ignacioq/Tapestree.jl) framework, with additional functionality currently available only in the `insane` branch of my fork.

# Installation

Until these changes are merged upstream, you need to install Tapestree directly from the fork:

```julia
using Pkg
Pkg.add(url="https://github.com/Jeremy-Andreoletti/Tapestree.jl", rev="insane")
```
Note: If you already had Tapestree installed, first run:

```julia
Pkg.rm("Tapestree")
```

# Directory Structure and Description
- `Calibration`: Validating the OBDD model.
- `Cetacea`: Diversification analyses.

# Contact Information
For any questions or replication issues, please reach out to:

Jérémy Andréoletti: jeremy.andreoletti@gmail.com