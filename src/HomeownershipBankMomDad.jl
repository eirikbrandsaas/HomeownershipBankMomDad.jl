
## Load packages
using Plots
using Printf        # For pretty printing
using DataFrames
using Interpolations
using CSV
using Sobol
using LinearAlgebra
using Random
using Statistics
using StatsBase     # Percentiles
using Latexify
using LaTeXStrings
using Dates # Allows you append todays date to filenames
using StableRNGs
using PanelDataTools
using BenchmarkTools

## Load files
include("Julia/struct.jl")
include("Julia/various_functions.jl")
include("Julia/VFI.jl")
include("Julia/decision_problem.jl")
include("Julia/fundamentals.jl")
include("Julia/simpanel.jl")
include("Julia/estimation.jl")
include("Julia/export_estimation.jl")

include("Julia/policyplots.jl")
include("Julia/quant_noalt.jl")
include("Julia/quant_endoprice.jl")
include("Julia/quant_blackwhite.jl")
include("Julia/quant_preferadjustment.jl")
include("Julia/quant_robustness.jl")

include("Julia/postest_plots.jl")
