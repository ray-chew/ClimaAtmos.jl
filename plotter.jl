import ClimaAnalysis
import CairoMakie
using Statistics
import ClimaComms
import ClimaAnalysis.Visualize as viz
import ClimaCore.InputOutput
import Plots # ClimaCorePlots # You may wish to use ClimaCorePlots for viz on unwrapped cubed-sphere
import ClimaCore.Fields

import ClimaAnalysis: slice, window
figres_x = 1000;
figres_y = 400;
days = 86400;

# User to define `basepath` and `dir` which contains `.nc` files 
# reduction = "average" to use default 
# 
basepath = "./"
dir = "output/baroclinic_wave_equil_debug/output_0074/"

simdir = ClimaAnalysis.SimDir(basepath*dir);
reduction = "average";
pfull = get(simdir; short_name = "pfull", reduction);
Δpfull = slice(pfull, z = 50.0) - slice(pfull, z=0.0)


function gen_all_surface_time_plots(var; ts = 0, te = 360, ti = 30, colorrange = (-3,3))
    more_kwargs =
        Dict(
            :plot => Dict(:colorrange => colorrange, :colormap => :balance),
            :axis => Dict(:dim_on_y => true, :limits => (nothing, nothing)))
    
    for ii in var.dims["time"]
        # ClimaAnalysis Plots
        fig = CairoMakie.Figure(size = (figres_x, figres_y))
        viz.plot!(fig, var; more_kwargs, time = ii)
        CairoMakie.save("./plt.png", fig)
    end
end