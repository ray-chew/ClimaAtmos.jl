using Revise, Infiltrator
import ClimaCore: Fields, Spaces, Operators, Utilities, Grids
import ClimaComms
using CUDA
import ClimaComms, ClimaComms.@import_required_backends
import ClimaAtmos as CA
import SciMLBase: step!

# ENV["CLIMACOMMS_DEVICE"]="CUDA" # If undefined defined at startup on `clima`
config_file = "./config/model_configs/baroclinic_wave_equil_debug.yml"
simulation = CA.AtmosSimulation(config_file)
CA.solve_atmos!(simulation)