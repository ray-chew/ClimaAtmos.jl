redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))
import ClimaAtmos as CA
import ClimaCore
import ClimaCore: DataLayouts, Fields, Geometry
import ClimaCore.Fields: Field, FieldVector, field_values
import ClimaCore.DataLayouts: AbstractData
import ClimaCore.Geometry: AxisTensor
import ClimaCore.Spaces: AbstractSpace
import ClimaComms
pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import Logging
using Test

import Random
Random.seed!(1234)

const device = ClimaComms.device()
const comms_ctx = ClimaComms.context(device)
ClimaComms.init(comms_ctx)

# This test checks that:

# 1. A simulation, saved to a checkpoint, is read back identically (up to some
#   tolerance and excluding those fields that are computed during the
#   calculation of the tendencies)
# 2. A simulation, saved to a previous checkpoint, and read back and evolved to
#   the same time is identical (up to some tolerance)
# 3. ClimaAtmos can automatically detect restarts

# Technical note:
#
# Test.jl really wants to give you a stacktrace for failing tests. This seems to be
# hardcoded in the package and not easy to change without defining a whole new
# AbstractTestSet. We don't want stacktraces, we just want to know which fields are
# different.
#
# For this reason, we don't use Test but just print to screen the differences.
# However, we still have to return an exit code with failure in case of the
# comparison fails. So, we have this global `SUCCESS` bool that is updated by
# the result of tests.
const SUCCESS::Base.RefValue{Bool} = Ref(true)

"""
    _error(arr1::AbstractArray, arr2::AbstractArray; ABS_TOL = 100eps(eltype(arr1)))

We compute the error in this way:
- when the absolute value is larger than ABS_TOL, we use the absolute error
- in the other cases, we compare the relative errors
"""
function _error(
    arr1::AbstractArray,
    arr2::AbstractArray;
    ABS_TOL = 100eps(eltype(arr1)),
)
    diff = abs.(arr1 .- arr2)
    denominator = abs.(arr1)
    error = ifelse.(denominator .> ABS_TOL, diff ./ denominator, diff)
    return error
end


"""
    compare(v1, v2; name = "", ignore = Set([:rc]))

Return whether `v1` and `v2` are the same (up to floating point errors).

`compare` walks through all the properties in `v1` and `v2` until it finds
that there are no more properties. At that point, `compare` tries to match the
resulting objects. When such objects are arrays with floating point, `compare`
defines a notion of `error` that is the following: when the absolute value is
less than `100eps(eltype)`, `error = absolute_error`, otherwise it is relative
error. The `error` is then compared against a tolerance.

Keyword arguments
=================

- `name` is used to collect the name of the property while we go recursively
  over all the properties. You can pass a base name.
- `ignore` is a collection of `Symbol`s that identify properties that are
  ignored when walking through the tree. This is useful for properties that
  are known to be different (e.g., `output_dir`).

`:rc` is some CUDA/CuArray internal object that we don't care about
"""
function compare(
    v1::T,
    v2::T;
    name = "",
    ignore = Set([:rc]),
) where {T <: Union{FieldVector, CA.AtmosCache, AbstractSpace}}
    pass = true
    return _compare(pass, v1, v2; name, ignore)
end

function _compare(pass, v1::T, v2::T; name, ignore) where {T}
    properties = filter(x -> !(x in ignore), propertynames(v1))
    if isempty(properties)
        pass &= _compare(v1, v2; name, ignore)
    else
        # Recursive case
        for p in properties
            pass &= _compare(
                pass,
                getproperty(v1, p),
                getproperty(v2, p);
                name = "$(name).$(p)",
                ignore,
            )
        end
    end
    return pass
end

function _compare(v1::T, v2::T; name, ignore) where {T}
    return print_maybe(v1 == v2, "$name differs")
end

function _compare(
    v1::T,
    v2::T;
    name,
    ignore,
) where {T <: Union{AbstractString, Symbol}}
    # What we can safely print without filling STDOUT
    return print_maybe(v1 == v2, "$name differs: $v1 vs $v2")
end

function _compare(v1::T, v2::T; name, ignore) where {T <: Number}
    # We check with triple equal so that we also catch NaNs being equal
    return print_maybe(v1 === v2, "$name differs: $v1 vs $v2")
end

function _compare(
    v1::T,
    v2::T;
    name,
    ignore,
) where {T <: Field{<:AbstractData{<:Real}}}
    return _compare(parent(v1), parent(v2); name, ignore)
end

function _compare(pass, v1::T, v2::T; name, ignore) where {T <: AbstractData}
    return pass && _compare(parent(v1), parent(v2); name, ignore)
end

function _compare(
    v1::AbstractArray{FT},
    v2::AbstractArray{FT};
    name,
    ignore,
) where {FT <: AbstractFloat}
    error = maximum(_error(v1, v2))
    return print_maybe(error <= 100eps(eltype(v1)), "$name error: $error")
end

function _compare(pass, v1::T1, v2::T2; name, ignore) where {T1, T2}
    error("v1 and v2 have different types")
end

function print_maybe(exp, what)
    exp || println(what)
    return exp
end

# Begin tests

# Disable all the @info statements that are produced when creating a simulation
Logging.disable_logging(Logging.Info)

if comms_ctx isa ClimaComms.SingletonCommsContext
    configurations = ["sphere", "box", "column"]
else
    configurations = ["sphere", "box"]
end

for configuration in configurations
    if configuration == "sphere"
        moistures = ["equil", "nonequil"]
        precips = ["0M", "1M"]
        topography = "Earth"
        turbconv_models = [nothing, "diagnostic_edmfx"]
        # turbconv_models = ["prognostic_edmfx"]
        radiations = [nothing]
    else
        moistures = ["equil"]
        precips = ["1M"]
        topography = "NoWarp"
        turbconv_models = ["diagnostic_edmfx"]
        radiations = [nothing]
    end

    for turbconv_mode in turbconv_models
        for radiation in radiations
            for moisture in moistures
                for precip in precips
                    if !isnothing(turbconv_mode)
                        # EDMF only supports equilibrium moisture
                        if occursin("edmfx", turbconv_mode)
                            moisture == "equil" || continue
                        end
                    end

                    println(
                        "config = $configuration $moisture $precip $topography $radiation",
                    )
                    # The `enable_bubble` case is broken for ClimaCore < 0.14.6, so we
                    # hard-code this to be always false for those versions
                    bubble = pkgversion(ClimaCore) > v"0.14.5"

                    # Make sure that all MPI processes agree on the output_loc
                    output_loc =
                        ClimaComms.iamroot(comms_ctx) ? mktempdir(pwd()) : ""
                    output_loc = ClimaComms.bcast(comms_ctx, output_loc)
                    ClimaComms.barrier(comms_ctx)

                    job_id = "restart"
                    test_dict = Dict(
                        "test_dycore_consistency" => true, # We will add NaNs to the cache, just to make sure
                        "check_nan_every" => 3,
                        "log_progress" => false,
                        "moist" => moisture,
                        "precip_model" => precip,
                        "config" => configuration,
                        "topography" => topography,
                        # "turbconv" => turbconv_mode,
                        "dt" => "1secs",
                        "bubble" => bubble,
                        "viscous_sponge" => true,
                        "rayleigh_sponge" => true,
                        "insolation" => "timevarying",
                        "rad" => radiation,
                        "dt_rad" => "1secs",
                        "surface_setup" => "DefaultMoninObukhov",
                        # "implicit_diffusion" => true,
                        "call_cloud_diagnostics_per_stage" => true,  # Needed to ensure that cloud variables are computed
                        "t_end" => "3secs",
                        "dt_save_state_to_disk" => "1secs",
                        "enable_diagnostics" => false,
                        "output_dir" => joinpath(output_loc, job_id),
                    )
                    more_ignore = Symbol[]

                    config = CA.AtmosConfig(test_dict; job_id, comms_ctx)

                    simulation = CA.get_simulation(config)
                    CA.solve_atmos!(simulation)

                    # Check re-importing the same state
                    restart_dir = simulation.output_dir
                    @test isfile(joinpath(restart_dir), "day0.3.hdf5")

                    # Reset random seed for RRTMGP
                    Random.seed!(1234)

                    println("    just reading data")
                    if turbconv_mode == "prognostic_edmf"
                        more_ignore = [:ᶠnh_pressure₃ʲs]
                    end

                    config_should_be_same = CA.AtmosConfig(
                        merge(test_dict, Dict("detect_restart_file" => true));
                        job_id,
                        comms_ctx,
                    )

                    simulation_restarted =
                        CA.get_simulation(config_should_be_same)

                    SUCCESS[] &= compare(
                        simulation.integrator.u,
                        simulation_restarted.integrator.u;
                        name = "integrator.u",
                    )
                    SUCCESS[] &= compare(
                        axes(simulation.integrator.u.c),
                        axes(simulation_restarted.integrator.u.c);
                        name = "space",
                    )
                    SUCCESS[] &= compare(
                        simulation.integrator.p,
                        simulation_restarted.integrator.p;
                        name = "integrator.p",
                        ignore = Set([
                            :ghost_buffer,
                            :hyperdiffusion_ghost_buffer,
                            :scratch,
                            :output_dir,
                            :ghost_buffer,
                            # Computed in tendencies (which are not computed in this case)
                            :hyperdiff,
                            :precipitation,
                            # rc is some CUDA/CuArray internal object that we don't care about
                            :rc,
                            # Config-specific
                            more_ignore...,
                        ]),
                    )

                    # Check re-importing from previous state and advancing one step
                    println("    reading and simulating")
                    # Reset random seed for RRTMGP
                    Random.seed!(1234)

                    restart_file =
                        joinpath(simulation.output_dir, "day0.2.hdf5")
                    @test isfile(joinpath(restart_dir), "day0.2.hdf5")
                    # Restart from specific file
                    config2 = CA.AtmosConfig(
                        merge(test_dict, Dict("restart_file" => restart_file));
                        job_id,
                        comms_ctx,
                    )

                    simulation_restarted2 = CA.get_simulation(config2)
                    CA.fill_with_nans!(simulation_restarted2.integrator.p)

                    CA.solve_atmos!(simulation_restarted2)
                    SUCCESS[] &= compare(
                        simulation.integrator.u,
                        simulation_restarted2.integrator.u;
                        name = "integrator.u",
                    )
                    SUCCESS[] &= compare(
                        simulation.integrator.p,
                        simulation_restarted2.integrator.p;
                        name = "integrator.p",
                        ignore = Set([
                            :scratch,
                            :output_dir,
                            :ghost_buffer,
                            :hyperdiffusion_ghost_buffer,
                            :rc,
                        ]),
                    )
                end
            end
        end
    end
end

# Ensure that we have the correct exit code
@test SUCCESS[]
