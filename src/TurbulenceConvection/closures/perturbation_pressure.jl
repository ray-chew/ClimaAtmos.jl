"""
    compute_nh_pressure!(
        state::State,
        grid::Grid,
        edmf::EDMFModel,
        surf,
    )

Computes the
    - perturbation pressure gradient
    - (perturbation?) pressure advection
    - (perturbation?) pressure drag

for all updrafts, following [He2020](@cite), given:

 - `state`: state
 - `grid`: grid
 - `edmf`: EDMF model
 - `surf`: `SurfaceFluxes.SurfaceFluxConditions`
"""
function compute_nh_pressure!(state::State, grid::Grid, edmf::EDMFModel, surf)

    FT = float_type(state)
    N_up = n_updrafts(edmf)
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)

    Ifc = CCO.InterpolateF2C()
    wvec = CC.Geometry.WVector
    w_bcs =
        (; bottom = CCO.SetValue(wvec(FT(0))), top = CCO.SetValue(wvec(FT(0))))
    ∇ = CCO.DivergenceC2F(; w_bcs...)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    aux_gm_f = face_aux_grid_mean(state)
    aux_en_f = face_aux_environment(state)
    ρ_f = aux_gm_f.ρ
    plume_scale_height = ntuple(N_up) do i
        compute_plume_scale_height(grid, state, edmf.H_up_min, i)
    end

    # Note: Independence of aspect ratio hardcoded in implementation.
    α₂_asp_ratio² = FT(0)
    press_model_params = pressure_model_params(edmf)
    α_b = press_model_params.α_b
    α_a = press_model_params.α_a
    α_d = press_model_params.α_d

    @inbounds for i in 1:N_up
        # pressure
        b_up = aux_up_f[i].buoy
        w_up = prog_up_f[i].w
        H_up = plume_scale_height[i]
        w_en = aux_en_f.w

        nh_pressure = aux_up_f[i].nh_pressure

        @. nh_pressure =
            -α_b * b_up +
            α_a * wcomponent(CCG.WVector(w_up)) * ∇(Ifc(CCG.WVector(w_up))) -
            α_d *
            (wcomponent(CCG.WVector(w_up - w_en))) *
            abs(wcomponent(CCG.WVector(w_up - w_en))) / H_up
    end
    return nothing
end
