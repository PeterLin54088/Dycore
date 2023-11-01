#####################################################################################################
using JGCM
using PyPlot
using LinearAlgebra
#####################################################################################################

#####################################################################################################
function modified_UV_Grid_From_Vor_Div!(mesh::Spectral_Spherical_Mesh, vor::Array{ComplexF64,3},  div::Array{ComplexF64,3}, grid_u::Array{Float64,3}, grid_v::Array{Float64,3})
    nd = mesh.nd
    @assert(size(vor)[3] == nd || size(vor)[3] == 1)
    if (size(vor)[3] == nd)
        spherical_ucos, spherical_vcos = mesh.spherical_d1, mesh.spherical_d2
    else
        spherical_ucos, spherical_vcos = mesh.spherical_ds1, mesh.spherical_ds2
    end  
    Compute_Ucos_Vcos_From_Vor_Div!(mesh, vor, div, spherical_ucos, spherical_vcos)
    
    Trans_Spherical_To_Grid!(mesh, spherical_ucos, grid_u)
    Trans_Spherical_To_Grid!(mesh, spherical_vcos, grid_v)
    
    cosθ = mesh.cosθ
    
    # Without cosine weighting
    
    # Divide_By_Cos!(cosθ, grid_u)
    # Divide_By_Cos!(cosθ, grid_v)
    
end
#####################################################################################################

function Shallow_Water_Main()
    """
    """
    
    name = "Shallow_Water"
    num_fourier, nθ, nd = 85, 128, 1
    num_spherical = num_fourier + 1
    nλ = 2nθ
    
    radius = 6371.0e3
    omega = 7.292e-5
    day_to_sec = 86400
    fric_damp_time  = 20.0 * day_to_sec
    therm_damp_time = 10.0 * day_to_sec
    
    h_0             = 1.0e4
    
    kappa_m = 1.0 / fric_damp_time
    kappa_t = 1.0 / therm_damp_time
    sea_level_ps_ref = 1.0e5
    
    # Initialize mesh
    mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nθ, nλ, nd, radius)
    θc, λc = mesh.θc,  mesh.λc
    cosθ, sinθ = mesh.cosθ, mesh.sinθ
    vert_coord = Vert_Coordinate(nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts", sea_level_ps_ref)
    
  
    # Initialize atmo_data
    atmo_data = Atmo_Data(name, nλ, nθ, nd, false, false, false, false, sinθ, radius, omega)
    
    # Initialize integrator
    damping_order = 4
    damping_coef = 1.e-04
    robert_coef  = 0.04 
    
    implicit_coef = 0.5
    
    day_to_sec = 86400
    start_time = 0
    end_time = 2*day_to_sec
    Δt = 600
    spinup_day = 1
    init_step = true
    
    integrator = Filtered_Leapfrog(robert_coef, 
    damping_order, damping_coef, mesh.laplacian_eig,
    implicit_coef,
    Δt, init_step, start_time, end_time)
    
    # Data Visualization
    op_man = Output_Manager(mesh, vert_coord, start_time, end_time, spinup_day)
    
    # Initialize data
    dyn_data = Dyn_Data(name, num_fourier, num_spherical, nλ, nθ, nd)
     
    # Initialize vorticity and divergence & uv wind
    """
    The model operates vorticity and divergence in spherical coordinate, not in grid space
    1. assign not-weighting vorticity and divergence 
    2. calculate not-weighting u and v wind
    3. re-weighting by self-defined function
    4. calculate vorticity and divergence in spherical coordinate
    """
    grid_vor, grid_div = dyn_data.grid_vor, dyn_data.grid_div
    spe_vor_c, spe_div_c = dyn_data.spe_vor_c, dyn_data.spe_div_c
    grid_u, grid_v = dyn_data.grid_u_c, dyn_data.grid_v_c
    
    # Temporary variable
    temp_grid_vor, temp_grid_div, temp_spe_vor, temp_spe_div = zeros(Float64, size(grid_vor)), zeros(Float64, size(grid_div)), zeros(ComplexF64, size(spe_vor_c)), zeros(ComplexF64, size(spe_div_c))
    
    # Vorticity strip
    xx = reshape(ones(nθ)' .* λc, (nλ,nθ,nd))
    yy = reshape(θc' .* ones(nλ), (nλ,nθ,nd))
    temp_grid_vor .= 4e-5 * exp.(-(((yy .- 19*pi/180)/(3*pi/180)).^2))
    temp_grid_div .= 0
    
    # Non-weighting UV wind & re-weighting
    Trans_Grid_To_Spherical!(mesh, temp_grid_vor, temp_spe_vor)
    Trans_Grid_To_Spherical!(mesh, temp_grid_div, temp_spe_div)
    modified_UV_Grid_From_Vor_Div!(mesh, temp_spe_vor, temp_spe_div, grid_u, grid_v)
    
    # spherical-weighted vorticity and divergence
    Vor_Div_From_Grid_UV!(mesh, grid_u, grid_v, spe_vor_c, spe_div_c) 
    Trans_Spherical_To_Grid!(mesh, spe_vor_c,  grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_c,  grid_div)
    
        # height
    grid_h, spe_h_c = dyn_data.grid_ps_c, dyn_data.spe_lnps_c
    grid_absvor = dyn_data.grid_absvor
    
    grid_vor .= 0
    grid_u .= 10
    
    Compute_Abs_Vor!(grid_vor, atmo_data.coriolis, grid_absvor)
    grid_δu = dyn_data.grid_δu
    grid_δv = dyn_data.grid_δv
    spe_δvor = dyn_data.spe_δvor
    spe_δdiv = dyn_data.spe_δdiv
    grid_δu .+=   grid_absvor .* grid_v
    grid_δv .+=  -grid_absvor .* grid_u
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)
    
    grid_kin, spe_kin = dyn_data.grid_d_full1, dyn_data.spe_d1
    spe_energy = dyn_data.spe_energy
    grid_kin .= 0.5 * (grid_u.^2 + grid_v.^2)
    Trans_Grid_To_Spherical!(mesh, grid_kin, spe_kin)
    Apply_InverseLaplacian!(mesh, spe_δdiv)
    spe_h_c .= spe_δdiv - spe_kin
    # spe_h_c[1,1] = 0
    Trans_Spherical_To_Grid!(mesh, spe_h_c,  grid_h)
    println(grid_h)
    PyPlot.figure(figsize = (16,9), dpi = 160)
    PyPlot.contourf(grid_h[:,:,1]' , levels = LinRange(-9800,9800,21), cmap = "bwr")
    PyPlot.xlim(0,255)
    PyPlot.xticks([0,63,127,191,255])
    PyPlot.ylim(0,127)
    PyPlot.yticks([0,31,63,95,127])
    PyPlot.savefig("test.png")
#     """
#     balanced height using thermal wind relation is required
#     """
#     grid_h .= h_0
#     Trans_Grid_To_Spherical!(mesh, grid_h, spe_h_c)
    
    
        
#     h_eq = zeros(Float64, nλ, nθ, 1)
#     h_eq .= h_0
  
#     NT =  Int64(end_time / Δt)
    
#     Update_Output_Init!(op_man, dyn_data, integrator.time)
    
#     Shallow_Water_Physics!(dyn_data, kappa_m, kappa_t, h_eq)
#     Shallow_Water_Dynamics!(mesh, atmo_data, h_0, dyn_data, integrator)
#     Update_Init_Step!(integrator)
#     integrator.time += Δt
#     Update_Output!(op_man, dyn_data, integrator.time)

    
#     for i = 2:NT
#         Shallow_Water_Physics!(dyn_data, kappa_m, kappa_t, h_eq)
#         Shallow_Water_Dynamics!(mesh, atmo_data, h_0, dyn_data, integrator)
#         integrator.time += Δt 
#         Update_Output!(op_man, dyn_data, integrator.time)
#         # PyPlot.figure(figsize = (16,9), dpi = 160)
#         # PyPlot.contourf(grid_u[:,:,1]', levels = LinRange(-10,10,21), extend = "both", cmap = "bwr")
#         # PyPlot.savefig("u_$i")
#         # PyPlot.figure(figsize = (16,9), dpi = 160)
#         # PyPlot.contourf(grid_vor[:,:,1]', levels = LinRange(-4e-5,4e-5,21), extend = "both", cmap = "bwr")
#         # PyPlot.savefig("vor_$i")
#         if (integrator.time%day_to_sec == 0)
#             @info (÷(integrator.time, day_to_sec))
#         end
#     end
#     return op_man
  end