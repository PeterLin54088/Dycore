#################### >>> Import package >>> ####################
using JGCM
using PyPlot
#################### <<< Import package <<< ####################



#################### >>> Function overwrite >>> ####################
function modified_UV_Grid_From_Vor_Div!(mesh::Spectral_Spherical_Mesh, vor::Array{ComplexF64,3},  div::Array{ComplexF64,3}, grid_u::Array{Float64,3}, grid_v::Array{Float64,3})
    """
    Copy from original JGCM, but modified the cosine-weighted part.
    
    The original function assumes that vort and div is cosine-weighted, and 
    after compute ucos vcos, it'll divide cos.
    By assigning non-weighted vort and div, and not to divide cosine, the 
    modified function can output correct result.
    """
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
#################### <<< Function overwrite <<< ####################



#################### >>> Main >>> ####################
function Shallow_Water_Main()
    """
    """
    day_to_sec = 86400
    
    
    
    #################### >>> Model initialization >>> ####################
    """
    Construct necessary objects to instantiate model
    1. mesh        <-   Spectral_Spherical_Mesh
    2. atmo_data   <-   Atmo_Data
    3. dyn_data    <-   Dyn_Data
    4. integrator  <-   Filtered_Leapfrog
    5. op_man      <-   Output_Manager
    """
    model_name = "Shallow_Water"
    nλ, nθ, nd = 128, 64, 1
    num_fourier = floor(Int64, nθ*(2/3))
    num_spherical = num_fourier + 1
    radius = 6371.0e3
    omega = 7.292e-5
    # Initialize mesh
    mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical,
                                   nθ, nλ, nd, 
                                   radius)
    θc, λc = mesh.θc,  mesh.λc
    cosθ, sinθ = mesh.cosθ, mesh.sinθ
    sea_level_ps_ref = 1.0e5
    vert_coord = Vert_Coordinate(nλ, nθ, nd, 
                                 "even_sigma", 
                                 "simmons_and_burridge", 
                                 "second_centered_wts", 
                                 sea_level_ps_ref)
    # Initialize atmo_data
    atmo_data = Atmo_Data(model_name, 
                          nλ, nθ, nd, 
                          false, false, false, false, 
                          sinθ, radius, omega)
    # Initialize dyn_data
    dyn_data = Dyn_Data(model_name,
                        num_fourier, num_spherical,
                        nλ, nθ, nd)
    # Initialize integrator
        # Robert-Asselin filter parameter
    damping_order = 10
    damping_coef = 1.0e-10
    robert_coef  = 0.04
        # Leapfrog parameter
    init_step = true
    implicit_coef = 0.5
        # Spinup period (day)
    spinup_day = 1
        # Time
    start_time = 0
    end_time = 30*day_to_sec
    Δt = 600
    NT =  Int64(end_time / Δt)
    
    integrator = Filtered_Leapfrog(robert_coef, damping_order, damping_coef, mesh.laplacian_eig,
                                   implicit_coef, Δt, init_step, start_time, end_time)
    # Data Visualization
    op_man = Output_Manager(mesh, 
                            vert_coord, 
                            start_time, end_time, 
                            spinup_day)
    #################### <<< Model initialization <<< ####################
    
    
    
    #################### >>> Field initialization >>> ####################
    """
    Some facts should know before modify this block
    1.The governing equations are
      ∂vor/∂t = ∇ ((f + vor) V)  := δvor
      ∂div/∂t = ∇ × ((f + vor) V) - ∇^2(E + H) := δdiv
      ∂H/∂t = -H div - ∇H V  := δH
    2.Method {Vor_Div_From_Grid_UV!} and {UV_Grid_From_Vor_Div!} are truncated,
      which produces gibbs phenomenon in fixed boundary
    ***background vorticity and divergence***
    By default, the vorticity and divergence in spectral space are cosine-weighted,
    and if we need to assign vorticity and divergence in grid space, we have to 
    manually weight them.
    An alternative method is that, calculate non-weighted vorticity and divergence,
    but add cosine-weight in our overwrited transfer function
    1. assign not-weighting vorticity and divergence 
    2. calculate not-weighting u and v wind
    3. re-weighting by self-defined function
    4. calculate vorticity and divergence in spherical coordinate
    
    ***background balanced height***
    In original SWE, the coriolis and pressure should reach equilibrium, and in this
    model, this relation is presented within divergence equation, describe as
    geopotential = invlap(absolute vorticity advection) - kinetic energy 
    ***perturbed vorticity and divergence***
    
    ***perturbed height***
    
    ***perturbed balanced wind***
    """
    grid_u, grid_v = dyn_data.grid_u_c, dyn_data.grid_v_c
    grid_vor, grid_div = dyn_data.grid_vor, dyn_data.grid_div
    spe_vor_c, spe_div_c = dyn_data.spe_vor_c, dyn_data.spe_div_c
    grid_h, spe_h_c = dyn_data.grid_ps_c, dyn_data.spe_lnps_c
    
    grid_absvor = dyn_data.grid_absvor
    grid_kin, spe_kin = dyn_data.grid_d_full1, dyn_data.spe_d1
    spe_energy = dyn_data.spe_energy
    
    grid_δu = dyn_data.grid_δu
    grid_δv = dyn_data.grid_δv
    spe_δvor = dyn_data.spe_δvor
    spe_δdiv = dyn_data.spe_δdiv
    
    MEAN_HEIGHT = 2.0e3 * 9.81 # Geopotential height
    
    # Background vort & div (u & v)
    grid_u_b, grid_v_b = zeros(Float64, size(grid_u)), zeros(Float64, size(grid_v))
    grid_vor_b, grid_div_b = zeros(Float64, size(grid_vor)), zeros(Float64, size(grid_div))
    spe_vor_b, spe_div_b = zeros(ComplexF64, size(spe_vor_c)), zeros(ComplexF64, size(spe_div_c))
        # Vorticity strip
    xx = reshape(ones(nθ)' .* λc, (nλ,nθ,nd))
    yy = reshape(θc' .* ones(nλ), (nλ,nθ,nd))
    grid_vor_b .= 4e-5 * exp.(-(((yy .- 19*pi/180)/(3*pi/180)).^2))
    grid_div_b .= 0
        # Cosine-weighting removal
    Trans_Grid_To_Spherical!(mesh, grid_vor_b, spe_vor_b)
    Trans_Grid_To_Spherical!(mesh, grid_div_b, spe_div_b)
    modified_UV_Grid_From_Vor_Div!(mesh, spe_vor_b, spe_div_b, grid_u_b, grid_v_b)
        # Smoothen signal at boundary
    grid_u_b .= grid_u_b .* exp.(-(((yy .- 19*pi/180)/(40*pi/180)).^2))
    Vor_Div_From_Grid_UV!(mesh, grid_u_b, grid_v_b, spe_vor_b, spe_div_b)
    UV_Grid_From_Vor_Div!(mesh, spe_vor_b, spe_div_b, grid_u_b, grid_v_b)
    Trans_Spherical_To_Grid!(mesh, spe_vor_b,  grid_vor_b)
    Trans_Spherical_To_Grid!(mesh, spe_div_b,  grid_div_b)
    
    UV_Grid_From_Vor_Div!(mesh, spe_vor_b, spe_div_b, grid_u_b, grid_v_b)
    
    println("background vorticity amplitude : ", maximum(grid_vor_b))
    # Background balanced height
    grid_h_b, spe_h_b = zeros(Float64, size(grid_h)), zeros(ComplexF64, size(spe_h_c))  
        # Absolute vorticity advection
    Compute_Abs_Vor!(grid_vor_b, atmo_data.coriolis, grid_absvor)
    grid_δu .=   grid_absvor .* grid_v_b
    grid_δv .=  -grid_absvor .* grid_u_b
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)
        # Kinetic energy
    Apply_InverseLaplacian!(mesh, spe_δdiv)
    grid_kin .= 0.5 * (grid_u_b.^2 + grid_v_b.^2)
    Trans_Grid_To_Spherical!(mesh, grid_kin, spe_kin)
        # Potential energy
    spe_h_b .= spe_δdiv - spe_kin
    println("background geopotential anomaly adjustment : ", spe_h_b[1,1])
    spe_h_b[1,1] = MEAN_HEIGHT # set mean height
    Trans_Spherical_To_Grid!(mesh, spe_h_b,  grid_h_b)
    
    # Perturbed vort & div (u & v)
    grid_u_p, grid_v_p = zeros(Float64, size(grid_u)), zeros(Float64, size(grid_v))
    grid_vor_p, grid_div_p = zeros(Float64, size(grid_vor)), zeros(Float64, size(grid_div))
    spe_vor_p, spe_div_p = zeros(ComplexF64, size(spe_vor_c)), zeros(ComplexF64, size(spe_div_c))
        # Vorticity blob
    xx = reshape(ones(nθ)' .* λc, (nλ,nθ,nd))
    yy = reshape(θc' .* ones(nλ), (nλ,nθ,nd))
    grid_vor_p .= 3e-6 * sin.(15*xx) .* exp.(-(((yy .- 20*pi/180)/(5*pi/180)).^2))
    grid_div_p .= 0
        # Cosine-weighting removal
    Trans_Grid_To_Spherical!(mesh, grid_vor_p, spe_vor_p)
    Trans_Grid_To_Spherical!(mesh, grid_div_p, spe_div_p)
    modified_UV_Grid_From_Vor_Div!(mesh, spe_vor_p, spe_div_p, grid_u_p, grid_v_p)
    Vor_Div_From_Grid_UV!(mesh, grid_u_p, grid_v_p, spe_vor_p, spe_div_p)
    Trans_Spherical_To_Grid!(mesh, spe_vor_p,  grid_vor_p)
    Trans_Spherical_To_Grid!(mesh, spe_div_p,  grid_div_p)

    # PyPlot.figure(figsize = (16,9), dpi = 160)
    # PyPlot.contourf(grid_vor_p[:,:,1]', levels = 32, extend = "both", cmap = "bwr")
    # PyPlot.savefig("u.png")
    
    # Perturbed balanced height
    grid_h_p, spe_h_p = zeros(Float64, size(grid_h)), zeros(ComplexF64, size(spe_h_c))
        # Absolute vorticity advection
    Compute_Abs_Vor!(grid_vor_p, atmo_data.coriolis, grid_absvor)
    grid_δu .=   grid_absvor .* grid_v_p
    grid_δv .=  -grid_absvor .* grid_u_p
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)
        # Kinetic energy
    Apply_InverseLaplacian!(mesh, spe_δdiv)
    grid_kin .= 0.5 * (grid_u_p.^2 + grid_v_p.^2)
    Trans_Grid_To_Spherical!(mesh, grid_kin, spe_kin)
        # Potential energy
    spe_h_p .= spe_δdiv - spe_kin
    println("perturbed geopotential distance is : ",spe_h_p[1,1])
    spe_h_p[1,1] = 0 # perturbed must be zero-mean
    Trans_Spherical_To_Grid!(mesh, spe_h_p,  grid_h_p)
    
    # Reference height
    h_eq = zeros(Float64, nλ, nθ, 1)
    h_eq .= MEAN_HEIGHT
    # Complete field
    grid_u .= grid_u_b + grid_u_p
    grid_v .= grid_v_b + grid_v_p
    grid_vor .= grid_vor_b + grid_vor_p
    grid_div .= grid_div_b + grid_div_p
    grid_h .= grid_h_b + grid_h_p
    spe_vor_c .= spe_vor_b + spe_vor_p
    spe_div_c .= spe_div_b + spe_div_p
    spe_h_c .= spe_h_b + spe_h_p
    
    # grid_u .= grid_u_b
    # grid_v .= grid_v_b
    # grid_vor .= grid_vor_b
    # grid_div .= grid_div_b
    # grid_h .= grid_h_b
    # spe_vor_c .= spe_vor_b
    # spe_div_c .= spe_div_b
    # spe_h_c .= spe_h_b
    # Output initial field
    Update_Output_Init!(op_man, dyn_data, integrator.time)
    #################### <<< Field initialization <<< ####################
    
    
    
    #################### >>> Parameter initialization >>> ####################
    # Linear damping
    fric_damp_time  = 1e5 * day_to_sec
    therm_damp_time = 1e5 * day_to_sec
    kappa_m = 1.0 / fric_damp_time
    kappa_t = 1.0 / therm_damp_time
    #################### <<< Parameter initialization <<< ####################
    
    
    
    #################### >>> Simulation process >>> ####################
    Shallow_Water_Physics!(dyn_data, kappa_m, kappa_t, h_eq)
    Shallow_Water_Dynamics!(mesh, atmo_data, MEAN_HEIGHT, dyn_data, integrator)
    Update_Init_Step!(integrator)
    integrator.time += Δt
    Update_Output!(op_man, dyn_data, integrator.time)
    for i = 2:NT
        Shallow_Water_Physics!(dyn_data, kappa_m, kappa_t, h_eq)
        Shallow_Water_Dynamics!(mesh, atmo_data, MEAN_HEIGHT, dyn_data, integrator)
        integrator.time += Δt 
        Update_Output!(op_man, dyn_data, integrator.time)
        # PyPlot.figure(figsize = (16,9), dpi = 160)
        # PyPlot.contourf(grid_u[:,:,1]', levels = LinRange(-10,10,21), extend = "both", cmap = "bwr")
        # PyPlot.savefig("u_$i")
        # PyPlot.figure(figsize = (16,9), dpi = 160)
        # PyPlot.contourf(grid_vor[:,:,1]', levels = LinRange(-4e-5,4e-5,21), extend = "both", cmap = "bwr")
        # PyPlot.savefig("vor_$i")
        if (integrator.time%day_to_sec == 0)
            @info (÷(integrator.time, day_to_sec))
        end
    end
    #################### <<< Simulation process <<< ####################
    return op_man
  end
#################### <<< Main <<< ####################