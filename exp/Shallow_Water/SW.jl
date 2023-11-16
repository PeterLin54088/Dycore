#################### >>> Import package >>> ####################
using JGCM
#################### <<< Import package <<< ####################



#################### >>> Function rewrite >>> ####################
function modified_UV_Grid_From_Vor_Div!(mesh::Spectral_Spherical_Mesh, vor::Array{ComplexF64,3},  div::Array{ComplexF64,3}, grid_u::Array{Float64,3}, grid_v::Array{Float64,3})
    """
    Reference from https://www.gfdl.noaa.gov/idealized-spectral-models-quickstart/
    
    **********************************************************************************************
    subroutine_uv_grid_from_vor_div : 
    Takes as input vorticity and divergence in the spectral domain and returns (u, v) on the grid. 
    It first computes (cos(θ)u, cos(θ)v), next calculate vor and div in the spectral domain, 
    then transforms back to the grid and divides by cos(θ). 
    
    subroutine_vor_div_from_uv_grid : 
    Revert process in subroutine_uv_grid_from_vor_div
    **********************************************************************************************
    
    Above are official documentation and instruction about how to implement transformation
    in UV and Vor&Div, however, the original author does not follow the instruction and only
    provide uv_grid as the only legal access to the initial field.
    
    This self-defined function takes as input vorticity and divergence "without cosine-weighted",
    and can output correct UV in grid space.
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

function Background_Vorticity_Strip!(mesh::Spectral_Spherical_Mesh, 
                                     atmo_data::Atmo_Data,
                                     grid_u::Array{Float64,3}, grid_v::Array{Float64,3}, 
                                     spe_vor_b::Array{ComplexF64,3}, spe_div_b::Array{ComplexF64,3},
                                     spe_h_b::Array{ComplexF64,3},
                                     vor_pert::Float64, vor_latitude::Float64, vor_width::Float64,
                                     Mean_Height::Float64)
    """
    Generate a vorticity strip following Suhas & Boos
    Also, generate the corresponding balanced height via geostrophic balance(?)
    """
    nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd
    θc, λc = mesh.θc, mesh.λc
    deg_to_rad = pi/180
    
    # Vorticity strip
        # Assign vor and div
    grid_u_b, grid_v_b = zeros(Float64, size(grid_u)), zeros(Float64, size(grid_v)) 
    lons = reshape(ones(nθ)' .* λc, (nλ,nθ,nd))
    lats = reshape(θc' .* ones(nλ), (nλ,nθ,nd))
    grid_vor_b = vor_pert * exp.(-((lats .- vor_latitude*deg_to_rad)/(vor_width*deg_to_rad)).^2)
    grid_div_b = 0.0 .* lons
        # Corresponding wind field
        # Cosine-weighting removal
    Trans_Grid_To_Spherical!(mesh, grid_vor_b, spe_vor_b)
    Trans_Grid_To_Spherical!(mesh, grid_div_b, spe_div_b)
    modified_UV_Grid_From_Vor_Div!(mesh, spe_vor_b, spe_div_b, grid_u_b, grid_v_b)
        # Avoid gibbs phenomenon (associated Legendre polynomials)
    grid_u_b .= grid_u_b .* exp.(-((lats .- vor_latitude*deg_to_rad)/(10*vor_width*deg_to_rad)).^2)
    Vor_Div_From_Grid_UV!(mesh, grid_u_b, grid_v_b, spe_vor_b, spe_div_b)
    Trans_Spherical_To_Grid!(mesh, spe_vor_b,  grid_vor_b)
    Trans_Spherical_To_Grid!(mesh, spe_div_b,  grid_div_b)
    
    # Background balanced height
        # Temporary variable
    grid_absvor = zeros(Float64, size(grid_vor_b))
    spe_δvor = zeros(ComplexF64, size(spe_vor_b))
    spe_δdiv = zeros(ComplexF64, size(spe_div_b))
    spe_kin = zeros(ComplexF64, size(spe_vor_b))
        # Absolute vorticity advection
    Compute_Abs_Vor!(grid_vor_b, atmo_data.coriolis, grid_absvor)
    grid_δu =  grid_absvor .* grid_v_b
    grid_δv = -grid_absvor .* grid_u_b
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)
        # Kinetic energy
    Apply_InverseLaplacian!(mesh, spe_δdiv)
    grid_kin = 0.5 * (grid_u_b.^2 + grid_v_b.^2)
    Trans_Grid_To_Spherical!(mesh, grid_kin, spe_kin)
        # Potential energy
    spe_h_b .= spe_δdiv - spe_kin
    spe_h_b[1,1] = Mean_Height
end

function Perturbed_Vorticity_Blobs!(mesh::Spectral_Spherical_Mesh,
                                    atmo_data::Atmo_Data,
                                    grid_u::Array{Float64,3}, grid_v::Array{Float64,3}, 
                                    spe_vor_p::Array{ComplexF64,3}, spe_div_p::Array{ComplexF64,3},
                                    spe_h_p::Array{ComplexF64,3},
                                    vor_pert::Float64, vor_latitude::Float64, vor_width::Float64,
                                    Mean_Height::Float64)
    """
    Generate small-amplitude sinusoidal vorticity anomalies following Suhas & Boos
    Also, generate the corresponding balanced height via geostrophic balance(?)
    """
     
    nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd
    θc, λc = mesh.θc, mesh.λc
    deg_to_rad = pi/180
    
    # Vorticity blobs
        # Assign vor and div
    grid_u_p, grid_v_p = zeros(Float64, size(grid_u)), zeros(Float64, size(grid_v))
    lons = reshape(ones(nθ)' .* λc, (nλ,nθ,nd))
    lats = reshape(θc' .* ones(nλ), (nλ,nθ,nd))
    grid_vor_p = vor_pert * sin.(15*lons) .* exp.(-((lats .- vor_latitude*deg_to_rad)/(vor_width*deg_to_rad)).^2)
    grid_div_p = 0.0 .* lons
        # Corresponding wind field
        # Cosine-weighting removal
    Trans_Grid_To_Spherical!(mesh, grid_vor_p, spe_vor_p)
    Trans_Grid_To_Spherical!(mesh, grid_div_p, spe_div_p)
    modified_UV_Grid_From_Vor_Div!(mesh, spe_vor_p, spe_div_p, grid_u_p, grid_v_p)
    Vor_Div_From_Grid_UV!(mesh, grid_u_p, grid_v_p, spe_vor_p, spe_div_p)
    
    # Background balanced height
        # Temporary variable
    grid_absvor = zeros(Float64, size(grid_vor_p))
    spe_δvor = zeros(ComplexF64, size(spe_vor_p))
    spe_δdiv = zeros(ComplexF64, size(spe_div_p))
    spe_kin = zeros(ComplexF64, size(spe_vor_p))
        # Absolute vorticity advection
    Compute_Abs_Vor!(grid_vor_p, atmo_data.coriolis, grid_absvor)
    grid_δu =  grid_absvor .* grid_v_p
    grid_δv = -grid_absvor .* grid_u_p
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)
        # Kinetic energy
    Apply_InverseLaplacian!(mesh, spe_δdiv)
    grid_kin = 0.5 * (grid_u_p.^2 + grid_v_p.^2)
    Trans_Grid_To_Spherical!(mesh, grid_kin, spe_kin)
        # Potential energy
    spe_h_p .= spe_δdiv - spe_kin
    spe_h_p[1,1] = Mean_Height
end
#################### <<< Function rewrite <<< ####################



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
    nλ, nθ, nd = 256, 128, 1
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
        # Hyper-viscosity parameter
    damping_order = 8
    damping_coef = 1.0e-5
        # Robert-Asselin filter parameter
    robert_coef  = 0.04
        # Leapfrog parameter
    init_step = true
    implicit_coef = 0.5 # centered
        # Spinup period (day)
    spinup_day = 1
        # Time
    start_time = 0
    end_time = 10*day_to_sec
    Δt = 3600 # CFL condition needed
    NT = Int64(end_time / Δt)
    
    integrator = Filtered_Leapfrog(robert_coef, 
                                   damping_order, damping_coef, 
                                   mesh.laplacian_eig, implicit_coef, 
                                   Δt, init_step, start_time, end_time)
    # Data Visualization
    op_man = Output_Manager(mesh, 
                            vert_coord, 
                            start_time, end_time, spinup_day)
    #################### <<< Model initialization <<< ####################
    
    
    
    #################### >>> Field initialization >>> ####################
    """
    The governing equations are
        ∂vor/∂t = ∇ ((f + vor) V)  := δvor
        ∂div/∂t = ∇ × ((f + vor) V) - ∇^2(E + H) := δdiv
        ∂H/∂t = -H div - ∇H V  := δH
    Gibbs phenomenon can occur when signals near fixed boundary (ex. poles) are not zeros
    """
    grid_u, grid_v = dyn_data.grid_u_c, dyn_data.grid_v_c
    grid_vor, spe_vor_c = dyn_data.grid_vor, dyn_data.spe_vor_c
    grid_div, spe_div_c = dyn_data.grid_div, dyn_data.spe_div_c
    grid_h, spe_h_c = dyn_data.grid_ps_c, dyn_data.spe_lnps_c
    
    # Background field
    spe_vor_b, spe_div_b = zeros(ComplexF64, size(spe_vor_c)), zeros(ComplexF64, size(spe_div_c))
    spe_h_b = zeros(ComplexF64, size(spe_h_c))
    
    #-------------------------------------------------------------------------------------#
    vor_pert = 5.37e-5
    vor_latitude = 19.0
    vor_width = 3.0
    MEAN_HEIGHT = 1.0e3 * 9.81
    #-------------------------------------------------------------------------------------#
    
    Background_Vorticity_Strip!(mesh, 
                                atmo_data,
                                grid_u, grid_v,
                                spe_vor_b, spe_div_b,
                                spe_h_b,
                                vor_pert, vor_latitude, vor_width,
                                MEAN_HEIGHT)
    
    # Perturbed vort & div (u & v)
    spe_vor_p, spe_div_p = zeros(ComplexF64, size(spe_vor_c)), zeros(ComplexF64, size(spe_div_c))
    spe_h_p = zeros(ComplexF64, size(spe_h_c))
    
    #-------------------------------------------------------------------------------------#
    vor_pert = 1e-7
    vor_latitude = 20.0
    vor_width = 5.0
    ANOMALY_HEIGHT = 0.0e3 * 9.81
    #-------------------------------------------------------------------------------------#
    
    Perturbed_Vorticity_Blobs!(mesh, 
                               atmo_data,
                               grid_u, grid_v,
                               spe_vor_p, spe_div_p,
                               spe_h_p,
                               vor_pert, vor_latitude, vor_width,
                               ANOMALY_HEIGHT)
   
    # Reference height
    h_eq = zeros(Float64, nλ, nθ, 1)
    h_eq .= MEAN_HEIGHT
    
    # Complete field
    spe_vor_c .= spe_vor_b + spe_vor_p
    spe_div_c .= spe_div_b + spe_div_p
    spe_h_c .= spe_h_b + spe_h_p
    
    Trans_Spherical_To_Grid!(mesh, spe_vor_c, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_c, grid_div)
    Trans_Spherical_To_Grid!(mesh, spe_h_c, grid_h)
    UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
    
    println(repeat("###", 30))
    println("Initial       vorticity maximum : ", round(maximum(grid_vor); digits = 9), "  (1/s)")
    println("Initial      zonal wind maximum : ", round(maximum(grid_u); digits = 4), "  (m/s)")
    println("Initial meridional wind maximum : ", round(maximum(grid_v); digits = 4), "  (m/s)")
    println(repeat("###", 30))
    # Output initial field
    Update_Output_Init!(op_man, dyn_data, integrator.time)
    #################### <<< Field initialization <<< ####################
    
    
    
    #################### >>> Parameter initialization >>> ####################
    """
    Linear damping
    Following Suhas & Boos, the physical damping is not included
    """
    
    # fric_damp_time  = 20.0 * day_to_sec
    # therm_damp_time = 10.0 * day_to_sec
    kappa_m = 0.0
    kappa_t = 0.0
    #################### <<< Parameter initialization <<< ####################
    
    
    
    #################### >>> Simulation process >>> ####################
    Shallow_Water_Physics!(dyn_data, kappa_m, kappa_t, h_eq)
    Shallow_Water_Dynamics!(mesh, atmo_data, dyn_data, integrator, MEAN_HEIGHT)
    Update_Init_Step!(integrator)
    integrator.time += Δt
    Update_Output!(op_man, dyn_data, integrator.time)
    for i = 2:NT
        Shallow_Water_Physics!(dyn_data, kappa_m, kappa_t, h_eq)
        Shallow_Water_Dynamics!(mesh, atmo_data, dyn_data, integrator, MEAN_HEIGHT)
        integrator.time += Δt 
        Update_Output!(op_man, dyn_data, integrator.time)
        if (integrator.time%day_to_sec == 0)
            @info (÷(integrator.time, day_to_sec))
        end
    end
    #################### <<< Simulation process <<< ####################
    return op_man
  end
#################### <<< Main <<< ####################