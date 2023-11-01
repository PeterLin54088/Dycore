export Output_Manager, Update_Output!, Update_Output_Init!, Finalize_Output!
export Lat_Lon_Pcolormesh, Zonal_Mean, Sigma_Zonal_Mean_Pcolormesh, Sigma_Zonal_Mean_Contourf
                             
mutable struct Output_Manager
    nλ::Int64
    nθ::Int64
    nd::Int64
    n_day::Int64
    
    day_to_sec::Int64
    start_time::Int64
    end_time::Int64
    current_time::Int64
    spinup_day::Int64

    λc::Array{Float64, 1}
    θc::Array{Float64, 1}
    σc::Array{Float64, 1}
  
    # nλ × nθ × nd × n_day
    # The average is (start, end], namely it does not include the first snapshot.
    u_daily_mean::Array{Float64, 4}
    v_daily_mean::Array{Float64, 4}
    h_daily_mean::Array{Float64, 4}
    
    # n_day
    n_daily_mean::Array{Float64, 1}
    
end

function Output_Manager(mesh::Spectral_Spherical_Mesh, vert_coord::Vert_Coordinate, start_time::Int64, end_time::Int64, spinup_day::Int64)
    nλ = mesh.nλ
    nθ = mesh.nθ
    nd = mesh.nd
    
    day_to_sec = 86400
    current_time = start_time

    λc = mesh.λc
    θc = mesh.θc

    #todo definition of sigma coordinate
    bk = vert_coord.bk
    σc = (bk[2:nd+1] + bk[1:nd])/2.0
    
    n_day = Int64((end_time - start_time)/day_to_sec) + 1
    
    # nλ × nθ × nd × n_day
    # The average is (start, end], namely it does not include the first snapshot.
    u_daily_mean = zeros(Float64, nλ, nθ, nd, n_day)
    v_daily_mean = zeros(Float64, nλ, nθ, nd, n_day)
    h_daily_mean = zeros(Float64, nλ, nθ, 1, n_day)
    n_daily_mean = zeros(Float64, n_day)
    
    Output_Manager(nλ, nθ, nd, n_day,
    day_to_sec, start_time, end_time, current_time, spinup_day,
    λc, θc, σc,
    u_daily_mean, v_daily_mean, h_daily_mean,
    n_daily_mean)
end

function Update_Output_Init!(output_manager::Output_Manager, dyn_data::Dyn_Data, current_time::Int64)
    @assert(current_time == output_manager.current_time)
    day_to_sec, start_time, n_day = output_manager.day_to_sec, output_manager.start_time, output_manager.n_day

    u_daily_mean, v_daily_mean, h_daily_mean = output_manager.u_daily_mean, output_manager.v_daily_mean, output_manager.h_daily_mean
    n_daily_mean = output_manager.n_daily_mean


    i_day = 1

    # if(i_day > n_day)
    #     @info "Warning: i_day > n_day in Output_Manager:Update!"
    #     return 
    # end
    
    u_daily_mean[:,:,:,i_day] .+= dyn_data.grid_u_c
    v_daily_mean[:,:,:,i_day] .+= dyn_data.grid_v_c
    h_daily_mean[:,:,1,i_day] .+= dyn_data.grid_ps_c[:,:,1]
    n_daily_mean[i_day] += 1
end

function Update_Output!(output_manager::Output_Manager, dyn_data::Dyn_Data, current_time::Int64)
    @assert(current_time > output_manager.current_time)
    output_manager.current_time = current_time
    day_to_sec, start_time, n_day = output_manager.day_to_sec, output_manager.start_time, output_manager.n_day

    u_daily_mean, v_daily_mean, h_daily_mean = output_manager.u_daily_mean, output_manager.v_daily_mean, output_manager.h_daily_mean
    n_daily_mean = output_manager.n_daily_mean


    i_day = div(current_time - start_time - 1, day_to_sec) + 1 + 1

    # if(i_day > n_day)
    #     @info "Warning: i_day > n_day in Output_Manager:Update!"
    #     return 
    # end
    
    u_daily_mean[:,:,:,i_day] .+= dyn_data.grid_u_c
    v_daily_mean[:,:,:,i_day] .+= dyn_data.grid_v_c
    h_daily_mean[:,:,1,i_day] .+= dyn_data.grid_ps_c[:,:,1]
    n_daily_mean[i_day] += 1
end

function Finalize_Output!(output_manager::Output_Manager, save_file_name::String = "None")

    n_day = output_manager.n_day

    u_daily_mean, v_daily_mean, h_daily_mean = output_manager.u_daily_mean, output_manager.v_daily_mean, output_manager.h_daily_mean
    n_daily_mean = output_manager.n_daily_mean
    
    for i_day = 1:n_day
        u_daily_mean[:,:,:,i_day] ./= n_daily_mean[i_day]
        v_daily_mean[:,:,:,i_day] ./= n_daily_mean[i_day]
        h_daily_mean[:,:,1,i_day] ./= n_daily_mean[i_day]
        n_daily_mean[i_day] = 1.0
    end
       
    if save_file_name != "None"
        @save save_file_name u_daily_mean v_daily_mean h_daily_mean
    end

end





function Sigma_Zonal_Mean_Contourf(output_manager::Output_Manager, save_file_pref::String)
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'

    t_zonal_mean, t_eq_zonal_mean, u_zonal_mean, v_zonal_mean = output_manager.t_zonal_mean, output_manager.t_eq_zonal_mean, output_manager.u_zonal_mean, output_manager.v_zonal_mean
   
    PyPlot.contourf(X, Y, t_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_T.png")
    PyPlot.close("all")

    PyPlot.contourf(X, Y, t_eq_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_Teq.png")
    PyPlot.close("all")

    PyPlot.contourf(X, Y, u_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_U.png")
    PyPlot.close("all")

    PyPlot.contourf(X, Y, v_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_V.png")
    PyPlot.close("all")
    
    
end


function Sigma_Zonal_Mean_Pcolormesh(output_manager::Output_Manager, save_file_pref::String)
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'

    t_zonal_mean, u_zonal_mean, v_zonal_mean = output_manager.t_zonal_mean, output_manager.u_zonal_mean, output_manager.v_zonal_mean
   
    PyPlot.pcolormesh(X, Y, t_zonal_mean, shading= "gouraud", cmap="viridis")
    PyPlot.gca().invert_yaxis()
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_T.png")
    PyPlot.close("all")

    PyPlot.pcolormesh(X, Y, u_zonal_mean, shading= "gouraud", cmap="viridis")
    PyPlot.gca().invert_yaxis()
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_U.png")
    PyPlot.close("all")

    PyPlot.pcolormesh(X, Y, v_zonal_mean, shading= "gouraud", cmap="viridis")
    PyPlot.gca().invert_yaxis()
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_V.png")
    PyPlot.close("all")
    
    
end


function Lat_Lon_Pcolormesh(output_manager::Output_Manager, grid_dat::Array{Float64,3}, level::Int64, save_file_name::String = "None")
    
    λc, θc = output_manager.λc, output_manager.θc
    nλ, nθ = length(λc), length(θc)
    λc_deg, θc_deg = λc*180/pi, θc*180/pi
    
    
    X,Y = repeat(λc_deg, 1, nθ), repeat(θc_deg, 1, nλ)'
    
    
    PyPlot.pcolormesh(X, Y, grid_dat[:,:,level], shading= "gouraud", cmap="viridis")
    PyPlot.axis("equal")
    PyPlot.colorbar()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end


function Lat_Lon_Pcolormesh(mesh::Spectral_Spherical_Mesh, grid_dat::Array{Float64,3}, level::Int64, save_file_name::String = "None")
    
    λc, θc = mesh.λc, mesh.θc
    nλ, nθ = length(λc), length(θc)
    λc_deg, θc_deg = λc*180/pi, θc*180/pi
    
    X,Y = repeat(λc_deg, 1, nθ), repeat(θc_deg, 1, nλ)'
    
    
    PyPlot.pcolormesh(X, Y, grid_dat[:,:,level], shading= "gouraud", cmap="viridis")
    PyPlot.axis("equal")
    PyPlot.colorbar()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end


function Zonal_Mean(grid_dat::Array{Float64,3})
    
    return dropdims(mean(grid_dat, dims=1), dims=1)
    
end


function Sigma_Zonal_Mean_Pcolormesh(output_manager::Output_Manager,
    zonal_mean_data::Array{Float64,2}, save_file_name::String = "None")
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'
    
    
    PyPlot.pcolormesh(X, Y, zonal_mean_data, shading= "gouraud", cmap="viridis")
    PyPlot.colorbar()
    PyPlot.gca().invert_yaxis()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end


function Sigma_Zonal_Mean_Contourf(output_manager::Output_Manager, 
    zonal_mean_data::Array{Float64,2}, save_file_name::String = "None")
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'
    
    PyPlot.contourf(X, Y, zonal_mean_data, levels = 10)
    PyPlot.gca().invert_yaxis()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end