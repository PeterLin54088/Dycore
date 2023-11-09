export Thermal_Wind_from_Height!

function Thermal_Wind_from_Height!(mesh::Spectral_Spherical_Mesh, height::Array{Float64,3}, grid_u::Array{Float64,3}, grid_v::Array{Float64,3})
    size = 10
    instantiate_matrix = 1.0* Matrix(I, size, size)
    println(instantiate_matrix)
    return nothing
end