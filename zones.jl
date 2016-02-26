# Strip detections based on visual zones (binocular, blind…)

function load_detections(detect_file::AbstractString)
    p, mask = h5read_particles(detect_file)
    detections = h5read(detect_file, "detections")
    p, mask, detections
end

function strip_non_binocular_detections!(detections::Array{Vector{UInt64},3},
                                         p::Matrix{State}, mask::BitMatrix,
                                         grid::Matrix{Vec2},
                                         binocular_angle::Real=deg2rad(31/2))
    @assert(0 ≤ binocular_angle ≤ π)
    strip_detections(detections, p, mask, grid, π - binocular_angle)
end

function strip_blind_detections!(detections::Array{Vector{UInt64},3},
                                 p::Matrix{State}, mask::BitMatrix,
                                 grid::Matrix{Vec2},
                                 blind_angle::Real=deg2rad(25/2))
    @assert(0 ≤ blind_angle ≤ π)
    strip_detections(detections, p, mask, grid, blind_angle)
end

"""
`strip_detections!` removes detections located within an angle
`angle_from_back` from the back of each particle.
"""
function strip_detections!(detections::Array{Vector{UInt64},3},
                           p::Matrix{State}, mask::BitMatrix,
                           grid::Matrix{Vec2},
                           angle_from_back::Real)
    K = size(p, 2)
    nx, ny = size(detections, 1, 2)
    for k in 1:k, j in 1:nx, i in 1:ny
        for l in unpack(detections[i,j,k])
            if abs(angle(grid[i,j] - p[l,k].pos, -p[l,k].vel)) ≤ angle_from_back
                unset(detections[i,j,k], l)
            end
        end
    end
end


# Plotting
# --------

function plot_zones(p::Matrix{State}, mask::BitMatrix,
                    detections_all::Array{Vector{UInt64},3},
                    detections_binoc::Array{Vector{UInt64},3},
                    detections_blind::Array{Vector{UInt64},3})
    k = 1
    nz = findnz(mask[:,k])

    draw(PDF("detect_all", 6inch, 4inch), plotdetect(p[nz,k], detections))
    draw(PDF("detect_binoc", 6inch, 4inch), plotdetect(p[nz,k], detections_binoc))
    draw(PDF("detect_blind", 6inch, 4inch), plotdetect(p[nz,k], detections_blind))
end