# Strip detections based on visual zones (binocular, blind…)

# function load_detections(detect_file::AbstractString)
#     p, mask = h5read_particles(detect_file)
#     detections = h5read(detect_file, "detections")
#     p, mask, detections
# end

function strip_non_binocular_detections!(detections::Array{BitVector,3},
                                         p::Matrix{State}, mask::BitMatrix,
                                         grid::Matrix{Vec2},
                                         binocular_angle::Real=deg2rad(31/2))
    @assert(0 ≤ binocular_angle ≤ π)
    strip_detections!(detections, p, mask, grid, π - binocular_angle)
end

function strip_blind_detections!(detections::Array{BitVector,3},
                                 p::Matrix{State}, mask::BitMatrix,
                                 grid::Matrix{Vec2},
                                 blind_angle::Real=deg2rad(25/2))
    @assert(0 ≤ blind_angle ≤ π)
    strip_detections!(detections, p, mask, grid, blind_angle)
end

"""
`strip_detections!` removes detections located within an angle
`angle_from_back` from the back of each particle.
"""
function strip_detections!(detections::Array{BitVector,3},
                           p::Matrix{State}, mask::BitMatrix,
                           grid::Matrix{Vec2},
                           angle_from_back::Real)
    K = size(p, 2)
    nx, ny = size(detections, 1, 2)
    for k in 1:K, j in 1:nx, i in 1:ny
        for l in find(detections[i,j,k])
            if abs(angle(p[l,k].pos - grid[i,j], p[l,k].vel)) ≤ angle_from_back
                detections[i,j,k][l] = false
            end
        end
    end
end


# Plotting
# --------

function plot_zones(p::Matrix{State}, mask::BitMatrix,
                    detections_all::Array{BitVector, 3},
                    detections_binoc::Array{BitVector, 3},
                    detections_blind::Array{BitVector, 3},
                    grid::Matrix{Vec2})
    k = 1
    nz = find(mask[:,k])

    draw(PDF("detect_all.pdf", 6inch, 4inch),
        plotdetect(p[nz,k], detections_all[:,:,k], grid))
    draw(PDF("detect_binoc.pdf", 6inch, 4inch),
        plotdetect(p[nz,k], detections_binoc[:,:,k], grid))
    draw(PDF("detect_blind.pdf", 6inch, 4inch),
        plotdetect(p[nz,k], detections_blind[:,:,k], grid))
end

function plot_zones_handpicked(p::Matrix{State}, mask::BitMatrix,
                               detections_all::Array{BitVector, 3},
                               detections_binoc::Array{BitVector, 3},
                               detections_blind::Array{BitVector, 3},
                               grid::Matrix{Vec2})
    for (k, name) in enumerate(["polarized", "milling2", "milling", "swarming"])
        nz = find(mask[:,k])
        m = maximum(map(countnz, detections_all[:,:,k]))
        draw(PDF("$(name)_all.pdf", 5inch, 5inch),
            plotdetect(p[nz,k], detections_all[:,:,k], grid, minvalue=0, maxvalue=m, title="full vision — $(name)"))
        draw(PDF("$(name)_binoc.pdf", 5inch, 5inch),
            plotdetect(p[nz,k], detections_binoc[:,:,k], grid, minvalue=0, title="binocular only — $(name)"))
        draw(PDF("$(name)_blind.pdf", 5inch, 5inch),
            plotdetect(p[nz,k], detections_blind[:,:,k], grid, minvalue=0, maxvalue=m, title="with blind zone — $(name)"))
    end
end