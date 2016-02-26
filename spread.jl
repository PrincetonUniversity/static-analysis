# Spread

#export run_spread

"""
`run_spread` computes for each point in the grid the expected cascade size.
"""
function run_spread(detect_file::AbstractString, net_file::AbstractString; α::Real = -0.2)
    p, mask = h5read_particles(detect_file)
    detections = h5read(detect_file, "detections")
    si = h5read(net_file, "social")

    N = size(p, 1) # max swarm size
    K = 2#size(p, 2) # replicates
    nx, ny = size(detections, 1, 2)

    # distance from edge
    dem = zeros(nx * ny, K)
    grid = [Vec2(x,y) for x in linspace(-50, 50, nx), y in linspace(-50, 50, ny)]
    println("Computing distance from edge…")
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = mask[:,k]
        dem[:,k] = dist_from_edge(α, p[nz,k], grid)
    end
    de_max = maximum(dem)
    println("\r100%")

    # distances, order parameters, dataframe
    println("Building DataFrame…")
    cols = [
        (:DistFromEdge,        Float64),
        (:InsideMaxDisk,       Bool),
        (:RelativeDir,         Float64),
        (:Swarm,               Int),
        (:SwarmPolarization,   Float64),
        (:SwarmRotation,       Float64),
        (:SwarmState,          Symbol),
        (:ExpectedCascadeSize, Float64),
    ]
    df = DataFrame([x[2]::DataType for x in cols], [x[1]::Symbol for x in cols], 0)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = find(mask[:,k])

        xmin, xmax = extrema([s.pos.x for s in p[nz,k]])
        ymin, ymax = extrema([s.pos.y for s in p[nz,k]])
        max = min(50 + xmin, 50 - xmax, 50 + ymin, 50 - ymax)

        # order parameters and state, relative direction
        dir = NaN
        op, or = order_parameters(p[nz,k])
        if op > 0.65 && or < 0.35
            state = :Polarized
            dir = relative_dir(p[nz,k], grid)
        elseif op < 0.35 && or > 0.65
            state = :Milling
        elseif op < 0.35 && or < 0.35
            state = :Swarming
        else
            state = :Transitioning
        end

        net = si[nz,nz,k]

        # compute cascade size
        ex = Matrix{Float64}(nx, ny)
        for j in 1:nx, i in 1:ny
            ex[i,j] = expected_cascade_size_brute_force(detections[i,j,k], net, nz)
        end

        # append to dataframe
        row = DataFrame(
            DistFromEdge        = dem[:,k],
            InsideMaxDisk       = dem[:,k] .<= max,
            RelativeDir         = dir,
            Swarm               = k,
            SwarmPolarization   = op,
            SwarmRotation       = or,
            SwarmState          = state,
            ExpectedCascadeSize = ex[:],
        )
        append!(df, row)
    end
    println("\r100%")

    return df
end


function expected_cascade_size_brute_force(detection, net, nz)
    N = 5
    n = 0.0

    I, J, V = findnz(net)
    order = sortperm(V, rev=true)
    edges = zip(I[order], J[order])
    visited = BitMatrix(size(net)...)

    @inbounds for k in 1:N
        r0 = intersect(unpack(detection), nz)
        visited[:] = false
        n += expected_cascade_size_brute_force(r0, visited, net, edges)
    end
    n / N
end

function expected_cascade_size_brute_force(r0::Vector{Int}, visited::BitMatrix, net::Matrix{Float64}, edges)
    oldvis = -1
    newvis = 0
    @inbounds while newvis > oldvis
        oldvis = newvis
        r1 = deepcopy(r0)
        for (i, j) in edges
            if !visited[i,j] && !(i in r0) && (j in r0)
                visited[i,j] = true
                newvis += 1
                if rand() < net[i,j]
                    push!(r1, i)
                end
            end
        end
        r0 = r1
    end
    length(r0)
end

function cart_prod(v::Vector{Float64}, d::Int)
    isempty(v) && return 0.0
    d == 1 && return sum(v)
    v[1] * cart_prod(v[2:end], d - 1) + cart_prod(v[2:end], d)
end


function run_cascade(detections, net, nz)
    nx, ny = size(detections)
    ex = Matrix{Float64}(nx, ny)
    for j in 1:nx, i in 1:ny
        ex[i,j] = expected_cascade_size_brute_force(detections[i,j], net, nz)
    end
    ex
end