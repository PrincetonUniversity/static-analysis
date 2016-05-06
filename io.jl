# Consolidate and reformat datasets for ellipswarm.

@step function get_net_file(p::Project)
    # TODO: run static from here
    path = p.conf[:net_file]
    ishdf5(path) || error("You must run `static` first.")
    path
end

@step function get_detect_file(p::Project)
    # TODO: run detect from here
    path = p.conf[:detect_file]
    ishdf5(path) || error("You must run `detect` first.")
    path
end

@step [get_net_file] function get_pos(p::Project)
    h5read_particles(p.step[get_net_file])
end


# Datasets
# --------

"makeSmallDatasets writes a subset of Brin and Colin's data in a single HDF5 file."
function makeSmallDatasets()
    K = 12004
    N = 101
    p = Matrix{State}(N, K)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        path = @sprintf("/Volumes/Couzin\ Lab/Simon/CTnet/1503/%d.mat", k)
        matopen(path) do file
            # 1 body length = 5 cm | 1 cm = 7.886 px
            X = read(file, "X") ./ (5 * 7.886)
            H = read(file, "H")
            nz = find(X[:,1])
            for i=1:N
                if i <= length(nz)
                    j = nz[i]
                    p[i,k] = State(Vec2(X[j,1], X[j,2]), atan2(H[j,2], H[j,1]))
                else
                    p[i,k] = State(Vec2(0, 0), 0)
                end
            end
        end
    end
    println("\r100%")
    h5write_particles(joinpath(data_path, "somedata.h5"), p)
end

"makeFullDatasets writes Brin and Colin's full data in a single HDF5 file."
function makeFullDatasets()
    root = "/Volumes/Couzin\ Lab/Simon/startle-data/"
    ls = filter(file -> startswith(file, 'T'), readdir(root))
    v = Vector{Matrix{State}}()
    M = collect(reverse(linspace(0, 100, length(ls)+1)))
    pop!(M)
    for dir in ls
        @printf("\r%3d%%", pop!(M))
        file = joinpath(root, dir, "fov_data.clf")
        detected = bitpack(h5read(file, "/fields/detected"))
        # 1 body length = 5 cm | 1 cm = 7.886 px
        x = h5read(file, "/fields/x") ./ (5 * 7.886)
        y = h5read(file, "/fields/y") ./ (5 * 7.886)
        vx = h5read(file, "/fields/heading_x")
        vy = h5read(file, "/fields/heading_y")
        K, N = size(x)
        p = Matrix{State}(N, K)
        for k=1:K
            x[k,:] -= mean(x[k,detected[k,:]])
            y[k,:] -= mean(y[k,detected[k,:]])
            nz = find(detected[k,:])
            for i=1:N
                if i <= length(nz)
                    j = nz[i]
                    p[i,k] = State(Vec2(x[k,j], y[k,j]), Vec2(vx[k,j], vy[k,j]))
                else
                    p[i,k] = State(Vec2(0, 0), Vec2(0, 0))
                end
            end
        end
        push!(v, p[:,[countnz(detected[k,:]) > 0 for k in 1:K]])
    end
    N, K = maximum([size(p, 1) for p in v]), sum([size(p, 2) for p in v])
    p = Matrix{State}(N, K)
    k = 1
    for i in eachindex(v)
        N, K = size(v[i])
        p[1:N,k:k+K-1] = v[i]
        p[N+1:end,k:k+K-1] = State(Vec2(0, 0), Vec2(0, 0))
        k += K
    end
    println("\r100%")
    h5write_particles(joinpath(data_path, "fulldata.h5"), p)
end

function makeHandpickedDatasets()
    p, mask = h5read_particles(joinpath(data_path, "fulldata.h5"))
    # polarized, open milling, closed milling, swarming
    selection = [2988, 2423, 3289, 1299]
    h5write_particles(joinpath(data_path, "handpickeddata.h5"), p[:,selection])
end

function makeUniformRandomDatasets(file::AbstractString; α::Real = -0.2)
    p, mask = h5read_particles(file)

    N = size(p, 1) # max swarm size
    K = size(p, 2) # replicates

    @inbounds for k in 1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = find(mask[:,k])

        # try to conserve polarization
        m = mean(p[nz,k])
        θ = angle(m.vel)
        r = norm(m.vel)
        #σ = acos(min(r, 1))
        σ = pi/2 * acos(min(sqrt(r), 1))
        xmin, xmax, ymin, ymax = bounds(p[nz,k])
        dx, dy = xmax - xmin, ymax - ymin

        # restict position to inside of largest α-shape
        outer, _, _, _ = alphashape(α, p[nz,k])
        poly = Vec2[p[nz[i],k].pos for i in outer[last(findmax(map(length, outer)))]]
        for n in nz
            while true
                pos = Vec2(xmin + dx * rand(), ymin + dy * rand())
                if inpolygon(pos, poly)
                    ϕ = θ + σ * randn()
                    vel = Vec2(cos(ϕ), sin(ϕ))
                    p[n,k] = State(pos, vel)
                    break
                end
            end
        end

    end

    println("\r100%")
    h5write_particles(joinpath(data_path, "fulldata_random.h5"), p)
end


function makeMitchellsBestCandidateRandomDatasets(file::AbstractString; α::Real = -0.2, num_candidates::Integer = 5, σ::Real = deg2rad(10))
    p, mask = h5read_particles(file)

    N = size(p, 1) # max swarm size
    K = size(p, 2) # replicates

    for k in 1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = find(mask[:,k])

        p[nz,k] = mbcRandom(p[nz,k], α=α, num_candidates=num_candidates, σ=σ)
    end

    println("\r100%")
    h5write_particles(joinpath(data_path, "fulldata_random.h5"), p)
end


function allweights(si, mask)
    r = Float64[]
    for k in 1:size(mask, 2)
        nz = mask[:,k]
        append!(r, si[nz,nz,k][:])
    end
    r
end

function mbcRandom(p::Vector{State};
        α::Real = -0.2, num_candidates::Integer = 5, σ::Real = deg2rad(10))
    xmin, xmax, ymin, ymax = bounds(p)
    dx, dy = xmax - xmin, ymax - ymin
    poly = mainshape(alphashape(α, p), [q.pos for q in p])
    q = similar(p)
    for i in eachindex(p)
        dmax = 0.0
        best = Vec2(0, 0)
        for j in 1:num_candidates
            pos = Vec2(xmin + dx * rand(), ymin + dy * rand())
            while !inpolygon(pos, poly)
                pos = Vec2(xmin + dx * rand(), ymin + dy * rand())
            end
            if i == 1
                best = pos
                break
            end
            d = dist(pos, nearest(q[1:i-1], pos).pos)
            if d > dmax
                dmax = d
                best = pos
            end
        end
        q[i] = State(best, nearest(p, best).vel + σ * randn())
    end
    q
end


function nearest(p::Vector{State}, q::Vec2)
    n = length(p)
    n > 0 || return State(q, Vec2(0,0))
    best = State(Vec2(0,0), Vec2(0,0))
    dmin = Inf
    @inbounds for i in eachindex(p)
        d = dist(p[i].pos, q)
        if d < dmin
            dmin = d
            best = p[i]
        end
    end
    best
end


# Bridson’s algorithm
function bridsonRandom(p::Vector{State}; α::Real = -0.2, R::Real = 0.7)
    d = R/√2
    max_candidates = 30

    L = 50
    n = Int(ceil(L / d))
    grid = Matrix{Nullable{Vec2}}(2n, 2n)
    fill!(grid, Nullable{Vec2}())
    idx(v::Vec2) = (n + 1 + Int(floor(v.x / d)), n + 1 + Int(floor(v.y / d)))

    function faraway(p::Vec2, grid)
        i0, j0 = idx(p)
        for j in -2:2, i in -2:2
            q = grid[i0 + i, j0 + j]
            if !isnull(q) && dist(p, get(q)) < R
                return false
            end
        end
        true
    end

    # try to conserve polarization
    m = mean(p)
    r, θ = norm(m.vel), angle(m.vel)
    σ = pi/2 * acos(min(sqrt(r), 1))
    function randvel()
        ϕ = θ + σ * randn()
        Vec2(cos(ϕ), sin(ϕ))
    end

    # restict position to inside of largest α-shape
    poly = mainshape(alphashape(α, p), [q.pos for q in p])

    xmin, xmax, ymin, ymax = bounds(p)
    dx, dy = xmax - xmin, ymax - ymin
    pos = Vec2(xmin + dx * rand(), ymin + dy * rand())
    while !inpolygon(pos, poly)
        pos = Vec2(xmin + dx * rand(), ymin + dy * rand())
    end
    active = [pos]
    grid[idx(pos)...] = pos

    while !isempty(active)
        i = rand(1:length(active))
        current = active[i]
        m = max_candidates
        while m > 0
            ϕ = 2pi * rand()
            r0 = sqrt(R^2 + 3R^2 * rand())
            pos = current + Vec2(r0 * cos(ϕ), r0 * sin(ϕ))
            if inpolygon(pos, poly) && faraway(pos, grid)
                break
            end
            m -= 1
        end
        if m > 0
            push!(active, pos)
            grid[idx(pos)...] = pos
        else
            deleteat!(active, i)
        end
    end

    [State(get(x), randvel()) for x in filter(x->!isnull(x), grid[:])]
end


# Utils
# -----

"state_dtype return the HDF5 compound datatype for a State (must be closed)."
function state_dtype()
    voffset = fieldoffsets(Vec2)
    vt = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(Vec2))
    HDF5.h5t_insert(vt, "X", voffset[1], HDF5.H5T_NATIVE_DOUBLE)
    HDF5.h5t_insert(vt, "Y", voffset[2], HDF5.H5T_NATIVE_DOUBLE)
    soffset = fieldoffsets(State)
    st = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(State))
    HDF5.h5t_insert(st, "Pos", soffset[1], vt)
    HDF5.h5t_insert(st, "Vel", soffset[2], vt)
    HDF5.h5t_close(vt)
    st
end

"h5read_particles reads a matrix of states from the HDF5 file at path."
function h5read_particles(path)
    file = HDF5.h5f_open(path, HDF5.H5F_ACC_RDONLY)
    dset = HDF5.h5d_open(file, "particles")
    fs = HDF5.h5d_get_space(dset)
    dims, maxdims = HDF5.h5s_get_simple_extent_dims(fs)
    HDF5.h5s_close(fs)
    p = Matrix{State}(dims...)
    st = state_dtype()
    HDF5.h5d_read(dset, st, p)
    HDF5.h5t_close(st)
    HDF5.h5d_close(dset)
    HDF5.h5f_close(file)
    mask = bitbroadcast(v -> abs(v.pos.x) > eps(), p)
    p, mask
end

"h5write_particles writes a matrix of states to the HDF5 file at path."
function h5write_particles(path, p::Array{State})
    file = HDF5.h5f_create(path)
    st = state_dtype()
    dims = [convert(HDF5.Hsize, i) for i in reverse(size(p))]
    fs = HDF5.h5s_create_simple(length(dims), dims, dims)
    dset = HDF5.h5d_create(file, "particles", st, fs)
    HDF5.h5s_close(fs)
    HDF5.h5d_write(dset, st, p)
    HDF5.h5t_close(st)
    HDF5.h5d_close(dset)
    HDF5.h5f_close(file)
    nothing
end
