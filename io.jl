# IO related functions

using HDF5
using MAT

"state_dtype return the HDF5 compound datatype for a State (must be closed)."
function state_dtype()
    vt = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(Vec2))
    HDF5.h5t_insert(vt, "X", 0, HDF5.H5T_NATIVE_DOUBLE)
    HDF5.h5t_insert(vt, "Y", 8, HDF5.H5T_NATIVE_DOUBLE)
    st = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(State))
    HDF5.h5t_insert(st, "Pos", 0, vt)
    HDF5.h5t_insert(st, "Vel", 16, vt)
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
    h5write_particles("data/colin2.h5", p)
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
    h5write_particles("data/fulldata.h5", p)
end

