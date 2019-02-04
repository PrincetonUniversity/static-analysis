# Predator model

function predator_project()
    Project("predator",
        α = -0.2,
        colormap = :plasma,
        theme = Dict(
            :minor_label_font => "CMU Serif",
            :minor_label_font_size => 10pt,
            :major_label_font => "CMU Serif",
            :major_label_font_size => 12pt,
            :key_title_font => "CMU Serif",
            :key_title_font_size => 12pt,
            :key_label_font => "CMU Serif",
            :key_label_font_size => 11pt,
            :plot_padding => 8pt,
        ),
    )
end

immutable Domain
    size::Float64
end

inside(p::Vec2, d::Domain) = 0 <= p.x < d.size && 0 <= p.y < d.size
wrap(p::Vec2, d::Domain) = mod(p, d.size)
unitvec(u::Vec2, v::Vec2, d::Domain) = vec(u, v, d) / dist(u, v, d)
dist(u::Vec2, v::Vec2, d::Domain) = norm(vec(u, v, d))

unitvec(u::Vec2) = u / norm(u)

function vec(u::Vec2, v::Vec2, d::Domain)
    x, y = v.x - u.x, v.y - u.y
    if 2x <= -d.size
        x += d.size
    elseif 2x > d.size
        x -= d.size
    end
    if 2y <= -d.size
        y += d.size
    elseif 2y > d.size
        y -= d.size
    end
    Vec2(x, y)
end

type Params
	N::Int
	m::Float64
	α::Float64
	β::Float64
	Cr::Float64
	lr::Float64
	Ca::Float64
	la::Float64
end


function grad_morse_potential(p::Vector{State}, u::State, Cr, lr, Ca, la, dom::Domain, skip=-1)
    w = Vec2(0, 0)
    @inbounds for i in eachindex(p)
        i == skip && continue
        v = p[i]
        d = dist(u.pos, v.pos, dom)
        U = (Cr * exp(-d / lr) / lr - Ca * exp(-d / la) / la) / d
        w += U * vec(u.pos, v.pos, dom)
    end
    w
end

function RK4(f, x, dt)
    k1 = f(x)
    k2 = f(x + dt / 2 * k1)
    k3 = f(x + dt / 2 * k2)
    k4 = f(x + dt * k3)
    x + dt * (k1 + 2k2 + 2k3 + k4) / 6
end

# Unused
function update!(p, q, m, α, β, Cr, lr, Ca, la, dt, dom::Domain)
    @inbounds for i in eachindex(p)
        u = p[i]
        ∇U = grad_morse_potential(p, u, Cr, lr, Ca, la, dom, i)
        f = v -> ((α - β * norm(v)^2) * v - ∇U) / m
        q[i] = State(wrap(u.pos + dt * u.vel, dom), RK4(f, u.vel, dt))
    end
end

function update_pred!(pred, pred_new, prey, pd::Params, dt, dom::Domain, follow)
	@inbounds for i in eachindex(pred)
		u = pred[i]
		if follow[i] == i
			∇U = grad_morse_potential(prey, u, 0, 1, pd.Ca, pd.la, dom)
		else
			j = follow[i]
			∇U = grad_morse_potential(pred[j:j], u, 4, 1, 3, 3, dom)
		end
		f = v -> ((pd.α - pd.β * norm(v)^2) * v - ∇U) / pd.m
		pred_new[i] = State(wrap(u.pos + dt * u.vel, dom), RK4(f, u.vel, dt))
	end
end

function update_prey!(prey, prey_new, pred, py::Params, dt, dom::Domain)
	@inbounds for i in eachindex(prey)
		u = prey[i]
		∇U = grad_morse_potential(prey, u, py.Cr, py.lr, py.Ca, py.la, dom, i)
		∇V = grad_morse_potential(pred, u, 10py.Cr, py.la, 0, 1, dom)
		f = v -> ((py.α - py.β * norm(v)^2) * v - ∇U - ∇V) / py.m
		prey_new[i] = State(wrap(u.pos + dt * u.vel, dom), RK4(f, u.vel, dt))
	end
end

function init!(p::Vector{State}, v0, dom::Domain)
    @inbounds for i in eachindex(p)
        θ = 2pi * rand()
        p[i] = State(dom.size * Vec2(rand(), rand()), v0 * Vec2(cos(θ), sin(θ)))
    end
end

function init!(pd::Params, py::Params, pred, prey, follow, p, dom::Domain)
    # Pred
    v0 = sqrt(pd.α / pd.β)
	@inbounds for i in eachindex(pred)
		θ = 2pi * rand()
		pred[i] = State(dom.size * Vec2(rand(), rand()), v0 * Vec2(cos(θ), sin(θ)))
		follow[i] = rand() < p ? i : 0
	end
	if !isempty(pred) && all(follow .== 0)
		follow[1] = 1
	end
	for i in randperm(size(pred, 1))
		if follow[i] == 0
			decided = find(follow)
			@assert !isempty(decided)
			j = decided[rand(1:length(decided))]
			while follow[j] != j
				j = follow[j]
			end
			follow[j] = i
			follow[i] = i
		end
	end

    # Prey
    v0 = sqrt(py.α / py.β)
	@inbounds for i in eachindex(prey)
		θ = 2pi * rand()
		prey[i] = State(dom.size * Vec2(rand(), rand()), v0 * Vec2(cos(θ), sin(θ)))
	end
end

function alignment(p::Vector{State})
    v = Vec2(0, 0)
    @inbounds for i in eachindex(p)
        v += unitvec(p[i].vel)
    end
    norm(v) / length(p)
end

function pplot(p::Vector{State}, dom::Domain)
    box = UnitBox(0, 0, dom.size, dom.size)
    λ = dom.size/400 # velocity scale factor (velocity shown as tail length)
    swarm = compose(context(units=box),
        (context(), fill(colorant"black"), circle([m.pos.x for m in p], [m.pos.y for m in p], [2px])),
        (context(), stroke(colorant"black"), line([[(m.pos.x, m.pos.y), (m.pos.x-λ*m.vel.x, m.pos.y-λ*m.vel.y)] for m in p]))
    )
    compose(context(), swarm, rectangle(), fill(colorant"snow"))
end

# Sample params: N=200, m=1, α=1.6, β=0.5, Cr=2, lr=0.5, Ca=0.5, la=2

pred = Params(10, 2, 0.8, 0.2, 0, 1, 3, 3)
prey = Params(200, 1, 1.2, 0.37, 2, 1, 1, 4)

function sim(pred::Params, prey::Params; steps=1000, dt=0.1, dom=Domain(100), p=0.5, α=-0.2, d=5.0, video=false, verbose=true)
    pred_old = Vector{State}(pred.N)
    prey_old = Vector{State}(prey.N)
    pred_new = Vector{State}(pred.N)
    prey_new = Vector{State}(prey.N)

    pred_risk = zeros(pred.N, steps)
    prey_risk = zeros(prey.N, steps)

    follow = Vector{Int}(pred.N)

    init!(pred, prey, pred_old, prey_old, follow, p, dom)

    vidpath = "/Users/sleblanc/data/swarm-parameter-scan/vid"
    video && try rm(vidpath, recursive=true) end
    video && mkdir(vidpath)

    for k=1:steps
        verbose && @printf("\r% 5.1f%%", 100(k - 1) / steps)
        video && draw(PNG(@sprintf("%s/vid%08d.png", vidpath, k), 6inch, 6inch), pplot(prey_old, dom))

        risk!(k, pred_risk, prey_risk, pred_old, prey_old, α, d, dom)

        update_pred!(pred_old, pred_new, prey_old, pred, dt, dom, follow)
        update_prey!(prey_old, prey_new, pred_old, prey, dt, dom)
        pred_old, pred_new = pred_new, pred_old
        prey_old, prey_new = prey_new, prey_old
    end

    # video && run(pipeline(`ffmpeg -y -i vid/vid%08d.png -an -pix_fmt yuv420p -vcodec h264 -q 1 -r 24 vid.mp4`, stdout=DevNull, stderr=DevNull))
    verbose && println("\r100.0%")

    pred_risk, prey_risk, follow
end


function runsims()
    N = 5
    p = 0.1:0.1:0.9
    riskpd = zeros(N, length(p))
    riskpy = zeros(N, length(p))
    chains = zeros(10, length(p))
    for i in eachindex(p)
        println("# p = ", p[i])
        for k in 1:N
            rpd, rpy, follow = sim(pred, prey, p=p[i])
            riskpd[k,i] = mean(rpd)
            riskpy[k,i] = mean(rpy)
            chains[:,i] += chainstats(follow)
            println("  - risk pred = ", riskpd[k,i])
            println("  - risk prey = ", riskpy[k,i])
        end
        chains[:,i] /= N
    end
    riskpd, riskpy, chains
end


function risk!(k, pred_risk, prey_risk, pred, prey, α, d, dom::Domain)
    shape = alphashape(α, prey)
    nschool_per_pred = zeros(Int, length(pred))
    @inbounds for n in eachindex(shape.outer)
        poly = Vec2[prey[i].pos for i in shape.outer[n]]
        idx = Bool[inpolygon(p.pos, poly) for p in prey]
        A = area(Vec2[prey[i].pos for i in shape.outer[n]])
        P = falses(length(pred))
        for j in eachindex(pred)
            for i in find(idx)
                if dist(prey[i].pos, pred[j].pos, dom) < d
                    P[j] = true
                    break
                end
            end
        end
        if any(P)
            N = countnz(P)
            nschool_per_pred[P] += 1
            pred_risk[P,k] += N / A
            prey_risk[idx,k] = N / A
        end
    end
    @inbounds for i in eachindex(pred)
        if nschool_per_pred[i] > 0
            pred_risk[i,k] = pred_risk[i,k] / nschool_per_pred[i]
        end
    end

    # Degen (TODO) and solo (done)
    # @inbounds for i in eachindex(shape.solo)
    #     for j in eachindex(pred)
    #         if dist(prey[i].pos, pred[j].pos, dom) < d
    #             pred_risk[j,k] = Inf
    #             prey_risk[i,k] = Inf
    #             break
    #         end
    #     end
    # end
end


function chainstats(follow)
    depth(follow, i) = follow[i] == i ? 1 : 1 + depth(follow, follow[i])
    len = zeros(Int, length(follow))
    for i in eachindex(follow)
        len[i] = depth(follow, i)
    end
    hist = zeros(Int, length(follow))
    for i in eachindex(len)
        hist[len[i]] += 1
        if len[i] > 1
            hist[len[i]-1] -= 1
        end
    end
    hist
end


# function area(coords)
#     b = coords[end]
#     v = 0.0
#     @inbounds for i in eachindex(coords)
#         a, b = b, coords[i]
#         v += a.y * b.x - a.x * b.y
#     end
#     v / 2
# end

function area(poly::Vector{Vec2})
    xmin, xmax = extrema(Float64[p.x for p in poly])
    ymin, ymax = extrema(Float64[p.y for p in poly])
    n = 1000
    k = @parallel (+) for i in 1:n
        p = Vec2(xmin + rand() * (xmax - xmin), ymin + rand() * (ymax - ymin))
        Int(inpolygon(p, poly))
    end
    k / n
end

# id = groups(X, r) categorizes each of n points into groups of maximum
# interindividual distance r. X is an (n x 2) matrix containing the x and y
# coordinates of the n points. id is a vector of length n containing the
# group id of each point. Solitary individuals all get id 0.
function groups(s::Vector{State}, r, dom::Domain)
    @assert(r > 0, "r must be a positive scalar")

    id = zeros(Int, length(s))
    nid = 1
    @inbounds for i=1:length(s)
        for j=i+1:length(s)
            d = dist(s[i].p, s[j].p, dom)
            if d <= r
                if id[i] > 0 && id[j] > 0 && id[i] != id[j]
                    # merge i's group and j's group
                    x = min(id[i], id[j])
                    y = max(id[i], id[j])
                    id[id .== y] = x
                    id[id .> y] -= 1
                    nid -= 1
                elseif id[i] > 0
                    # attach j to i's group
                    id[j] = id[i]
                elseif id[j] > 0
                    # attach i to j's group
                    id[i] = id[j]
                else
                    # form a new group
                    id[i] = nid
                    id[j] = nid
                    nid = nid + 1
                end
            end
        end
    end
    id
end

function groups(s::Matrix{State}, r, dom::Domain)
    hist = zeros(size(s))
    for k=1:size(s, 2)
        id = groups(s[:,k], r, dom)
        for i=unique(id), j=1:size(s, 1)
            c = count(x -> x == i, id)
            if c > 0
                hist[c,k] += 1
            end
        end
    end
    hist ./ size(s, 2)
end

function plotgroups(pred::Vector{State}, prey::Vector{State}, follow, r, dom::Domain)
    box = UnitBox(0, 0, dom.size, dom.size)
    λ = 0.4 # velocity scale factor (velocity shown as tail length)

    # Pred
    preds = context(units=box)
    for i=1:length(pred)
        m = pred[i]
        n = 1
        while follow[i] != i
            i = follow[i]
            n += 1
        end
        col = Color.HSV(0, 1, 0.9^n)
        preds = compose(preds,
            (context(), fill(col), circle(m.p.x, m.p.y, 4px)),
            (context(), stroke(col), linewidth(3px), line([(m.p.x, m.p.y), (m.p.x-1.5λ*m.v.x, m.p.y-1.5λ*m.v.y)]))
        )
    end

    # Prey
    id = groups(prey, r, dom)
    preys = context(units=box)
    n = maximum(id)
    for i in unique(id)
        p = prey[id .== i]
        col = Color.HSV(360i/n, 1, 0.7)
        preys = compose(preys,
            (context(), fill(col), circle([m.p.x for m in p], [m.p.y for m in p], [2px])),
            (context(), stroke(col), line([[(m.p.x, m.p.y), (m.p.x-λ*m.v.x, m.p.y-λ*m.v.y)] for m in p]))
        )
    end
    compose(context(), preds, preys, rectangle(), fill(colorant"snow"))
end

