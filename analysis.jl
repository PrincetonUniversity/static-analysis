# Analyze ellipswarm static data

using HDF5
using DataFrames
using Gadfly
using Compose
using Colors
using MAT

"Vec2 is a 2D vector."
immutable Vec2
	x::Float64
	y::Float64
end

"A State contains the position and orientation of a particle."
immutable State
	pos::Vec2
	dir::Float64
end

import Base: +, -, *, /
for op in [:+, :-, :*, :/]
	@eval $op(u::Vec2, v::Vec2) = Vec2($op(u.x, v.x), $op(u.y, v.y))
	@eval $op(x::Real, u::Vec2) = Vec2($op(x, u.x), $op(x, u.y))
	@eval $op(u::Vec2, x::Real) = Vec2($op(u.x, x), $op(u.y, x))
end

import Base: norm, dot
norm(u::Vec2) = hypot(u.x, u.y)
dot(u::Vec2, v::Vec2) = u.x * v.x + u.y * v.y

import Base: mean
function mean(p::Vector{State})
	n = length(p)
	mx, my, dx, dy = 0.0, 0.0, 0.0, 0.0
	for i in 1:n
		mx += p[i].pos.x
		my += p[i].pos.y
		dx += cos(p[i].dir)
		dy += sin(p[i].dir)
	end
	State(Vec2(mx / n, my / n), atan2(dy, dx))
end


"run makes a plot of personal information and social influence for groups of particles."
function run(file::AbstractString, α::Real = -0.5)

	p = h5read_particles(file)
	pi = h5read(file, "personal")
	si = h5read(file, "social")

	N = size(p, 1) # swarm size
	K = size(p, 2) # replicates

	# distance from back
	println("Computing distance from back…")
	db = zeros(N, K)
	for k=1:K
		n = countnz([abs(v.pos.x) > eps() for v in p[:,k]])
		db[1:n,k] = dist_from_back(p[1:n,k])
	end


	# distance from edge
	println("Computing distance from edge…")
	de = zeros(N, K)
	for k=1:K
		@printf("\r%3d%%", 100(k-1) / K)
		n = countnz([abs(v.pos.x) > eps() for v in p[:,k]])
		de[1:n,k] = dist_from_edge(α, p[1:n,k])
	end
	println("\r100%")

	# dataframe
	println("Building DataFrame…")
	df = DataFrame([Float64, Float64, Float64, Symbol], [:DistFromBack, :DistFromEdge, :Info, :Kind], 0)
	for k=1:K
		n = countnz([abs(v.pos.x) > eps() for v in p[:,k]])
		w = si[1:n,1:n,k]
		msi = cc(w .* (w .> 1e-4)) # FIXME
		append!(df, DataFrame(DistFromBack=db[1:n,k], DistFromEdge=de[1:n,k], Info=pi[1:n,k], Kind=:Personal))
		append!(df, DataFrame(DistFromBack=db[1:n,k], DistFromEdge=de[1:n,k], Info=msi, Kind=:Social))
	end

	# bins
	nbins = 11
	db_edges = linspace(0, 1, nbins + 1)
	de_edges = linspace(0, nextfloat(maximum(de)), nbins + 1)
	db_centers = (db_edges[1:end-1] + db_edges[2:end]) / 2
	de_centers = (de_edges[1:end-1] + de_edges[2:end]) / 2
	db_pe_bins = [Float64[] for _ in 1:nbins]
	db_so_bins = [Float64[] for _ in 1:nbins]
	de_pe_bins = [Float64[] for _ in 1:nbins]
	de_so_bins = [Float64[] for _ in 1:nbins]
	for i=1:size(df, 1)
		d = df[i,:DistFromBack]
		for j=1:nbins
			if db_edges[j] <= d < db_edges[j+1]
				if df[i,:Kind] == :Personal
					push!(db_pe_bins[j], df[i,:Info])
				else
					push!(db_so_bins[j], df[i,:Info])
				end
				break
			end
		end
		d = df[i,:DistFromEdge]
		for j=1:nbins
			if de_edges[j] <= d < de_edges[j+1]
				if df[i,:Kind] == :Personal
					push!(de_pe_bins[j], df[i,:Info])
				else
					push!(de_so_bins[j], df[i,:Info])
				end
				break
			end
		end
	end
	m_db_pe_bins = map(mean, db_pe_bins)
	m_db_so_bins = map(mean, db_so_bins)
	m_de_pe_bins = map(mean, de_pe_bins)
	m_de_so_bins = map(mean, de_so_bins)
	s_db_pe_bins = map(x -> std(x) / sqrt(length(x)), db_pe_bins)
	s_db_so_bins = map(x -> std(x) / sqrt(length(x)), db_so_bins)
	s_de_pe_bins = map(x -> std(x) / sqrt(length(x)), de_pe_bins)
	s_de_so_bins = map(x -> std(x) / sqrt(length(x)), de_so_bins)
	dfb = DataFrame(DistFromBack=db_centers, MeanInfo=m_db_pe_bins, LowInfoSE=m_db_pe_bins-s_db_pe_bins, HighInfoSE=m_db_pe_bins+s_db_pe_bins, Kind=:Personal)
	append!(dfb, DataFrame(DistFromBack=db_centers, MeanInfo=m_db_so_bins, LowInfoSE=m_db_so_bins-s_db_so_bins, HighInfoSE=m_db_so_bins+s_db_so_bins, Kind=:Social))
	dfe = DataFrame(DistFromEdge=de_centers, MeanInfo=m_de_pe_bins, LowInfoSE=m_de_pe_bins-s_de_pe_bins, HighInfoSE=m_de_pe_bins+s_de_pe_bins, Kind=:Personal)
	append!(dfe, DataFrame(DistFromEdge=de_centers, MeanInfo=m_de_so_bins, LowInfoSE=m_de_so_bins-s_de_so_bins, HighInfoSE=m_de_so_bins+s_de_so_bins, Kind=:Social))

	# normalize units
	# for kind in [:Personal, :Social]
	# 	m, s = mean_and_std(df[df[:Kind].==kind,:Info])
	# 	df[df[:Kind].==kind,:Info] -= m
	# 	df[df[:Kind].==kind,:Info] /= s
	# end
	for kind in [:Personal, :Social]
		id = dfb[:Kind].==kind
		m, s = mean_and_std(dfb[id,:MeanInfo])
		dfb[id,:MeanInfo] -= m
		dfb[id,:MeanInfo] /= s
		dfb[id,:LowInfoSE] -= m
		dfb[id,:LowInfoSE] /= s
		dfb[id,:HighInfoSE] -= m
		dfb[id,:HighInfoSE] /= s
		id = dfe[:Kind].==kind
		m, s = mean_and_std(dfe[id & !isnan(dfe[:,:MeanInfo]),:MeanInfo])
		dfe[id,:MeanInfo] -= m
		dfe[id,:MeanInfo] /= s
		dfe[id,:LowInfoSE] -= m
		dfe[id,:LowInfoSE] /= s
		dfe[id,:HighInfoSE] -= m
		dfe[id,:HighInfoSE] /= s
	end

	return nothing, df, dfb, dfe

	# plot
	println("Building plot…")
	# p1 = plot(df, x=:DistFromBack, y=:Info, color=:Kind, Geom.smooth,
	# 	Guide.xlabel("Normalized distance from back"),
	# 	Guide.ylabel("Standardized units"))
	# p2 = plot(df, x=:DistFromEdge, y=:Info, color=:Kind, Geom.smooth,
	# 	Guide.xlabel("Distance from edge (body length)"),
	# 	Guide.ylabel("Standardized units"))
	p1 = plot(dfb, x=:DistFromBack, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Kind, Geom.point, Geom.errorbar,
		Guide.xlabel("Normalized distance from back"),
		Guide.ylabel("Standardized units"))
	p2 = plot(dfe, x=:DistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Kind, Geom.point, Geom.errorbar,
		Guide.xlabel("Distance from edge (body length)"),
		Guide.ylabel("Standardized units"))
	p = vstack(p1, p2)

	draw(PDF("info.pdf", 6inch, 8inch), p)

	p, df, dfb, dfe
end

function runscale(mode)
	println("> scale/data_0.5_$(mode).h5")
	_, _, dfb, dfe = run("scale/data_0.5_$(mode).h5")
	dfb[:Scale] = 0.5
	dfe[:Scale] = 0.5
	for scale in logspace(-log10(2), log10(2), 9)[2:end]
		println("> scale/data_$(scale)_$(mode).h5")
		_, _, dfb2, dfe2 = run("scale/data_$(scale)_$(mode).h5", -0.5 / scale)
		dfb2[:Scale] = scale
		dfe2[:Scale] = scale
		append!(dfb, dfb2)
		append!(dfe, dfe2)
	end
	writetable("scale/dfb_$(mode).csv", dfb)
	writetable("scale/dfe_$(mode).csv", dfe)
	dfb, dfe
end

function runsize(mode)
	println("> size/data_5.0_$(mode).h5")
	_, _, dfb, dfe = run("size/data_5.0_$(mode).h5")
	dfb[:Size] = 5.0
	dfe[:Size] = 5.0
	for size in 10.0:5.0:25.0
		println("> size/data_$(size)_$(mode).h5")
		_, _, dfb2, dfe2 = run("size/data_$(size)_$(mode).h5", -0.5)
		dfb2[:Size] = size
		dfe2[:Size] = size
		append!(dfb, dfb2)
		append!(dfe, dfe2)
	end
	writetable("size/dfb_$(mode).csv", dfb)
	writetable("size/dfe_$(mode).csv", dfe)
	dfb, dfe
end

function plotscale(mode)
	dfb = readtable("scale/dfb_$mode.csv")
	dfb[:Kind] = convert(DataArrays.DataArray{Symbol,1}, dfb[:Kind])
	dfe = readtable("scale/dfe_$mode.csv")
	dfe[:Kind] = convert(DataArrays.DataArray{Symbol,1}, dfe[:Kind])

	x0, x1 = extrema((linspace(0, 1, 12)[1:end-1] + linspace(0, 1, 12)[2:end])/2)
	dfe[:NormDistFromEdge] = 0.0
	for scale in unique(dfe[:Scale])
		id = dfe[:Scale].==scale
		min, max = extrema(dfe[id,:DistFromEdge])
		dfe[id,:NormDistFromEdge] = x0 + (dfe[id,:DistFromEdge] - min) * (x1 - x0) / (max - min)
	end

	p1 = plot(dfb[dfb[:Kind].==:Personal,:], x=:DistFromBack, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Normalized distance from back"),
		Guide.ylabel("Standardized units"),
		Guide.title("External visual field ($mode)"),
		Scale.color_continuous(minvalue=0.5, maxvalue=2.0))
	p2 = plot(dfe[dfe[:Kind].==:Personal,:], x=:DistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Distance from edge (body length)"),
		Guide.ylabel("Standardized units"),
		Guide.title("External visual field ($mode)"),
		Scale.color_continuous(minvalue=0.5, maxvalue=2.0))
	p3 = plot(dfe[dfe[:Kind].==:Personal,:], x=:NormDistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Normalized distance from edge"),
		Guide.ylabel("Standardized units"),
		Guide.title("External visual field ($mode)"),
		Scale.color_continuous(minvalue=0.5, maxvalue=2.0))
	p = vstack(p1, p2, p3)

	draw(PDF("scale_personal_$(mode).pdf", 6inch, 12inch), p)

	p1 = plot(dfb[dfb[:Kind].==:Social,:], x=:DistFromBack, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Normalized distance from back"),
		Guide.ylabel("Standardized units"),
		Guide.title("Social influence ($mode)"),
		Scale.color_continuous(minvalue=0.5, maxvalue=2.0))
	p2 = plot(dfe[dfe[:Kind].==:Social,:], x=:DistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Distance from edge (body length)"),
		Guide.ylabel("Standardized units"),
		Guide.title("Social influence ($mode)"),
		Scale.color_continuous(minvalue=0.5, maxvalue=2.0))
	p3 = plot(dfe[dfe[:Kind].==:Social,:], x=:NormDistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Normalized distance from edge"),
		Guide.ylabel("Standardized units"),
		Guide.title("Social influence ($mode)"),
		Scale.color_continuous(minvalue=0.5, maxvalue=2.0))
	p = vstack(p1, p2, p3)

	draw(PDF("scale_social_$(mode).pdf", 6inch, 12inch), p)
end

function plotsize(mode)
	dfb = readtable("size/dfb_$mode.csv")
	dfb[:Kind] = convert(DataArrays.DataArray{Symbol,1}, dfb[:Kind])
	dfe = readtable("size/dfe_$mode.csv")
	dfe[:Kind] = convert(DataArrays.DataArray{Symbol,1}, dfe[:Kind])

	x0, x1 = extrema((linspace(0, 1, 12)[1:end-1] + linspace(0, 1, 12)[2:end])/2)
	dfe[:NormDistFromEdge] = 0.0
	for size in unique(dfe[:Size])
		id = dfe[:Size].==size
		min, max = extrema(dfe[id,:DistFromEdge])
		dfe[id,:NormDistFromEdge] = x0 + (dfe[id,:DistFromEdge] - min) * (x1 - x0) / (max - min)
	end

	p1 = plot(dfb[dfb[:Kind].==:Personal,:], x=:DistFromBack, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Normalized distance from back"),
		Guide.ylabel("Standardized units"),
		Guide.title("External visual field ($mode)"),
		Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
	p2 = plot(dfe[dfe[:Kind].==:Personal,:], x=:DistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Distance from edge (body length)"),
		Guide.ylabel("Standardized units"),
		Guide.title("External visual field ($mode)"),
		Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
	p3 = plot(dfe[dfe[:Kind].==:Personal,:], x=:NormDistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Normalized distance from edge"),
		Guide.ylabel("Standardized units"),
		Guide.title("External visual field ($mode)"),
		Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
	p = vstack(p1, p2, p3)

	draw(PDF("size_personal_$(mode).pdf", 6inch, 12inch), p)

	p1 = plot(dfb[dfb[:Kind].==:Social,:], x=:DistFromBack, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Normalized distance from back"),
		Guide.ylabel("Standardized units"),
		Guide.title("Social influence ($mode)"),
		Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
	p2 = plot(dfe[dfe[:Kind].==:Social,:], x=:DistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Distance from edge (body length)"),
		Guide.ylabel("Standardized units"),
		Guide.title("Social influence ($mode)"),
		Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
	p3 = plot(dfe[dfe[:Kind].==:Social,:], x=:NormDistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
		Guide.xlabel("Normalized distance from edge"),
		Guide.ylabel("Standardized units"),
		Guide.title("Social influence ($mode)"),
		Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
	p = vstack(p1, p2, p3)

	draw(PDF("size_social_$(mode).pdf", 6inch, 12inch), p)
end

"cc computes the local weighted directed clustering coefficient."
function cc(w::Matrix{Float64})
	# assuming w is square weight matrix
	n = size(w, 1)
	a = w .!= 0 # adjacency matrix
	CC = zeros(n)
	for i=1:n
		dg, db, C = 0.0, 0.0, 0.0
		for j=1:n
			for k=1:n
				C += (w[i,j] + w[j,i]) * (w[i,j] + w[k,i]) * (w[j,k] + w[k,j])
			end
			j != i || continue
			dg += a[i,j] + a[j,i]
			db += a[i,j] * a[j,i]
		end
		CC[i] = C == 0 ? 0 : C / 2(dg * (dg - 1) - 2db)
	end
	CC
end

"dist_from_back computes the distance of each particle to the back of the swarm."
function dist_from_back(p::Vector{State})
	db = similar(p, Float64)
	m = mean(p)
	for i in eachindex(p)
		@inbounds db[i] = dot(p[i].pos - m.pos, Vec2(cos(p[i].dir), sin(p[i].dir)))
	end
	min, max = extrema(db)
	(db - min) ./ (max - min)
end

"dist_from_edge computes the distance of each particle to the edge of the α-shape."
function dist_from_edge(α::Real, p::Vector{State})
	de = similar(p, Float64)
	outer, inner, degen, solo = alphashape(α, p)
	@inbounds for i in eachindex(p)
		if i in solo
			de[i] = 0
			continue
		end
		dmin = Inf
		m = p[i].pos
		for kind in (outer, inner, degen), v in kind, j=2:length(v)
			a, b = p[v[j-1]].pos, p[v[j]].pos
			u = b - a
			t = dot(m - a, u) / norm(u)
			d = norm(a + clamp(t, 0, 1) * u - m)
			dmin = min(dmin, d)
		end
		de[i] = dmin
	end
	de
end


"netplot plots a weighted directed network of oriented ellipses."
function netplot(p::Vector{State}, w::Matrix{Float64})
	xmin, xmax = extrema([v.pos.x for v in p])
	ymin, ymax = extrema([v.pos.y for v in p])
	box = UnitBox(xmin - 1, ymin - 1, xmax - xmin + 2, ymax - ymin + 2)
	c = context(units=box)
	p = []
	for i in eachindex(p)
		ctx = context(units=box, rotation=Rotation(p[i].dir, (p[i].pos.x, p[i].pos.y)))
		push!(p, (ctx, ellipse(p[i].pos.x - 0.8/2, p[i].pos.y, 1/2, 0.125/2)))
	end
	#c = netplot(x, y, w)
	e = []
	for j in eachindex(x), i in eachindex(y)
		push!(e, (c, line([(p[i].pos.x, p[i].pos.y), (p[j].pos.x, p[j].pos.y)]), stroke(RGBA(i<j, 0, i>j, w[i,j]))))
	end
	c = compose(c, e...)
	compose(c, p...)
end

"makeDatasets writes Brin and Colin's data in a single HDF5 file."
function makeDatasets()
	K = 12004
	N = 101
	p = Matrix{State}(N, K)
	for k=1:K
		@printf("\r%3d%%", 100(k-1) / K)
		path = @sprintf("/Volumes/Couzin\ Lab/CTnet/1503/%d.mat", k)
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
	h5write_particles("colin2.h5", p)
end

"alphashape computes the α-shape of a group of particles."
function alphashape(α::Real, p::Vector{State})
	n = length(p)
	r = abs(1/α)
	if α < 0
		incircle(x, xc, vc) = norm(x - xc) ≥ r
	elseif α == 0
		incircle(x, xc, vc) = dot(x - xc, vc) ≤ 0
	else
		incircle(x, xc, vc) = norm(x - xc) ≤ r
	end

	# find list of edges that are part of α-shape
	edges = Tuple{Int,Int}[]
	degen = Vector{Int}[]
	dist = zeros(n, n)
	@inbounds for i=1:n, j=i+1:n
		u = p[j].pos - p[i].pos
		d = norm(u)
		dist[i,j] = d

		# skip points that are too far away
		d ≤ 2r || continue

		# find center of circles
		h = sqrt(r^2 - (d/2)^2)
		vh = p[i].pos + u/2
		v = Vec2(-u.y, u.x) / d
		v0 = vh + h * v
		v1 = vh - h * v

		# skip edge if circle does not contain all points
		c0, c1 = true, true
		for k=1:n
			i != k != j || continue
			c0 &= incircle(p[k].pos, v0, v)
			c1 &= incircle(p[k].pos, v1, v)
			c1 || c0 || break
		end
		c1 && c0 && push!(degen, [i, j])
		c1 $ c0 && push!(edges, (i, j))
	end

	# find solitary particles
	solo = Int[]
	@inbounds for i=1:n
		if all(dist[i,i+1:end] .≥ 2r) && all(dist[1:i-1,i] .≥ 2r)
			push!(solo, i)
		end
	end

	# build list of points from list of edges
	outer = Vector{Int}[]
	@inbounds while !isempty(edges)
		e = pop!(edges)
		out = Int[e[1]]
		ix, vx = 1, e[2]
		while vx != out[1]
			push!(out, vx)
			i = findfirst(e -> vx in e, edges)
			e = edges[i]
			vx = e[1] == vx ? e[2] : e[1]
			deleteat!(edges, i)
		end
		push!(outer, out)
	end

	# distinguish inner and outer α-shapes
	m = length(outer)
	m > 0 || return outer, outer, degen, solo
	q = falses(m, m)
	@inbounds for i=1:m, j=1:m
		i != j || continue
		q[i,j] = inpolygon(p[outer[i][1]].pos, [v.pos for v in p[outer[j]]])
	end
	out = map(iseven, squeeze(sum(q, 2), 2))
	inner = outer[!out]
	outer = outer[out]
	outer, inner, degen, solo
end

"inpolygon tests wether a point is in a polygon."
function inpolygon(p::Vec2, poly::Vector{Vec2})
	wn = 0
	xd(e, p) = (e[2].x - e[1].x) * (p.y - e[1].y) - (p.x - e[1].x) * (e[2].y - e[1].y)
	push!(poly, poly[1]) # close polygon
	n = length(poly)
	for i in 1:n-1
		if poly[i].y <= p.y
			if poly[i+1].y > p.y && xd(poly[i:i+1], p) > 0
				wn += 1
			end
		else
			if poly[i+1].y <= p.y && xd(poly[i:i+1], p) < 0
				wn -= 1
			end
		end
	end
	wn != 0
end

"plotalpha plots the precomputed α-shape of a group of particles."
function plotalpha(outer, inner, degen, solo, p::Vector{State})
	layers = Layer[]
	for u in outer
		push!(u, u[1])
		push!(layers, layer(x=[v.pos.x for v in p[u]], y=[v.pos.y for v in p[u]], Geom.point, Geom.path, Theme(default_color=colorant"red")))
	end
	for u in inner
		push!(u, u[1])
		push!(layers, layer(x=[v.pos.x for v in p[u]], y=[v.pos.y for v in p[u]], Geom.point, Geom.path, Theme(default_color=colorant"orange")))
	end
	for u in degen
		push!(u, u[1])
		push!(layers, layer(x=[v.pos.x for v in p[u]], y=[v.pos.y for v in p[u]], Geom.point, Geom.path, Theme(default_color=colorant"gray")))
	end
	push!(layers, layer(x=[v.pos.x for v in p[solo]], y=[v.pos.y for v in p[solo]], Geom.point, Theme(default_color=colorant"pink")))
	plot(layers...)
	# plot(
	# 	[layer(x=[v.pos.x for v in p[[id;id[1]]]], y=x[[id;id[1]],2], Geom.point, Geom.path, Theme(default_color=colorant"red")) for id in outer]...,
	# 	[layer(x=x[[id;id[1]],1], y=x[[id;id[1]],2], Geom.point, Geom.path, Theme(default_color=colorant"orange")) for id in inner]...,
	# 	[layer(x=x[collect(id),1], y=x[collect(id),2], Geom.point, Geom.path, Theme(default_color=colorant"gray")) for id in degen]...,
	# 	layer(x=x[solo,1], y=x[solo,2], Geom.point, Theme(default_color=colorant"pink")),
	# 	layer(x=x[:,1], y=x[:,2], Geom.point)
	# )
end

"makevidalpha save a time sequence of α-shape plots"
function makevidalpha(α::Real, p::Matrix{State})
	K = size(x,2)
	for k=1:K
		@printf("\r%3d%%", 100(k-1) / K)
		n = countnz([v.pos.x for v in p[:,k]])
		outer, inner, degen, solo = alphashape(α, p[1:n,k])
		h = plotalpha(outer, inner, degen, solo, p[1:n,k])
		draw(PNG(@sprintf("vid/%05d.png", k), 6inch, 6inch), h)
	end
	println("\r100%")
end

"state_dtype return the HDF5 compound datatype for a State (must be closed)."
function state_dtype()
	vt = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(Vec2))
	HDF5.h5t_insert(vt, "X", 0, HDF5.H5T_NATIVE_DOUBLE)
	HDF5.h5t_insert(vt, "Y", 8, HDF5.H5T_NATIVE_DOUBLE)
	st = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(State))
	HDF5.h5t_insert(st, "Pos", 0, vt)
	HDF5.h5t_insert(st, "Dir", 16, HDF5.H5T_NATIVE_DOUBLE)
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
	p
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

nothing
