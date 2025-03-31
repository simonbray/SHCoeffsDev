module Geotools
using Random
using Distributions
using LinearAlgebra
using GeometryBasics

export plegendre, ylm
export spherical_samples, rtp2xyz, xyz2rtp
export ray_intersects_triangle
export barycentric_coordinates
export load_tmesh

function plegendre(l::Int64, m::Int64, x::Float64)
    pmm = 1.0
    if m > 0
        omx2 = (1.0 - x) * (1.0 + x)
        fact = 1.0
        for _ ∈ 1:m
            pmm *= omx2 * fact / (fact + 1.0)
            fact += 2.0
        end
    end
    pmm = sqrt((2 * m + 1) * pmm / (4.0 * π))
    if m & 1 != 0
        pmm = -pmm
    end
    if l == m
        return pmm
    else
        pmmp1 = x * sqrt(2.0 * m + 3.0) * pmm
        if l == m + 1
            return pmmp1
        else
            pll = 0.0
            oldfact = sqrt(2.0 * m + 3.0)
            for ll ∈ (m+2):l
                fact = sqrt((4.0 * ll * ll - 1.0) / (ll * ll - m * m))
                pll = (x * pmmp1 - pmm / oldfact) * fact
                oldfact = fact
                pmm = pmmp1
                pmmp1 = pll
            end
            return pll
        end
    end
end


function ylm(l::Int64, m::Int64, θ::Float64, ϕ::Float64)
    m_sqrt2 = 1.41421356237309504880
    if m < 0
        result = m_sqrt2 * sin(-m * ϕ) * plegendre(l, -m, cos(θ))
    elseif m == 0
        result = plegendre(l, m, cos(θ))
    else
        result = m_sqrt2 * cos(m * ϕ) * plegendre(l, m, cos(θ))
    end
    return (result)
end


function spherical_samples(n_samples::Int64, lmax::Int64)
    Random.seed!(123)

    n = Int64(ceil(sqrt(n_samples)))

    vec3d = Vector{Vector{Float64}}(undef, n^2)
    ylm_values = Vector{Vector{Float64}}(undef, n^2)

    i = 1
    for a ∈ 0:n-1
        for b ∈ 0:n-1
            x = (Float64(a) + rand(Uniform(0, 1))) / n
            y = (Float64(b) + rand(Uniform(0, 1))) / n

            θ = 2.0 * acos(sqrt(1.0 - x))
            ϕ = 2.0 * π * y

            j = 1
            band = zeros((lmax + 1)^2)
            for l ∈ 0:lmax
                for m ∈ -l:l
                    band[j] = ylm(l, m, θ, ϕ)
                    j += 1
                end
            end

            vec3d[i] = [1.0, θ, ϕ]
            ylm_values[i] = band
            i += 1
        end
    end
    return (vec3d, ylm_values)
end


function rtp2xyz(vec::Vector{Float64})
    x = vec[1] * sin(vec[2]) * cos(vec[3])
    y = vec[1] * sin(vec[2]) * sin(vec[3])
    z = vec[1] * cos(vec[2])

    return ([x, y, z])
end


function xyz2rtp(vec::Vector{Float64})
    r = sqrt(vec[1]^2 + vec[2]^2 + vec[3]^2)
    θ = acos(vec[3] / r)
    ϕ = atan(vec[2], vec[1])

    return ([r, θ, ϕ])
end


function ray_intersects_triangle(ray::Vector{Float32}, triangle::Vector{Vector{Float32}})
    intersection = zeros(3)

    ϵ = 0.0000001

    v0 = triangle[1]
    v1 = triangle[2]
    v2 = triangle[3]

    edge1 = v1 - v0
    edge2 = v2 - v0

    h = cross(ray, edge2)
    a = dot(edge1, h)

    if -ϵ < a < ϵ
        return (false, intersection)
    end

    f = 1.0 / a
    s = -v0
    u = f * dot(s, h)

    if u < 0.0 || u > 1.0
        return (false, intersection)
    end

    q = cross(s, edge1)
    v = f * dot(ray, q)

    if v < 0.0 || u + v > 1.0
        return (false, intersection)
    end

    t = f * dot(edge2, q)
    if t > ϵ
        intersection .= ray * t
        return (true, intersection)
    end

    return (false, intersection)
end


function barycentric_coordinates(triangle::Vector{Vector{Float32}}, P::Vector{Float64})
    A, B, C = triangle

    ABC = length(cross(B - A, C - A))
    CAP = length(cross(C - A, C - P))
    ABP = length(cross(A - B, A - P))

    u = CAP / ABC
    v = ABP / ABC

    return (u, v)
end


function load_tmesh(in_file_name::String; facetype=GLTriangleFace, pointtype=Point3f, normalstype=Vec3f)
    in_file = open(in_file_name, "r")
    lines = readlines(in_file)
    close(in_file)

    n_points = parse(Int64, lines[1])

    points = pointtype[]
    point_normals = normalstype[]
    faces = facetype[]

    point_esp = eltype(normalstype)[]

    # read the data
    for (i, line) in enumerate(lines)
        tokens = split(line)
        if size(tokens) == (9,)
            numbers = parse.(eltype(pointtype), tokens[1:8])
            push!(points, pointtype(numbers[1:3]))
            push!(point_normals, pointtype(numbers[4:6]))
            push!(point_esp, numbers[8])
        elseif line == "3 0"
            line = lines[i+1:i+3]
            push!(faces, NgonFace{3,eltype(facetype)}(reinterpret(ZeroIndex{UInt32}, parse.(UInt32, line))))
        end
    end

    # Centralize the mesh
    x = 0.0
    y = 0.0
    z = 0.0
    for point ∈ points
        x += point[1]
        y += point[2]
        z += point[3]
    end
    x /= n_points
    y /= n_points
    z /= n_points

    points_centred = pointtype[]
    for point ∈ points
        push!(points_centred, pointtype([point[1] - x, point[2] - y, point[3] - z]))
    end

    return Mesh(points_centred, faces; normal=point_normals, ESP=point_esp)
end

end