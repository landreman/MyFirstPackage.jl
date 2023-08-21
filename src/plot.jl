using GLMakie

function plot(curve::Curve)
    println("foobar2 !")
    #fig = Figure(resolution = (2560, 1500))
    fig = Figure(resolution = (2850, 1800))
    #ax = Axis3(fig[1, 1])
    ax = LScene(fig[1, 1])
    t = range(0, 2π, length=200)
    r, drdt, d2rdt2, d3rdt3 = position_vector(curve, t)
    @show t
    @show size(r)
    lines!(r)
    #lines!(r[1, :], r[2, :], r[3, :])

    nvects = 40
    δ = 0.2
    t_vects = range(0, 2π, length=nvects + 1)[1 : end - 1]
    r, drdt, d2rdt2, d3rdt3 = position_vector(curve, t_vects)
    @show t_vects
    for j in 1:nvects
        tangent = drdt[:, j] / norm(drdt[:, j])

        # Compute the alternative normal vector
        #pre_pre_normal = r[:, j]
        pre_pre_normal = [0, 0, 1.0]
        pre_normal = pre_pre_normal - dot(tangent, pre_pre_normal) * tangent
        normal = pre_normal / norm(pre_normal)
        
        binormal = cross(tangent, normal)

        @assert dot(tangent, tangent) ≈ 1
        @assert dot(normal, normal) ≈ 1
        @assert dot(binormal, binormal) ≈ 1
        @assert abs(dot(tangent, normal)) < 1e-12
        @assert abs(dot(tangent, binormal)) < 1e-12
        @assert abs(dot(normal, binormal)) < 1e-12

        r_start = r[:, j]
        lines!([r_start (r_start + δ * tangent)], color=:red)
        lines!([r_start (r_start + δ * normal)], color=:green)
        lines!([r_start (r_start + δ * binormal)], color=:blue)
    end

    #fig
    #return lines(r)
    #zoom!(ax, 2.0)

    # Both of the next 2 lines are required to zoom in:
    zoom!(ax.scene, cameracontrols(ax.scene), 0.5)
    update_cam!(ax.scene, cameracontrols(ax.scene))

    fig
end