#==========================================================================================
                                Using relevant packages
==========================================================================================#
using IntegralEquations
using LinearAlgebra
using Test
#==========================================================================================
                                    Creating Tests
==========================================================================================#
@testset "Curve Integration" begin
    n_elements = 100

    lin_circ    = IntegralEquations.mesh_circle(ContinuousCurveLinear(3),n_elements)
    lin_interps = IntegralEquations.interpolate_elements(lin_circ;n=4)

    quad_circ    = IntegralEquations.mesh_circle(ContinuousCurveQuadratic(3),n_elements)
    quad_interps = IntegralEquations.interpolate_elements(quad_circ;n=4)

    Slin  = 0.0
    Squad = 0.0
    for lin_interp ∈ lin_interps
        Slin  += sum(lin_interp.jacobian_mul_weights)
    end
    for quad_interp ∈ quad_interps
        Squad += sum(quad_interp.jacobian_mul_weights)
    end

    @test isapprox(Slin,2π,atol=5e-3)
    @test isapprox(Squad,2π,atol=1e-6)
end

@testset "BEM Assembly" begin
    n_elements = 100
    k = 1.0
    r = 1.0

    geometric_elements = (ContinuousCurveLinear(3),
                          ContinuousCurveQuadratic(3))
    physics_elements   = (DiscontinuousCurveConstant(3),
                          DiscontinuousCurveLinear(3),
                          DiscontinuousCurveQuadratic(3))

    for ge in geometric_elements, pe in physics_elements
        mesh = IntegralEquations.mesh_circle(ge,pe,n_elements;radius=r)
        F,G,C=assemble_parallel!(mesh,k,mesh.sources;n=4,progress=false,gOn=false);

        ϕ = angle.(mesh.sources[1,:] + im*mesh.sources[2,:])
        p = IntegralEquations.plane_wave_scattering_circle(ϕ,k*r,10)
        pB = (F + Diagonal(C .- 1.0))\exp.(im*k*mesh.sources[1,:])

        @test all(isapprox.(abs.(p),abs.(pB),atol=1e-3))
    end

    for ge in geometric_elements
        mesh = IntegralEquations.mesh_circle(ge,n_elements;radius=r)
        F,G,C=assemble_parallel!(mesh,k,mesh.sources;n=4,progress=false,gOn=false);

        ϕ = angle.(mesh.sources[1,:] + im*mesh.sources[2,:])
        p = IntegralEquations.plane_wave_scattering_circle(ϕ,k*r,10)
        pB = (F + Diagonal(C .- 1.0))\exp.(im*k*mesh.sources[1,:])

        @test all(isapprox.(abs.(p),abs.(pB),atol=1e-3))
    end

end
