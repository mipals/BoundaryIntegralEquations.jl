#==========================================================================================
                                Using relevant packages
==========================================================================================#
using Test
using IntegralEquations
using LinearAlgebra
using IterativeSolvers
#==========================================================================================
                            Testing mesh interpolation schemes
==========================================================================================#
@testset "Testing losses" begin
    mesh_file = "../examples/meshes/sphere_1m"
    geometry_orders = [:linear,:quadratic]
    #! Currently the dicontinuous elements does not work. Constant elements are also no-go.
    #! physics_orders  = [:linear,:quadratic,:disctrilinear,:disctriquadratic]
    physics_orders  = [:linear,:quadratic]
    for go in geometry_orders, po in physics_orders
        if go == :linear && po == :quadratic
            continue
        end
        #* Loading correct mesh
        mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=go,physics_order=po)
        #* Setting frequency
        freq = 100.0
        for fmm_on in [true false]
            #* Computing everything for the Kirchoffs' decomposition
            LGM = LossyGlobalOuter(mesh,freq;fmm_on=fmm_on,depth=1,n=3,progress=false)
            #* Defining velocity
            n = size(mesh.normals,2)
            v0 = [0.01*ones(n); zeros(2n)]
            #* Computing right-hand-side
            rhs = LGM.Ga*gmres(LGM.inner,(LGM.Dr*v0 - LGM.Nd'*gmres(LGM.Gv,LGM.Hv*v0)));
            #* Comuting the acoustical pressure by solving the linear system
            pa = gmres(LGM,rhs);
            #? Extracting the remaining modes from the acosutical pressure
            dpa  = -gmres(LGM.Ga,LGM.Ha*pa)     #! Normal derivative of acoustic pressure
            ph   = -LGM.tau_a/LGM.tau_h * pa    #! Thermal pressure
            dph  = -gmres(LGM.Gh,LGM.Hh*ph)     #! Normal derivative of thermal pressure
            v   = v0 - (LGM.phi_a*LGM.Dc*pa  +  #! Viscous velocity
                        LGM.phi_a*LGM.Nd*dpa +
                        LGM.phi_h*LGM.Dc*ph  +
                        LGM.phi_h*LGM.Nd*dph)
            dvn = -gmres(LGM.Gv,LGM.Hv*v)       #! normal derivative of viscous velocity
            #? Testing if constants are correct
            @test LGM.mu_a == LGM.phi_a - LGM.phi_h*LGM.tau_a/LGM.tau_h
            @test LGM.mu_h == LGM.phi_h*LGM.tau_a/LGM.tau_h
            #? All of these are solved quantities - So they're satisfied by construction
            @test LGM.Hh*ph ≈ -LGM.Gh*dph       #! Checking Acoustical BEM System
            @test LGM.Ha*pa ≈ -LGM.Ga*dpa       #! Checking Thermal BEM system
            @test LGM.Hv*v  ≈ -LGM.Gv*dvn       #! Checking viscous BEM system
            @test LGM.tau_a*pa ≈ -LGM.tau_h*ph  #! Checking Isothermal condition
            #? These are the extra conditions
            @test LGM.phi_a*(LGM.Dc*pa + LGM.Nd*dpa) + LGM.phi_h*(LGM.Dc*ph + LGM.Nd*dph) + v ≈ v0
            @test LGM.phi_a*LGM.Dc*pa + LGM.phi_a*LGM.Nd*dpa + LGM.phi_h*LGM.Dc*ph + LGM.phi_h*LGM.Nd*dph + v ≈ v0
            #! Currently the null-divergence is not satisfied using this approach.
            @test_broken LGM.Dr*v ≈ -LGM.Nd'*dvn
        end
    end
end
