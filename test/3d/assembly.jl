@testset "Triangular Assembly" begin
    freq  = 100.0                                   # Frequency                 [Hz]
    c     = 340.0                                   # Speed of sound            [m/s]
    k     = 2*π*freq/c                              # Wavenumber                [1/m]
    angles = [π/2 0.0]                              # Angles of incoming wave   [radians]
    radius = 1.0                                    # Radius of sphere_1m       [m]
    # Computing incident pressure
    # Relative path from the "runtests.jl" file (not the "meshing.jl" file)
    mesh_file = "../examples/meshes/sphere_1m"
    geometry_orders = [:linear,:quadratic]
    physics_orders = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic]
    errors = [0.037 0.0063 0.0073 0.0041 0.004;
              0.048 2.8e-5 0.0008 3.38e-5 4.9e-5]
    for (i,go) in enumerate(geometry_orders), (j,po) in enumerate(physics_orders)
        if go == :linear && po == :geometry
            continue
        end
        mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=go,physics_order=po)
        # The elements described by 'topology' should match the number of shape functions
        #===================================================================================
                                        Assembling
        ===================================================================================#
        Fp,Gp,Cp = assemble_parallel!(mesh,k,mesh.sources;progress=false)
        pI       = incoming_wave(angles,1.0,mesh.sources,k)
        #===================================================================================
                    Setting up a linear system and solving for the pressure
        ===================================================================================#
        Ap    = Fp + Diagonal(1.0 .- Cp)
        p_bem = gmres(Ap,pI)
        ## Computing analytical solution
        surface_angles  = acos.(-mesh.sources[1,:]/radius)
        p_analytical, _ = plane_wave_scattering_sphere(k,radius,1.0,surface_angles,1e-6)
        @test isapprox(norm(p_analytical - conj(p_bem))/norm(p_analytical),0.0,atol=errors[i,j])
    end
end

@testset "Quadrilateral Assembly" begin
    freq  = 100.0                                   # Frequency                 [Hz]
    c     = 340.0                                   # Speed of sound            [m/s]
    k     = 2*π*freq/c                              # Wavenumber                [1/m]
    angles = [π/2 0.0]                              # Angles of incoming wave   [radians]
    radius = 1.0                                    # Radius of sphere_1m       [m]
    # Computing incident pressure
    # Relative path from the "runtests.jl" file (not the "meshing.jl" file)
    mesh_file = "../examples/meshes/quad_sphere"
    geometry_orders = [:linear,:quadratic]
    physics_orders  = [:linear,:geometry,:discquadconstant,:discquadlinear,:discquadquadratic]
    errors = [0.012 0.012 0.014 0.009 0.008;
              0.004 2.8e-5 0.0019 0.003 3.1e-5]
    for (i,go) in enumerate(geometry_orders), (j,po) in enumerate(physics_orders)
        if go == :linear && po == :geometry
            continue
        end
        mesh = load3dQuadComsolMesh(mesh_file;geometry_order=go,physics_order=po)
        #===================================================================================
                                        Assembling
        ===================================================================================#
        Fp,Gp,Cp = assemble_parallel!(mesh,k,mesh.sources;progress=false)
        pI       = incoming_wave(angles,1.0,mesh.sources,k)
        #===================================================================================
                    Setting up a linear system and solving for the pressure
        ===================================================================================#
        Ap    = Fp + Diagonal(1.0 .- Cp)
        p_bem = gmres(Ap,pI)
        ## Computing analytical solution
        surface_angles  = acos.(-mesh.sources[1,:]/radius)
        p_analytical, _ = plane_wave_scattering_sphere(k,radius,1.0,surface_angles,1e-6)
        @test isapprox(norm(p_analytical - conj(p_bem))/norm(p_analytical),0.0,atol=errors[i,j])
    end
end
