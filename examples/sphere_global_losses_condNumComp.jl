#==========================================================================================
                            Adding Related Packages
==========================================================================================#
using LinearAlgebra
using BoundaryIntegralEquations
using JLD2
#==========================================================================================
                                Loading Mesh
==========================================================================================#
# Triangular Meshes
mesh_path = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes");
mesh_files = ["sphere_1m_coarser","sphere_1m_coarse","sphere_1m","sphere_1m_fine","sphere_1m_finer"]

for freq in [100.0,500.0,1000.0]

    for mesh_file in mesh_files

        tri_mesh_file = joinpath(mesh_path,mesh_file);
        mesh = load3dTriangularComsolMesh(tri_mesh_file);
        #==========================================================================================
                                        Setting up constants
        ==========================================================================================#
        depth = 1;
        xyzb = mesh.sources;
        n = size(xyzb,2);
        #===========================================================================================
                                Computing Condition Numbers
        ===========================================================================================#
        data_file = joinpath(mesh_file, "results_$(n)DOFs_freq$(Int(freq)).jld2");
        println("Computing setup: $(data_file)")

        # Computing matrices
        LGM = LossyGlobalOuter(mesh,freq;fmm_on=false,depth=depth,progress=false);
        condGa  = cond(LGM.Ga);
        println("Conditon Number of Ga (n=$(n),freq=$(Int(freq))): $(condGa)")
        condHa  = cond(LGM.Ha);
        println("Conditon Number of Ha (n=$(n),freq=$(Int(freq))): $(condHa)")
        condGh  = cond(Matrix(LGM.Gh));
        println("Conditon Number of Gh (n=$(n),freq=$(Int(freq))): $(condGh)")
        condHh  = cond(Matrix(LGM.Hh));
        println("Conditon Number of Hh (n=$(n),freq=$(Int(freq))): $(condHh)")
        condGv  = cond(Matrix(LGM.Gv[1:n,1:n]));
        println("Conditon Number of Gv (n=$(n),freq=$(Int(freq))): $(condGv)")
        condHv  = cond(Matrix(LGM.Hv[1:n,1:n]));
        println("Conditon Number of Hv (n=$(n),freq=$(Int(freq))): $(condHv)")

        # Creating the different systems
        F = BoundaryIntegralEquations._full1(LGM)
        condF1  = cond(F);
        println("Conditon Number of 1x1 (n=$(n),freq=$(Int(freq))): $(condF1)")
        F  = BoundaryIntegralEquations._full4(LGM);
        condF4  = cond(F);
        println("Conditon Number of 4x4 (n=$(n),freq=$(Int(freq))): $(condF4)")
        F  = BoundaryIntegralEquations._full10(LGM);
        condF10 = cond(F);
        println("Conditon Number of 10x10 (n=$(n),freq=$(Int(freq))): $(condF10)")

        # Saving data
        jldsave("results_$(n)DOFs_freq$(Int(freq)).jld2", 
                condF1=condF1,
                condF4=condF4,
                condF10=condF10,
                condGa=condGa,
                condHa=condHa,
                condGh=condGh,
                condHh=condHh,
                condGv=condGv,
                condHv=condHv
                );
    end
end