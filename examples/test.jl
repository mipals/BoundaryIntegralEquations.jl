using SparseArrays

Dt1 = BB.Dt₁
Dt2 = BB.Dt₂

n1 = nnz(Dt1)
n2 = nnz(Dt2)

element_connections, source_connections = IntegralEquations.connected_sources(mesh)

Dx,Dy,Dz = IntegralEquations.global_coordinate_shape_function_derivative(mesh);
Dt,Ds    = IntegralEquations.shape_function_derivatives(mesh);

@time Dt,Ds    = IntegralEquations.shape_function_derivatives(mesh);
@time Dx,Dy,Dz = IntegralEquations.global_coordinate_shape_function_derivative(mesh);


using Pkg
Pkg.add(PackageSpec(url="https://github.com/flatironinstitute/FMM3D.git", subdir="julia"))
