using IntegralEquations
using LinearAlgebra
using Plots

meshFile = "examples/meshes/sphere_1m_extremely_fine"
mesh = load3dTriangularComsolMesh(meshFile;geometryType="TriLinear")


element_interpolations = IntegralEquations.interpolate_elements(mesh)

S = 0.0
for i = 1:length(element_interpolations)
    S += dot(element_interpolations[i].jacobians,mesh.shape_function.weights)
end
S


@time IntegralEquations.interpolate_elements(mesh);


# quad_mesh = load3dTriangularComsolMesh(meshFile;geometryType="TriLinear")
coordinates,topology,entities = read_comsol_mesh("examples/meshes/quad_cylinder",QuadrilateralQuadratic9(3,3))
quadMesh = load3dQuadComsolMesh("examples/meshes/quad_cylinder")

Ql = QuadrilateralLinear(3,3)
IntegralEquations.set_nodal_interpolation!(Ql)
Ql4 = QuadrilateralLinear4(3,3)
IntegralEquations.set_nodal_interpolation!(Ql4)

Qq = QuadrilateralQuadratic(3,3)
IntegralEquations.set_nodal_interpolation!(Qq)
Qq9 = QuadrilateralQuadratic9(3,3)
IntegralEquations.set_nodal_interpolation!(Qq9)


tmp = IntegralEquations.get_element_normals(Qq9, coordinates, topology)