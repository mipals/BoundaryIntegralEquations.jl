using IntegralEquations
using LinearAlgebra
using Plots

meshFile = "examples/meshes/cylinder"
mesh = load3dTriangularComsolMesh(meshFile)
mesh = load3dTriangularComsolMesh(meshFile;geometryType="TriLinear")
mesh = load3dTriangularComsolMesh(meshFile;physicsType="linear")

mesh = load3dTriangularComsolMesh(meshFile;physicsType="disctriconstant")
mesh = load3dTriangularComsolMesh(meshFile;physicsType="disctrilinear")
mesh = load3dTriangularComsolMesh(meshFile;physicsType="disctriquadratic")


element_interpolations = IntegralEquations.interpolate_elements(mesh)

S = 0.0
for i = 1:length(element_interpolations)
    S += dot(element_interpolations[i].jacobians,mesh.shape_function.weights)
end
S

@time IntegralEquations.interpolate_elements(mesh);


# quad_mesh = load3dTriangularComsolMesh(meshFile;geometryType="TriLinear")
quad_mesh = "examples/meshes/quad_cylinder"
coordinates,topology,entities = read_comsol_mesh(quad_mesh,QuadrilateralQuadratic9(3,3))
quadMesh = load3dQuadComsolMesh(quad_mesh)
quadMeshgl = load3dQuadComsolMesh(quad_mesh;geometryType="quadlinear")
quadMeshpl = load3dQuadComsolMesh(quad_mesh;physicsType="linear")

quadMeshDiscCon = load3dQuadComsolMesh(quad_mesh;physicsType="discquadconstant")
quadMeshDiscCon = load3dQuadComsolMesh(quad_mesh;physicsType="discquadlinear")
quadMeshDiscCon = load3dQuadComsolMesh(quad_mesh;physicsType="discquadquadratic")

@time IntegralEquations.interpolate_elements(quadMesh);
TMP = IntegralEquations.interpolate_elements(quadMesh)

Ql = QuadrilateralLinear(3,3)
IntegralEquations.set_nodal_interpolation!(Ql)
Ql4 = QuadrilateralLinear4(3,3)
IntegralEquations.set_nodal_interpolation!(Ql4)

Qq = QuadrilateralQuadratic(3,3)
IntegralEquations.set_nodal_interpolation!(Qq)
Qq9 = QuadrilateralQuadratic9(3,3)
IntegralEquations.set_nodal_interpolation!(Qq9)


normals = mesh.normals
tangentX = similar(normals)
tangentY = similar(normals)
IntegralEquations.tangents!(normals,tangentX,tangentY)

