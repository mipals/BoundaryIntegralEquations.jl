```@meta
EditURL = "../../../examples/3d_cube_wave_anechoic.jl"
```

# Cube with anechoic condition (Interior)
# Importing related packages

````julia
using BoundaryIntegralEquations # For BIEs
using IterativeSolvers          # For gmres
using LinearAlgebra             # For Diagonal
using Plots                     # For 2d plots
using Meshes                    # For 3d mesh plots
import GLMakie as wgl
````

# Setting up constants

````julia
frequency = 54.59;                              # Frequency                [Hz]
c  = 343;                                       # Speed up sound           [m/s]
ρ₀ = 1.21;                                      # Mean density             [kg/m^3]
Z₀ = ρ₀*c;                                      # Characteristic impedance [Rayl]
P₀ = 1.0;                                       # Pressure of plane wave   [m/s]
l  = 1.0;                                       # Length of cube           [m]
k  = 2π*frequency/c;                            # Computing the wavenumber
````

# Loading Mesh
Loading and visualizing the triangular cube mesh

````julia
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","1m_cube_extra_coarse");
physics_orders  = [:linear,:geometry,:disctriconstant,:disctrilinear,:disctriquadratic];
mesh = load3dTriangularComsolMesh(mesh_file;geometry_order=:linear, physics_order=physics_orders[5])
bc_pre = [1] .- 1; # The .-1 is due to COMSOL 0-indexing of exported entities
bc_ane = [6] .- 1; # The .-1 is due to COMSOL 0-indexing of exported entities
````

Creating simple meshes for full mesh and boundary condition

````julia
simple_mesh = create_bc_simple_mesh(mesh,[bc_pre; bc_ane],false);
simple_pre  = create_bc_simple_mesh(mesh,bc_pre);
simple_ana  = create_bc_simple_mesh(mesh,bc_ane);
````

We now plot the mesh with the pressure condition shown in red and the anechoic condition shown in blue

````julia
fig = viz(simple_pre;showfacets=true,color=:red)
viz!(simple_mesh;showfacets=true,alpha=0.1)
viz!(simple_ana;showfacets=true,color=:blue)
````

![](3d_cube_wave_anechoic_viz.png)
# Analytical Solution
The analytical description of the interior pressure of a planewave is given by
```math
 p(x) = P_0\exp(-\mathrm{i}kx),
```
where ``P_0`` is the magnitude of the planewave.
We now compute the analytical expression is computed at the points

````julia
x_ana = collect(0.0:0.01:1)
p_analytical = P₀*exp.(-im*k*x_ana);
````

## Solution using the dense BEM
We start by solving the BEM system using dense matrices. For this we need to first assemble the matrices

````julia
F,G,C = assemble_parallel!(mesh,k,mesh.sources,n=4,m=4,progress=false);
H  = Diagonal(C) - F; # Adding the integral free term
````

Now we find the correct indicies

````julia
bc_pre = 0 .== mesh.entities;
bc_ane = 5 .== mesh.entities;
bc1 = sort(unique(mesh.physics_topology[:,bc_pre]));
bc2 = sort(unique(mesh.physics_topology[:,bc_ane]));
````

First we define the ``p=1`` bc at ``x=0``

````julia
ps = zeros(ComplexF64,length(C));
ps[bc1] .= 1.0;
bs = -(H*ps); # Computing right-hand side
````

Then we define the anechoic condition at ``x=1``.

````julia
Z = zeros(ComplexF64,length(C));
Z[bc1] .=  1.0;  # Unknown velocities
Z[bc2] .= -im*k; # Impedance condition
````

````julia
Hmask = ones(ComplexF64,length(C));
Hmask[bc1] .= 0.0; # Known pressures
A = H*Diagonal(Hmask) + G*Diagonal(-Z);
````

Using this we can solve for the surface pressure

````julia
z_bem = gmres(A,bs);
````

Extracting correct surface pressures

````julia
p_bem = copy(z_bem);
p_bem[bc1] .= 1.0;
````

Extracting correct surface velocities

````julia
v_bem = zeros(ComplexF64,length(C));
v_bem[bc1] = z_bem[bc1];
v_bem[bc2] = p_bem[bc2]*(-im*k);
````

## Solution using the FMM-BEM
First we define the two operators

````julia
Gf = FMMGOperator(mesh,k,depth=2);
Ff = FMMFOperator(mesh,k,depth=2);
Hf = 0.5I - Ff;
````

Then we compute the right-hand side and the linear system similar to before. Note that the system matrix will be represented by a collection of linear maps

````julia
bf = -(Hf*ps);
Af = Hf*Diagonal(Hmask) + Gf*Diagonal(-Z)
````

````
1344×1344 LinearMaps.LinearCombination{ComplexF64} with 2 maps:
  1344×1344 LinearMaps.CompositeMap{ComplexF64} with 2 maps:
    1344×1344 LinearMaps.WrappedMap{ComplexF64} of
      1344×1344 LinearAlgebra.Diagonal{ComplexF64, Vector{ComplexF64}}
    1344×1344 LinearMaps.LinearCombination{ComplexF64} with 2 maps:
      1344×1344 LinearMaps.UniformScalingMap{Float64} with scaling factor: 0.5
      1344×1344 LinearMaps.ScaledMap{ComplexF64} with scale: -1 of
        1344×1344 BoundaryIntegralEquations.FMMFOperator{ComplexF64}
  1344×1344 LinearMaps.CompositeMap{ComplexF64} with 2 maps:
    1344×1344 LinearMaps.WrappedMap{ComplexF64} of
      1344×1344 LinearAlgebra.Diagonal{ComplexF64, Vector{ComplexF64}}
    1344×1344 BoundaryIntegralEquations.FMMGOperator{ComplexF64}
````

Similary to before we can now solve the linear system using an iterative scheme

````julia
z_fmm = gmres(Af,bf);
````

Finally we extract the pressure (and assert the pressures at boundary 1)

````julia
p_fmm = copy(z_fmm);
p_fmm[bc1] .= 1.0;
````

as well as the surface velocities

````julia
v_fmm = zeros(ComplexF64,length(C));
v_fmm[bc1] = z_fmm[bc1];
v_fmm[bc2] = p_fmm[bc2]*(-im*k);
````

# Field evaluations
In order to compare the results with the analytical solution we must use the found surface pressures to compute pressure at the interior field points. We therefore start by creating a matrix of points (as columns). Thi is done by chosing ``N`` linearly spaced points defined by ``(x,0.5,0.5)`` with ``x\in[0.1,0.9]``.

````julia
N = 20;
x = collect(0.1:(0.9-0.1)/(N-1):0.9);
y = 0.5ones(N);
z = 0.5ones(N);
X_fieldpoints = [x'; y'; z'];
````

Using the described field-points we can assemble the (dense) matrices for the field point evaluation

````julia
Fp,Gp,Cp = assemble_parallel!(mesh,k,X_fieldpoints,n=2,m=2,progress=false);
````

Using the matrices we can evalute the pressure at the field points using the previously found surface pressure

````julia
p_field  = Fp*p_bem + Gp*v_bem;
pf_field = Fp*p_fmm + Gp*v_fmm;
````

Plotting the field point pressure it can be seen that all 3 methods find the correct pressures

````julia
plot(x_ana,real.(p_analytical/P₀),  label="Re-Analytical")
scatter!(x,real.(p_field/P₀),   label="Re-Full", markershape=:rect)
scatter!(x,real.(pf_field/P₀),  label="Re-FMM")
plot!(x_ana,imag.(p_analytical/P₀), label="Im-Analytical")
scatter!(x,imag.(p_field/P₀),   label="Im-Full", markershape=:rect)
scatter!(x,imag.(pf_field/P₀),  label="Im-FMM")
xlabel!("r/a")
````

```@raw html
<?xml version="1.0" encoding="utf-8"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="600" height="400" viewBox="0 0 2400 1600">
<defs>
  <clipPath id="clip720">
    <rect x="0" y="0" width="2400" height="1600"/>
  </clipPath>
</defs>
<path clip-path="url(#clip720)" d="M0 1600 L2400 1600 L2400 0 L0 0  Z" fill="#ffffff" fill-rule="evenodd" fill-opacity="1"/>
<defs>
  <clipPath id="clip721">
    <rect x="480" y="0" width="1681" height="1600"/>
  </clipPath>
</defs>
<path clip-path="url(#clip720)" d="M192.941 1423.18 L2352.76 1423.18 L2352.76 47.2441 L192.941 47.2441  Z" fill="#ffffff" fill-rule="evenodd" fill-opacity="1"/>
<defs>
  <clipPath id="clip722">
    <rect x="192" y="47" width="2161" height="1377"/>
  </clipPath>
</defs>
<polyline clip-path="url(#clip722)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="254.067,1423.18 254.067,47.2441 "/>
<polyline clip-path="url(#clip722)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="763.458,1423.18 763.458,47.2441 "/>
<polyline clip-path="url(#clip722)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="1272.85,1423.18 1272.85,47.2441 "/>
<polyline clip-path="url(#clip722)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="1782.24,1423.18 1782.24,47.2441 "/>
<polyline clip-path="url(#clip722)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="2291.63,1423.18 2291.63,47.2441 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="192.941,1423.18 2352.76,1423.18 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="254.067,1423.18 254.067,1404.28 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="763.458,1423.18 763.458,1404.28 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="1272.85,1423.18 1272.85,1404.28 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="1782.24,1423.18 1782.24,1404.28 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="2291.63,1423.18 2291.63,1404.28 "/>
<path clip-path="url(#clip720)" d="M216.371 1454.1 Q212.76 1454.1 210.931 1457.66 Q209.126 1461.2 209.126 1468.33 Q209.126 1475.44 210.931 1479.01 Q212.76 1482.55 216.371 1482.55 Q220.005 1482.55 221.811 1479.01 Q223.639 1475.44 223.639 1468.33 Q223.639 1461.2 221.811 1457.66 Q220.005 1454.1 216.371 1454.1 M216.371 1450.39 Q222.181 1450.39 225.237 1455 Q228.315 1459.58 228.315 1468.33 Q228.315 1477.06 225.237 1481.67 Q222.181 1486.25 216.371 1486.25 Q210.561 1486.25 207.482 1481.67 Q204.427 1477.06 204.427 1468.33 Q204.427 1459.58 207.482 1455 Q210.561 1450.39 216.371 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M236.533 1479.7 L241.417 1479.7 L241.417 1485.58 L236.533 1485.58 L236.533 1479.7 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M261.602 1454.1 Q257.991 1454.1 256.162 1457.66 Q254.357 1461.2 254.357 1468.33 Q254.357 1475.44 256.162 1479.01 Q257.991 1482.55 261.602 1482.55 Q265.236 1482.55 267.042 1479.01 Q268.871 1475.44 268.871 1468.33 Q268.871 1461.2 267.042 1457.66 Q265.236 1454.1 261.602 1454.1 M261.602 1450.39 Q267.412 1450.39 270.468 1455 Q273.547 1459.58 273.547 1468.33 Q273.547 1477.06 270.468 1481.67 Q267.412 1486.25 261.602 1486.25 Q255.792 1486.25 252.713 1481.67 Q249.658 1477.06 249.658 1468.33 Q249.658 1459.58 252.713 1455 Q255.792 1450.39 261.602 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M291.764 1454.1 Q288.153 1454.1 286.324 1457.66 Q284.519 1461.2 284.519 1468.33 Q284.519 1475.44 286.324 1479.01 Q288.153 1482.55 291.764 1482.55 Q295.398 1482.55 297.204 1479.01 Q299.033 1475.44 299.033 1468.33 Q299.033 1461.2 297.204 1457.66 Q295.398 1454.1 291.764 1454.1 M291.764 1450.39 Q297.574 1450.39 300.63 1455 Q303.708 1459.58 303.708 1468.33 Q303.708 1477.06 300.63 1481.67 Q297.574 1486.25 291.764 1486.25 Q285.954 1486.25 282.875 1481.67 Q279.82 1477.06 279.82 1468.33 Q279.82 1459.58 282.875 1455 Q285.954 1450.39 291.764 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M726.259 1454.1 Q722.648 1454.1 720.819 1457.66 Q719.014 1461.2 719.014 1468.33 Q719.014 1475.44 720.819 1479.01 Q722.648 1482.55 726.259 1482.55 Q729.893 1482.55 731.699 1479.01 Q733.527 1475.44 733.527 1468.33 Q733.527 1461.2 731.699 1457.66 Q729.893 1454.1 726.259 1454.1 M726.259 1450.39 Q732.069 1450.39 735.125 1455 Q738.203 1459.58 738.203 1468.33 Q738.203 1477.06 735.125 1481.67 Q732.069 1486.25 726.259 1486.25 Q720.449 1486.25 717.37 1481.67 Q714.315 1477.06 714.315 1468.33 Q714.315 1459.58 717.37 1455 Q720.449 1450.39 726.259 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M746.421 1479.7 L751.305 1479.7 L751.305 1485.58 L746.421 1485.58 L746.421 1479.7 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M765.518 1481.64 L781.837 1481.64 L781.837 1485.58 L759.893 1485.58 L759.893 1481.64 Q762.555 1478.89 767.138 1474.26 Q771.745 1469.61 772.925 1468.27 Q775.171 1465.74 776.05 1464.01 Q776.953 1462.25 776.953 1460.56 Q776.953 1457.8 775.009 1456.07 Q773.087 1454.33 769.986 1454.33 Q767.787 1454.33 765.333 1455.09 Q762.902 1455.86 760.125 1457.41 L760.125 1452.69 Q762.949 1451.55 765.402 1450.97 Q767.856 1450.39 769.893 1450.39 Q775.263 1450.39 778.458 1453.08 Q781.652 1455.77 781.652 1460.26 Q781.652 1462.39 780.842 1464.31 Q780.055 1466.2 777.949 1468.8 Q777.37 1469.47 774.268 1472.69 Q771.166 1475.88 765.518 1481.64 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M791.698 1451.02 L810.055 1451.02 L810.055 1454.96 L795.981 1454.96 L795.981 1463.43 Q796.999 1463.08 798.018 1462.92 Q799.036 1462.73 800.055 1462.73 Q805.842 1462.73 809.222 1465.9 Q812.601 1469.08 812.601 1474.49 Q812.601 1480.07 809.129 1483.17 Q805.657 1486.25 799.337 1486.25 Q797.161 1486.25 794.893 1485.88 Q792.648 1485.51 790.24 1484.77 L790.24 1480.07 Q792.323 1481.2 794.546 1481.76 Q796.768 1482.32 799.245 1482.32 Q803.249 1482.32 805.587 1480.21 Q807.925 1478.1 807.925 1474.49 Q807.925 1470.88 805.587 1468.77 Q803.249 1466.67 799.245 1466.67 Q797.37 1466.67 795.495 1467.08 Q793.643 1467.5 791.698 1468.38 L791.698 1451.02 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1235.15 1454.1 Q1231.54 1454.1 1229.71 1457.66 Q1227.91 1461.2 1227.91 1468.33 Q1227.91 1475.44 1229.71 1479.01 Q1231.54 1482.55 1235.15 1482.55 Q1238.79 1482.55 1240.59 1479.01 Q1242.42 1475.44 1242.42 1468.33 Q1242.42 1461.2 1240.59 1457.66 Q1238.79 1454.1 1235.15 1454.1 M1235.15 1450.39 Q1240.96 1450.39 1244.02 1455 Q1247.1 1459.58 1247.1 1468.33 Q1247.1 1477.06 1244.02 1481.67 Q1240.96 1486.25 1235.15 1486.25 Q1229.34 1486.25 1226.26 1481.67 Q1223.21 1477.06 1223.21 1468.33 Q1223.21 1459.58 1226.26 1455 Q1229.34 1450.39 1235.15 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1255.31 1479.7 L1260.2 1479.7 L1260.2 1485.58 L1255.31 1485.58 L1255.31 1479.7 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1270.43 1451.02 L1288.79 1451.02 L1288.79 1454.96 L1274.71 1454.96 L1274.71 1463.43 Q1275.73 1463.08 1276.75 1462.92 Q1277.77 1462.73 1278.79 1462.73 Q1284.57 1462.73 1287.95 1465.9 Q1291.33 1469.08 1291.33 1474.49 Q1291.33 1480.07 1287.86 1483.17 Q1284.39 1486.25 1278.07 1486.25 Q1275.89 1486.25 1273.62 1485.88 Q1271.38 1485.51 1268.97 1484.77 L1268.97 1480.07 Q1271.05 1481.2 1273.28 1481.76 Q1275.5 1482.32 1277.98 1482.32 Q1281.98 1482.32 1284.32 1480.21 Q1286.66 1478.1 1286.66 1474.49 Q1286.66 1470.88 1284.32 1468.77 Q1281.98 1466.67 1277.98 1466.67 Q1276.1 1466.67 1274.23 1467.08 Q1272.37 1467.5 1270.43 1468.38 L1270.43 1451.02 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1310.54 1454.1 Q1306.93 1454.1 1305.11 1457.66 Q1303.3 1461.2 1303.3 1468.33 Q1303.3 1475.44 1305.11 1479.01 Q1306.93 1482.55 1310.54 1482.55 Q1314.18 1482.55 1315.98 1479.01 Q1317.81 1475.44 1317.81 1468.33 Q1317.81 1461.2 1315.98 1457.66 Q1314.18 1454.1 1310.54 1454.1 M1310.54 1450.39 Q1316.35 1450.39 1319.41 1455 Q1322.49 1459.58 1322.49 1468.33 Q1322.49 1477.06 1319.41 1481.67 Q1316.35 1486.25 1310.54 1486.25 Q1304.73 1486.25 1301.66 1481.67 Q1298.6 1477.06 1298.6 1468.33 Q1298.6 1459.58 1301.66 1455 Q1304.73 1450.39 1310.54 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1745.04 1454.1 Q1741.43 1454.1 1739.6 1457.66 Q1737.79 1461.2 1737.79 1468.33 Q1737.79 1475.44 1739.6 1479.01 Q1741.43 1482.55 1745.04 1482.55 Q1748.67 1482.55 1750.48 1479.01 Q1752.31 1475.44 1752.31 1468.33 Q1752.31 1461.2 1750.48 1457.66 Q1748.67 1454.1 1745.04 1454.1 M1745.04 1450.39 Q1750.85 1450.39 1753.91 1455 Q1756.98 1459.58 1756.98 1468.33 Q1756.98 1477.06 1753.91 1481.67 Q1750.85 1486.25 1745.04 1486.25 Q1739.23 1486.25 1736.15 1481.67 Q1733.1 1477.06 1733.1 1468.33 Q1733.1 1459.58 1736.15 1455 Q1739.23 1450.39 1745.04 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1765.2 1479.7 L1770.09 1479.7 L1770.09 1485.58 L1765.2 1485.58 L1765.2 1479.7 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1779.09 1451.02 L1801.31 1451.02 L1801.31 1453.01 L1788.77 1485.58 L1783.88 1485.58 L1795.69 1454.96 L1779.09 1454.96 L1779.09 1451.02 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1810.48 1451.02 L1828.84 1451.02 L1828.84 1454.96 L1814.76 1454.96 L1814.76 1463.43 Q1815.78 1463.08 1816.8 1462.92 Q1817.82 1462.73 1818.84 1462.73 Q1824.62 1462.73 1828 1465.9 Q1831.38 1469.08 1831.38 1474.49 Q1831.38 1480.07 1827.91 1483.17 Q1824.44 1486.25 1818.12 1486.25 Q1815.94 1486.25 1813.67 1485.88 Q1811.43 1485.51 1809.02 1484.77 L1809.02 1480.07 Q1811.1 1481.2 1813.33 1481.76 Q1815.55 1482.32 1818.03 1482.32 Q1822.03 1482.32 1824.37 1480.21 Q1826.71 1478.1 1826.71 1474.49 Q1826.71 1470.88 1824.37 1468.77 Q1822.03 1466.67 1818.03 1466.67 Q1816.15 1466.67 1814.28 1467.08 Q1812.42 1467.5 1810.48 1468.38 L1810.48 1451.02 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M2243.7 1481.64 L2251.34 1481.64 L2251.34 1455.28 L2243.03 1456.95 L2243.03 1452.69 L2251.29 1451.02 L2255.97 1451.02 L2255.97 1481.64 L2263.61 1481.64 L2263.61 1485.58 L2243.7 1485.58 L2243.7 1481.64 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M2273.05 1479.7 L2277.94 1479.7 L2277.94 1485.58 L2273.05 1485.58 L2273.05 1479.7 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M2298.12 1454.1 Q2294.51 1454.1 2292.68 1457.66 Q2290.88 1461.2 2290.88 1468.33 Q2290.88 1475.44 2292.68 1479.01 Q2294.51 1482.55 2298.12 1482.55 Q2301.76 1482.55 2303.56 1479.01 Q2305.39 1475.44 2305.39 1468.33 Q2305.39 1461.2 2303.56 1457.66 Q2301.76 1454.1 2298.12 1454.1 M2298.12 1450.39 Q2303.93 1450.39 2306.99 1455 Q2310.07 1459.58 2310.07 1468.33 Q2310.07 1477.06 2306.99 1481.67 Q2303.93 1486.25 2298.12 1486.25 Q2292.31 1486.25 2289.23 1481.67 Q2286.18 1477.06 2286.18 1468.33 Q2286.18 1459.58 2289.23 1455 Q2292.31 1450.39 2298.12 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M2328.28 1454.1 Q2324.67 1454.1 2322.84 1457.66 Q2321.04 1461.2 2321.04 1468.33 Q2321.04 1475.44 2322.84 1479.01 Q2324.67 1482.55 2328.28 1482.55 Q2331.92 1482.55 2333.72 1479.01 Q2335.55 1475.44 2335.55 1468.33 Q2335.55 1461.2 2333.72 1457.66 Q2331.92 1454.1 2328.28 1454.1 M2328.28 1450.39 Q2334.09 1450.39 2337.15 1455 Q2340.23 1459.58 2340.23 1468.33 Q2340.23 1477.06 2337.15 1481.67 Q2334.09 1486.25 2328.28 1486.25 Q2322.47 1486.25 2319.4 1481.67 Q2316.34 1477.06 2316.34 1468.33 Q2316.34 1459.58 2319.4 1455 Q2322.47 1450.39 2328.28 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1255.29 1537.87 Q1254.31 1537.3 1253.13 1537.04 Q1251.98 1536.76 1250.58 1536.76 Q1245.62 1536.76 1242.95 1540 Q1240.3 1543.22 1240.3 1549.27 L1240.3 1568.04 L1234.42 1568.04 L1234.42 1532.4 L1240.3 1532.4 L1240.3 1537.93 Q1242.15 1534.69 1245.11 1533.13 Q1248.07 1531.54 1252.3 1531.54 Q1252.91 1531.54 1253.64 1531.63 Q1254.37 1531.7 1255.26 1531.85 L1255.29 1537.87 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1271.85 1520.52 L1277.26 1520.52 L1260.71 1574.09 L1255.29 1574.09 L1271.85 1520.52 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M1299.6 1550.12 Q1292.5 1550.12 1289.77 1551.75 Q1287.03 1553.37 1287.03 1557.29 Q1287.03 1560.4 1289.06 1562.25 Q1291.13 1564.07 1294.67 1564.07 Q1299.54 1564.07 1302.46 1560.63 Q1305.42 1557.16 1305.42 1551.43 L1305.42 1550.12 L1299.6 1550.12 M1311.28 1547.71 L1311.28 1568.04 L1305.42 1568.04 L1305.42 1562.63 Q1303.42 1565.88 1300.43 1567.44 Q1297.44 1568.97 1293.11 1568.97 Q1287.63 1568.97 1284.39 1565.91 Q1281.17 1562.82 1281.17 1557.67 Q1281.17 1551.65 1285.18 1548.6 Q1289.22 1545.54 1297.21 1545.54 L1305.42 1545.54 L1305.42 1544.97 Q1305.42 1540.93 1302.75 1538.73 Q1300.11 1536.5 1295.3 1536.5 Q1292.25 1536.5 1289.35 1537.23 Q1286.45 1537.97 1283.78 1539.43 L1283.78 1534.02 Q1287 1532.78 1290.02 1532.17 Q1293.04 1531.54 1295.91 1531.54 Q1303.64 1531.54 1307.46 1535.55 Q1311.28 1539.56 1311.28 1547.71 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><polyline clip-path="url(#clip722)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="192.941,1143.54 2352.76,1143.54 "/>
<polyline clip-path="url(#clip722)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="192.941,791.086 2352.76,791.086 "/>
<polyline clip-path="url(#clip722)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="192.941,438.636 2352.76,438.636 "/>
<polyline clip-path="url(#clip722)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="192.941,86.1857 2352.76,86.1857 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="192.941,1423.18 192.941,47.2441 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="192.941,1143.54 211.838,1143.54 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="192.941,791.086 211.838,791.086 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="192.941,438.636 211.838,438.636 "/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="192.941,86.1857 211.838,86.1857 "/>
<path clip-path="url(#clip720)" d="M50.9921 1143.99 L80.6679 1143.99 L80.6679 1147.92 L50.9921 1147.92 L50.9921 1143.99 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M100.76 1129.33 Q97.1493 1129.33 95.3206 1132.9 Q93.515 1136.44 93.515 1143.57 Q93.515 1150.68 95.3206 1154.24 Q97.1493 1157.78 100.76 1157.78 Q104.395 1157.78 106.2 1154.24 Q108.029 1150.68 108.029 1143.57 Q108.029 1136.44 106.2 1132.9 Q104.395 1129.33 100.76 1129.33 M100.76 1125.63 Q106.571 1125.63 109.626 1130.24 Q112.705 1134.82 112.705 1143.57 Q112.705 1152.3 109.626 1156.9 Q106.571 1161.49 100.76 1161.49 Q94.9502 1161.49 91.8715 1156.9 Q88.816 1152.3 88.816 1143.57 Q88.816 1134.82 91.8715 1130.24 Q94.9502 1125.63 100.76 1125.63 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M120.922 1154.94 L125.807 1154.94 L125.807 1160.82 L120.922 1160.82 L120.922 1154.94 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M136.038 1126.26 L154.394 1126.26 L154.394 1130.19 L140.32 1130.19 L140.32 1138.66 Q141.339 1138.32 142.357 1138.15 Q143.376 1137.97 144.394 1137.97 Q150.181 1137.97 153.561 1141.14 Q156.941 1144.31 156.941 1149.73 Q156.941 1155.31 153.468 1158.41 Q149.996 1161.49 143.677 1161.49 Q141.501 1161.49 139.232 1161.12 Q136.987 1160.75 134.58 1160.01 L134.58 1155.31 Q136.663 1156.44 138.885 1157 Q141.107 1157.55 143.584 1157.55 Q147.589 1157.55 149.927 1155.45 Q152.265 1153.34 152.265 1149.73 Q152.265 1146.12 149.927 1144.01 Q147.589 1141.9 143.584 1141.9 Q141.709 1141.9 139.834 1142.32 Q137.982 1142.74 136.038 1143.62 L136.038 1126.26 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M99.765 776.885 Q96.1539 776.885 94.3252 780.45 Q92.5197 783.991 92.5197 791.121 Q92.5197 798.227 94.3252 801.792 Q96.1539 805.334 99.765 805.334 Q103.399 805.334 105.205 801.792 Q107.033 798.227 107.033 791.121 Q107.033 783.991 105.205 780.45 Q103.399 776.885 99.765 776.885 M99.765 773.181 Q105.575 773.181 108.631 777.788 Q111.709 782.371 111.709 791.121 Q111.709 799.848 108.631 804.454 Q105.575 809.037 99.765 809.037 Q93.9549 809.037 90.8762 804.454 Q87.8206 799.848 87.8206 791.121 Q87.8206 782.371 90.8762 777.788 Q93.9549 773.181 99.765 773.181 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M119.927 802.486 L124.811 802.486 L124.811 808.366 L119.927 808.366 L119.927 802.486 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M144.996 776.885 Q141.385 776.885 139.556 780.45 Q137.751 783.991 137.751 791.121 Q137.751 798.227 139.556 801.792 Q141.385 805.334 144.996 805.334 Q148.63 805.334 150.436 801.792 Q152.265 798.227 152.265 791.121 Q152.265 783.991 150.436 780.45 Q148.63 776.885 144.996 776.885 M144.996 773.181 Q150.806 773.181 153.862 777.788 Q156.941 782.371 156.941 791.121 Q156.941 799.848 153.862 804.454 Q150.806 809.037 144.996 809.037 Q139.186 809.037 136.107 804.454 Q133.052 799.848 133.052 791.121 Q133.052 782.371 136.107 777.788 Q139.186 773.181 144.996 773.181 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M100.76 424.435 Q97.1493 424.435 95.3206 427.999 Q93.515 431.541 93.515 438.671 Q93.515 445.777 95.3206 449.342 Q97.1493 452.883 100.76 452.883 Q104.395 452.883 106.2 449.342 Q108.029 445.777 108.029 438.671 Q108.029 431.541 106.2 427.999 Q104.395 424.435 100.76 424.435 M100.76 420.731 Q106.571 420.731 109.626 425.337 Q112.705 429.921 112.705 438.671 Q112.705 447.397 109.626 452.004 Q106.571 456.587 100.76 456.587 Q94.9502 456.587 91.8715 452.004 Q88.816 447.397 88.816 438.671 Q88.816 429.921 91.8715 425.337 Q94.9502 420.731 100.76 420.731 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M120.922 450.036 L125.807 450.036 L125.807 455.916 L120.922 455.916 L120.922 450.036 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M136.038 421.356 L154.394 421.356 L154.394 425.291 L140.32 425.291 L140.32 433.763 Q141.339 433.416 142.357 433.254 Q143.376 433.069 144.394 433.069 Q150.181 433.069 153.561 436.24 Q156.941 439.411 156.941 444.828 Q156.941 450.407 153.468 453.508 Q149.996 456.587 143.677 456.587 Q141.501 456.587 139.232 456.217 Q136.987 455.846 134.58 455.106 L134.58 450.407 Q136.663 451.541 138.885 452.096 Q141.107 452.652 143.584 452.652 Q147.589 452.652 149.927 450.546 Q152.265 448.439 152.265 444.828 Q152.265 441.217 149.927 439.11 Q147.589 437.004 143.584 437.004 Q141.709 437.004 139.834 437.421 Q137.982 437.837 136.038 438.717 L136.038 421.356 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M90.5752 99.5305 L98.2141 99.5305 L98.2141 73.1649 L89.904 74.8316 L89.904 70.5723 L98.1678 68.9057 L102.844 68.9057 L102.844 99.5305 L110.483 99.5305 L110.483 103.466 L90.5752 103.466 L90.5752 99.5305 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M119.927 97.5861 L124.811 97.5861 L124.811 103.466 L119.927 103.466 L119.927 97.5861 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M144.996 71.9844 Q141.385 71.9844 139.556 75.5492 Q137.751 79.0908 137.751 86.2204 Q137.751 93.3268 139.556 96.8916 Q141.385 100.433 144.996 100.433 Q148.63 100.433 150.436 96.8916 Q152.265 93.3268 152.265 86.2204 Q152.265 79.0908 150.436 75.5492 Q148.63 71.9844 144.996 71.9844 M144.996 68.2807 Q150.806 68.2807 153.862 72.8871 Q156.941 77.4704 156.941 86.2204 Q156.941 94.9472 153.862 99.5537 Q150.806 104.137 144.996 104.137 Q139.186 104.137 136.107 99.5537 Q133.052 94.9472 133.052 86.2204 Q133.052 77.4704 136.107 72.8871 Q139.186 68.2807 144.996 68.2807 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><polyline clip-path="url(#clip722)" style="stroke:#009af9; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="254.067,86.1857 274.443,86.2209 294.819,86.3266 315.194,86.5029 335.57,86.7495 355.946,87.0666 376.321,87.4541 396.697,87.912 417.072,88.4401 437.448,89.0386 457.824,89.7072 478.199,90.446 498.575,91.2548 518.95,92.1337 539.326,93.0824 559.702,94.1009 580.077,95.1891 600.453,96.3469 620.829,97.5742 641.204,98.8708 661.58,100.237 681.955,101.672 702.331,103.175 722.707,104.748 743.082,106.389 763.458,108.099 783.833,109.877 804.209,111.723 824.585,113.638 844.96,115.619 865.336,117.669 885.712,119.786 906.087,121.969 926.463,124.22 946.838,126.538 967.214,128.922 987.59,131.372 1007.97,133.888 1028.34,136.47 1048.72,139.117 1069.09,141.83 1089.47,144.607 1109.84,147.449 1130.22,150.355 1150.59,153.326 1170.97,156.36 1191.35,159.458 1211.72,162.619 1232.1,165.843 1252.47,169.129 1272.85,172.477 1293.22,175.888 1313.6,179.359 1333.98,182.893 1354.35,186.486 1374.73,190.141 1395.1,193.855 1415.48,197.629 1435.85,201.463 1456.23,205.355 1476.6,209.306 1496.98,213.315 1517.36,217.382 1537.73,221.506 1558.11,225.688 1578.48,229.926 1598.86,234.219 1619.23,238.569 1639.61,242.974 1659.98,247.434 1680.36,251.948 1700.74,256.516 1721.11,261.137 1741.49,265.811 1761.86,270.538 1782.24,275.317 1802.61,280.148 1822.99,285.03 1843.37,289.962 1863.74,294.944 1884.12,299.976 1904.49,305.057 1924.87,310.187 1945.24,315.365 1965.62,320.59 1985.99,325.862 2006.37,331.181 2026.75,336.546 2047.12,341.957 2067.5,347.412 2087.87,352.912 2108.25,358.455 2128.62,364.042 2149,369.671 2169.38,375.343 2189.75,381.056 2210.13,386.81 2230.5,392.605 2250.88,398.439 2271.25,404.313 2291.63,410.225 "/>
<path clip-path="url(#clip722)" d="M441.824 86.6717 L441.824 118.672 L473.824 118.672 L473.824 86.6717 L441.824 86.6717 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M527.616 81.1805 L527.616 113.18 L559.616 113.18 L559.616 81.1805 L527.616 81.1805 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M613.408 82.6566 L613.408 114.657 L645.408 114.657 L645.408 82.6566 L613.408 82.6566 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M699.2 87.2619 L699.2 119.262 L731.2 119.262 L731.2 87.2619 L699.2 87.2619 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M784.992 93.7597 L784.992 125.76 L816.992 125.76 L816.992 93.7597 L784.992 93.7597 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M870.784 101.71 L870.784 133.71 L902.784 133.71 L902.784 101.71 L870.784 101.71 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M956.576 110.937 L956.576 142.937 L988.576 142.937 L988.576 110.937 L956.576 110.937 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1042.37 121.372 L1042.37 153.372 L1074.37 153.372 L1074.37 121.372 L1042.37 121.372 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1128.16 132.987 L1128.16 164.987 L1160.16 164.987 L1160.16 132.987 L1128.16 132.987 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1213.95 145.74 L1213.95 177.74 L1245.95 177.74 L1245.95 145.74 L1213.95 145.74 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1299.74 159.614 L1299.74 191.614 L1331.74 191.614 L1331.74 159.614 L1299.74 159.614 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1385.54 174.589 L1385.54 206.589 L1417.54 206.589 L1417.54 174.589 L1385.54 174.589 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1471.33 190.636 L1471.33 222.636 L1503.33 222.636 L1503.33 190.636 L1471.33 190.636 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1557.12 207.73 L1557.12 239.73 L1589.12 239.73 L1589.12 207.73 L1557.12 207.73 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1642.91 225.852 L1642.91 257.852 L1674.91 257.852 L1674.91 225.852 L1642.91 225.852 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1728.7 244.989 L1728.7 276.989 L1760.7 276.989 L1760.7 244.989 L1728.7 244.989 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1814.5 265.177 L1814.5 297.177 L1846.5 297.177 L1846.5 265.177 L1814.5 265.177 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1900.29 286.555 L1900.29 318.555 L1932.29 318.555 L1932.29 286.555 L1900.29 286.555 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1986.08 309.58 L1986.08 341.58 L2018.08 341.58 L2018.08 309.58 L1986.08 309.58 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M2071.87 335.642 L2071.87 367.642 L2103.87 367.642 L2103.87 335.642 L2071.87 335.642 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="457.824" cy="103.295" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="543.616" cy="98.0966" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="629.408" cy="99.8696" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="715.2" cy="104.774" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="800.992" cy="111.571" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="886.784" cy="119.823" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="972.576" cy="129.356" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1058.37" cy="140.099" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1144.16" cy="152.016" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1229.95" cy="165.074" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1315.74" cy="179.255" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1401.54" cy="194.534" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1487.33" cy="210.884" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1573.12" cy="228.278" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1658.91" cy="246.693" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1744.7" cy="266.117" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1830.5" cy="286.58" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1916.29" cy="308.211" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="2002.08" cy="331.446" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="2087.87" cy="357.615" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<polyline clip-path="url(#clip722)" style="stroke:#c271d2; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="254.067,791.086 274.443,798.135 294.819,805.183 315.194,812.23 335.57,819.274 355.946,826.316 376.321,833.355 396.697,840.389 417.072,847.418 437.448,854.441 457.824,861.458 478.199,868.469 498.575,875.471 518.95,882.465 539.326,889.45 559.702,896.425 580.077,903.389 600.453,910.342 620.829,917.284 641.204,924.212 661.58,931.128 681.955,938.029 702.331,944.916 722.707,951.787 743.082,958.642 763.458,965.481 783.833,972.302 804.209,979.105 824.585,985.889 844.96,992.653 865.336,999.398 885.712,1006.12 906.087,1012.82 926.463,1019.5 946.838,1026.16 967.214,1032.79 987.59,1039.4 1007.97,1045.99 1028.34,1052.55 1048.72,1059.08 1069.09,1065.59 1089.47,1072.07 1109.84,1078.52 1130.22,1084.94 1150.59,1091.33 1170.97,1097.69 1191.35,1104.02 1211.72,1110.32 1232.1,1116.59 1252.47,1122.83 1272.85,1129.03 1293.22,1135.2 1313.6,1141.34 1333.98,1147.44 1354.35,1153.5 1374.73,1159.53 1395.1,1165.52 1415.48,1171.47 1435.85,1177.39 1456.23,1183.26 1476.6,1189.1 1496.98,1194.9 1517.36,1200.66 1537.73,1206.37 1558.11,1212.05 1578.48,1217.68 1598.86,1223.27 1619.23,1228.82 1639.61,1234.32 1659.98,1239.78 1680.36,1245.19 1700.74,1250.56 1721.11,1255.89 1741.49,1261.16 1761.86,1266.39 1782.24,1271.57 1802.61,1276.71 1822.99,1281.79 1843.37,1286.83 1863.74,1291.81 1884.12,1296.75 1904.49,1301.64 1924.87,1306.47 1945.24,1311.25 1965.62,1315.98 1985.99,1320.66 2006.37,1325.29 2026.75,1329.86 2047.12,1334.38 2067.5,1338.84 2087.87,1343.25 2108.25,1347.61 2128.62,1351.9 2149,1356.15 2169.38,1360.33 2189.75,1364.46 2210.13,1368.53 2230.5,1372.55 2250.88,1376.5 2271.25,1380.4 2291.63,1384.24 "/>
<path clip-path="url(#clip722)" d="M441.824 841.894 L441.824 873.894 L473.824 873.894 L473.824 841.894 L441.824 841.894 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M527.616 871.204 L527.616 903.204 L559.616 903.204 L559.616 871.204 L527.616 871.204 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M613.408 900.72 L613.408 932.72 L645.408 932.72 L645.408 900.72 L613.408 900.72 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M699.2 930.145 L699.2 962.145 L731.2 962.145 L731.2 930.145 L699.2 930.145 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M784.992 959.343 L784.992 991.343 L816.992 991.343 L816.992 959.343 L784.992 959.343 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M870.784 988.231 L870.784 1020.23 L902.784 1020.23 L902.784 988.231 L870.784 988.231 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M956.576 1016.8 L956.576 1048.8 L988.576 1048.8 L988.576 1016.8 L956.576 1016.8 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1042.37 1044.91 L1042.37 1076.91 L1074.37 1076.91 L1074.37 1044.91 L1042.37 1044.91 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1128.16 1072.5 L1128.16 1104.5 L1160.16 1104.5 L1160.16 1072.5 L1128.16 1072.5 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1213.95 1099.56 L1213.95 1131.56 L1245.95 1131.56 L1245.95 1099.56 L1213.95 1099.56 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1299.74 1126.04 L1299.74 1158.04 L1331.74 1158.04 L1331.74 1126.04 L1299.74 1126.04 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1385.54 1151.88 L1385.54 1183.88 L1417.54 1183.88 L1417.54 1151.88 L1385.54 1151.88 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1471.33 1177.05 L1471.33 1209.05 L1503.33 1209.05 L1503.33 1177.05 L1471.33 1177.05 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1557.12 1201.48 L1557.12 1233.48 L1589.12 1233.48 L1589.12 1201.48 L1557.12 1201.48 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1642.91 1225.11 L1642.91 1257.11 L1674.91 1257.11 L1674.91 1225.11 L1642.91 1225.11 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1728.7 1247.89 L1728.7 1279.89 L1760.7 1279.89 L1760.7 1247.89 L1728.7 1247.89 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1814.5 1269.64 L1814.5 1301.64 L1846.5 1301.64 L1846.5 1269.64 L1814.5 1269.64 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1900.29 1290.04 L1900.29 1322.04 L1932.29 1322.04 L1932.29 1290.04 L1900.29 1290.04 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M1986.08 1308.23 L1986.08 1340.23 L2018.08 1340.23 L2018.08 1308.23 L1986.08 1308.23 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip722)" d="M2071.87 1321.87 L2071.87 1353.87 L2103.87 1353.87 L2103.87 1321.87 L2071.87 1321.87 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="457.824" cy="861.754" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="543.616" cy="890.861" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="629.408" cy="920.033" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="715.2" cy="949.063" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="800.992" cy="977.848" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="886.784" cy="1006.32" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="972.576" cy="1034.41" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1058.37" cy="1062.08" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1144.16" cy="1089.25" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1229.95" cy="1115.9" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1315.74" cy="1141.96" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1401.54" cy="1167.4" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1487.33" cy="1192.16" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1573.12" cy="1216.19" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1658.91" cy="1239.45" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1744.7" cy="1261.84" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1830.5" cy="1283.21" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="1916.29" cy="1303.24" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="2002.08" cy="1321.08" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip722)" cx="2087.87" cy="1334.38" r="14.4" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip720)" d="M264.934 1377.32 L780.845 1377.32 L780.845 1014.44 L264.934 1014.44  Z" fill="#ffffff" fill-rule="evenodd" fill-opacity="1"/>
<polyline clip-path="url(#clip720)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="264.934,1377.32 780.845,1377.32 780.845,1014.44 264.934,1014.44 264.934,1377.32 "/>
<polyline clip-path="url(#clip720)" style="stroke:#009af9; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="288.932,1066.28 432.92,1066.28 "/>
<path clip-path="url(#clip720)" d="M473.307 1067.35 Q474.811 1067.86 476.224 1069.53 Q477.659 1071.19 479.094 1074.11 L483.839 1083.56 L478.816 1083.56 L474.395 1074.69 Q472.682 1071.22 471.061 1070.08 Q469.464 1068.95 466.687 1068.95 L461.594 1068.95 L461.594 1083.56 L456.918 1083.56 L456.918 1049 L467.474 1049 Q473.399 1049 476.316 1051.47 Q479.233 1053.95 479.233 1058.95 Q479.233 1062.21 477.705 1064.37 Q476.2 1066.52 473.307 1067.35 M461.594 1052.84 L461.594 1065.11 L467.474 1065.11 Q470.853 1065.11 472.566 1063.56 Q474.302 1061.98 474.302 1058.95 Q474.302 1055.92 472.566 1054.39 Q470.853 1052.84 467.474 1052.84 L461.594 1052.84 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M509.719 1069.53 L509.719 1071.61 L490.135 1071.61 Q490.413 1076.01 492.774 1078.32 Q495.159 1080.62 499.395 1080.62 Q501.848 1080.62 504.14 1080.01 Q506.455 1079.41 508.723 1078.21 L508.723 1082.24 Q506.432 1083.21 504.024 1083.72 Q501.617 1084.23 499.14 1084.23 Q492.936 1084.23 489.302 1080.62 Q485.691 1077 485.691 1070.85 Q485.691 1064.48 489.117 1060.75 Q492.566 1057 498.399 1057 Q503.631 1057 506.663 1060.38 Q509.719 1063.74 509.719 1069.53 M505.459 1068.28 Q505.413 1064.78 503.492 1062.7 Q501.594 1060.62 498.446 1060.62 Q494.881 1060.62 492.728 1062.63 Q490.598 1064.64 490.274 1068.3 L505.459 1068.28 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M514.557 1068.67 L527.033 1068.67 L527.033 1072.47 L514.557 1072.47 L514.557 1068.67 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M544.51 1053.6 L538.168 1070.8 L550.876 1070.8 L544.51 1053.6 M541.871 1049 L547.172 1049 L560.343 1083.56 L555.482 1083.56 L552.334 1074.69 L536.756 1074.69 L533.607 1083.56 L528.677 1083.56 L541.871 1049 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M586.755 1067.91 L586.755 1083.56 L582.496 1083.56 L582.496 1068.05 Q582.496 1064.37 581.061 1062.54 Q579.626 1060.71 576.755 1060.71 Q573.306 1060.71 571.316 1062.91 Q569.325 1065.11 569.325 1068.9 L569.325 1083.56 L565.042 1083.56 L565.042 1057.63 L569.325 1057.63 L569.325 1061.66 Q570.853 1059.32 572.913 1058.16 Q574.996 1057 577.704 1057 Q582.172 1057 584.464 1059.78 Q586.755 1062.54 586.755 1067.91 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M607.033 1070.52 Q601.871 1070.52 599.88 1071.7 Q597.89 1072.88 597.89 1075.73 Q597.89 1078 599.371 1079.34 Q600.876 1080.66 603.445 1080.66 Q606.987 1080.66 609.116 1078.16 Q611.269 1075.64 611.269 1071.47 L611.269 1070.52 L607.033 1070.52 M615.528 1068.76 L615.528 1083.56 L611.269 1083.56 L611.269 1079.62 Q609.811 1081.98 607.635 1083.12 Q605.459 1084.23 602.311 1084.23 Q598.329 1084.23 595.968 1082 Q593.63 1079.76 593.63 1076.01 Q593.63 1071.63 596.547 1069.41 Q599.487 1067.19 605.297 1067.19 L611.269 1067.19 L611.269 1066.77 Q611.269 1063.83 609.325 1062.24 Q607.403 1060.62 603.908 1060.62 Q601.686 1060.62 599.579 1061.15 Q597.473 1061.68 595.528 1062.75 L595.528 1058.81 Q597.866 1057.91 600.065 1057.47 Q602.264 1057 604.348 1057 Q609.973 1057 612.751 1059.92 Q615.528 1062.84 615.528 1068.76 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M624.301 1047.54 L628.561 1047.54 L628.561 1083.56 L624.301 1083.56 L624.301 1047.54 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M648.26 1085.96 Q646.454 1090.59 644.741 1092 Q643.028 1093.42 640.158 1093.42 L636.755 1093.42 L636.755 1089.85 L639.255 1089.85 Q641.014 1089.85 641.987 1089.02 Q642.959 1088.18 644.139 1085.08 L644.903 1083.14 L634.417 1057.63 L638.931 1057.63 L647.033 1077.91 L655.135 1057.63 L659.648 1057.63 L648.26 1085.96 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M669.741 1050.27 L669.741 1057.63 L678.514 1057.63 L678.514 1060.94 L669.741 1060.94 L669.741 1075.01 Q669.741 1078.18 670.597 1079.09 Q671.477 1079.99 674.139 1079.99 L678.514 1079.99 L678.514 1083.56 L674.139 1083.56 Q669.209 1083.56 667.334 1081.73 Q665.459 1079.87 665.459 1075.01 L665.459 1060.94 L662.334 1060.94 L662.334 1057.63 L665.459 1057.63 L665.459 1050.27 L669.741 1050.27 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M684.116 1057.63 L688.375 1057.63 L688.375 1083.56 L684.116 1083.56 L684.116 1057.63 M684.116 1047.54 L688.375 1047.54 L688.375 1052.93 L684.116 1052.93 L684.116 1047.54 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M715.944 1058.62 L715.944 1062.61 Q714.139 1061.61 712.31 1061.12 Q710.505 1060.62 708.653 1060.62 Q704.509 1060.62 702.218 1063.25 Q699.926 1065.87 699.926 1070.62 Q699.926 1075.36 702.218 1078 Q704.509 1080.62 708.653 1080.62 Q710.505 1080.62 712.31 1080.13 Q714.139 1079.62 715.944 1078.62 L715.944 1082.56 Q714.162 1083.39 712.241 1083.81 Q710.343 1084.23 708.19 1084.23 Q702.333 1084.23 698.884 1080.55 Q695.435 1076.87 695.435 1070.62 Q695.435 1064.27 698.907 1060.64 Q702.403 1057 708.468 1057 Q710.435 1057 712.31 1057.42 Q714.185 1057.81 715.944 1058.62 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M735.134 1070.52 Q729.972 1070.52 727.981 1071.7 Q725.991 1072.88 725.991 1075.73 Q725.991 1078 727.472 1079.34 Q728.977 1080.66 731.546 1080.66 Q735.088 1080.66 737.217 1078.16 Q739.37 1075.64 739.37 1071.47 L739.37 1070.52 L735.134 1070.52 M743.629 1068.76 L743.629 1083.56 L739.37 1083.56 L739.37 1079.62 Q737.912 1081.98 735.736 1083.12 Q733.56 1084.23 730.412 1084.23 Q726.43 1084.23 724.069 1082 Q721.731 1079.76 721.731 1076.01 Q721.731 1071.63 724.648 1069.41 Q727.588 1067.19 733.398 1067.19 L739.37 1067.19 L739.37 1066.77 Q739.37 1063.83 737.426 1062.24 Q735.505 1060.62 732.009 1060.62 Q729.787 1060.62 727.68 1061.15 Q725.574 1061.68 723.63 1062.75 L723.63 1058.81 Q725.968 1057.91 728.167 1057.47 Q730.366 1057 732.449 1057 Q738.074 1057 740.852 1059.92 Q743.629 1062.84 743.629 1068.76 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M752.403 1047.54 L756.662 1047.54 L756.662 1083.56 L752.403 1083.56 L752.403 1047.54 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M338.171 1095.36 L338.171 1140.87 L383.682 1140.87 L383.682 1095.36 L338.171 1095.36 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="4.55111"/>
<path clip-path="url(#clip720)" d="M473.307 1119.19 Q474.811 1119.7 476.224 1121.37 Q477.659 1123.03 479.094 1125.95 L483.839 1135.4 L478.816 1135.4 L474.395 1126.53 Q472.682 1123.06 471.061 1121.92 Q469.464 1120.79 466.687 1120.79 L461.594 1120.79 L461.594 1135.4 L456.918 1135.4 L456.918 1100.84 L467.474 1100.84 Q473.399 1100.84 476.316 1103.31 Q479.233 1105.79 479.233 1110.79 Q479.233 1114.05 477.705 1116.21 Q476.2 1118.36 473.307 1119.19 M461.594 1104.68 L461.594 1116.95 L467.474 1116.95 Q470.853 1116.95 472.566 1115.4 Q474.302 1113.82 474.302 1110.79 Q474.302 1107.76 472.566 1106.23 Q470.853 1104.68 467.474 1104.68 L461.594 1104.68 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M509.719 1121.37 L509.719 1123.45 L490.135 1123.45 Q490.413 1127.85 492.774 1130.16 Q495.159 1132.46 499.395 1132.46 Q501.848 1132.46 504.14 1131.85 Q506.455 1131.25 508.723 1130.05 L508.723 1134.08 Q506.432 1135.05 504.024 1135.56 Q501.617 1136.07 499.14 1136.07 Q492.936 1136.07 489.302 1132.46 Q485.691 1128.84 485.691 1122.69 Q485.691 1116.32 489.117 1112.59 Q492.566 1108.84 498.399 1108.84 Q503.631 1108.84 506.663 1112.22 Q509.719 1115.58 509.719 1121.37 M505.459 1120.12 Q505.413 1116.62 503.492 1114.54 Q501.594 1112.46 498.446 1112.46 Q494.881 1112.46 492.728 1114.47 Q490.598 1116.48 490.274 1120.14 L505.459 1120.12 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M514.557 1120.51 L527.033 1120.51 L527.033 1124.31 L514.557 1124.31 L514.557 1120.51 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M534.001 1100.84 L553.862 1100.84 L553.862 1104.77 L538.677 1104.77 L538.677 1114.96 L552.381 1114.96 L552.381 1118.89 L538.677 1118.89 L538.677 1135.4 L534.001 1135.4 L534.001 1100.84 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M558.052 1125.16 L558.052 1109.47 L562.311 1109.47 L562.311 1125 Q562.311 1128.68 563.746 1130.53 Q565.181 1132.36 568.052 1132.36 Q571.501 1132.36 573.492 1130.16 Q575.505 1127.96 575.505 1124.17 L575.505 1109.47 L579.765 1109.47 L579.765 1135.4 L575.505 1135.4 L575.505 1131.41 Q573.954 1133.77 571.894 1134.93 Q569.857 1136.07 567.149 1136.07 Q562.681 1136.07 560.367 1133.29 Q558.052 1130.51 558.052 1125.16 M568.769 1108.84 L568.769 1108.84 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M588.538 1099.38 L592.797 1099.38 L592.797 1135.4 L588.538 1135.4 L588.538 1099.38 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M601.709 1099.38 L605.968 1099.38 L605.968 1135.4 L601.709 1135.4 L601.709 1099.38 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><circle clip-path="url(#clip720)" cx="360.926" cy="1169.96" r="20.48" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="4.55111"/>
<path clip-path="url(#clip720)" d="M473.307 1171.03 Q474.811 1171.54 476.224 1173.21 Q477.659 1174.87 479.094 1177.79 L483.839 1187.24 L478.816 1187.24 L474.395 1178.37 Q472.682 1174.9 471.061 1173.76 Q469.464 1172.63 466.687 1172.63 L461.594 1172.63 L461.594 1187.24 L456.918 1187.24 L456.918 1152.68 L467.474 1152.68 Q473.399 1152.68 476.316 1155.15 Q479.233 1157.63 479.233 1162.63 Q479.233 1165.89 477.705 1168.05 Q476.2 1170.2 473.307 1171.03 M461.594 1156.52 L461.594 1168.79 L467.474 1168.79 Q470.853 1168.79 472.566 1167.24 Q474.302 1165.66 474.302 1162.63 Q474.302 1159.6 472.566 1158.07 Q470.853 1156.52 467.474 1156.52 L461.594 1156.52 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M509.719 1173.21 L509.719 1175.29 L490.135 1175.29 Q490.413 1179.69 492.774 1182 Q495.159 1184.3 499.395 1184.3 Q501.848 1184.3 504.14 1183.69 Q506.455 1183.09 508.723 1181.89 L508.723 1185.92 Q506.432 1186.89 504.024 1187.4 Q501.617 1187.91 499.14 1187.91 Q492.936 1187.91 489.302 1184.3 Q485.691 1180.68 485.691 1174.53 Q485.691 1168.16 489.117 1164.43 Q492.566 1160.68 498.399 1160.68 Q503.631 1160.68 506.663 1164.06 Q509.719 1167.42 509.719 1173.21 M505.459 1171.96 Q505.413 1168.46 503.492 1166.38 Q501.594 1164.3 498.446 1164.3 Q494.881 1164.3 492.728 1166.31 Q490.598 1168.32 490.274 1171.98 L505.459 1171.96 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M514.557 1172.35 L527.033 1172.35 L527.033 1176.15 L514.557 1176.15 L514.557 1172.35 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M534.001 1152.68 L553.862 1152.68 L553.862 1156.61 L538.677 1156.61 L538.677 1166.8 L552.381 1166.8 L552.381 1170.73 L538.677 1170.73 L538.677 1187.24 L534.001 1187.24 L534.001 1152.68 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M561.269 1152.68 L568.237 1152.68 L577.056 1176.19 L585.922 1152.68 L592.89 1152.68 L592.89 1187.24 L588.329 1187.24 L588.329 1156.89 L579.417 1180.59 L574.718 1180.59 L565.806 1156.89 L565.806 1187.24 L561.269 1187.24 L561.269 1152.68 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M602.172 1152.68 L609.139 1152.68 L617.959 1176.19 L626.825 1152.68 L633.792 1152.68 L633.792 1187.24 L629.232 1187.24 L629.232 1156.89 L620.32 1180.59 L615.621 1180.59 L606.709 1156.89 L606.709 1187.24 L602.172 1187.24 L602.172 1152.68 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><polyline clip-path="url(#clip720)" style="stroke:#c271d2; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="288.932,1221.8 432.92,1221.8 "/>
<path clip-path="url(#clip720)" d="M456.918 1204.52 L461.594 1204.52 L461.594 1239.08 L456.918 1239.08 L456.918 1204.52 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M490.899 1218.13 Q492.497 1215.26 494.719 1213.89 Q496.941 1212.52 499.95 1212.52 Q504.001 1212.52 506.2 1215.37 Q508.399 1218.2 508.399 1223.43 L508.399 1239.08 L504.117 1239.08 L504.117 1223.57 Q504.117 1219.84 502.797 1218.03 Q501.478 1216.23 498.77 1216.23 Q495.46 1216.23 493.538 1218.43 Q491.617 1220.63 491.617 1224.42 L491.617 1239.08 L487.335 1239.08 L487.335 1223.57 Q487.335 1219.82 486.015 1218.03 Q484.696 1216.23 481.941 1216.23 Q478.677 1216.23 476.756 1218.45 Q474.835 1220.65 474.835 1224.42 L474.835 1239.08 L470.552 1239.08 L470.552 1213.15 L474.835 1213.15 L474.835 1217.18 Q476.293 1214.79 478.33 1213.66 Q480.367 1212.52 483.168 1212.52 Q485.992 1212.52 487.96 1213.96 Q489.95 1215.39 490.899 1218.13 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M514.742 1224.19 L527.219 1224.19 L527.219 1227.99 L514.742 1227.99 L514.742 1224.19 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M544.695 1209.12 L538.353 1226.32 L551.061 1226.32 L544.695 1209.12 M542.056 1204.52 L547.357 1204.52 L560.529 1239.08 L555.668 1239.08 L552.519 1230.21 L536.941 1230.21 L533.793 1239.08 L528.862 1239.08 L542.056 1204.52 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M586.941 1223.43 L586.941 1239.08 L582.681 1239.08 L582.681 1223.57 Q582.681 1219.89 581.246 1218.06 Q579.811 1216.23 576.941 1216.23 Q573.492 1216.23 571.501 1218.43 Q569.51 1220.63 569.51 1224.42 L569.51 1239.08 L565.228 1239.08 L565.228 1213.15 L569.51 1213.15 L569.51 1217.18 Q571.038 1214.84 573.098 1213.68 Q575.181 1212.52 577.89 1212.52 Q582.357 1212.52 584.649 1215.3 Q586.941 1218.06 586.941 1223.43 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M607.218 1226.04 Q602.056 1226.04 600.065 1227.22 Q598.075 1228.4 598.075 1231.25 Q598.075 1233.52 599.556 1234.86 Q601.061 1236.18 603.63 1236.18 Q607.172 1236.18 609.301 1233.68 Q611.454 1231.16 611.454 1226.99 L611.454 1226.04 L607.218 1226.04 M615.714 1224.28 L615.714 1239.08 L611.454 1239.08 L611.454 1235.14 Q609.996 1237.5 607.82 1238.64 Q605.644 1239.75 602.496 1239.75 Q598.515 1239.75 596.153 1237.52 Q593.815 1235.28 593.815 1231.53 Q593.815 1227.15 596.732 1224.93 Q599.672 1222.71 605.482 1222.71 L611.454 1222.71 L611.454 1222.29 Q611.454 1219.35 609.51 1217.76 Q607.589 1216.14 604.093 1216.14 Q601.871 1216.14 599.765 1216.67 Q597.658 1217.2 595.714 1218.27 L595.714 1214.33 Q598.052 1213.43 600.251 1212.99 Q602.45 1212.52 604.533 1212.52 Q610.158 1212.52 612.936 1215.44 Q615.714 1218.36 615.714 1224.28 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M624.487 1203.06 L628.746 1203.06 L628.746 1239.08 L624.487 1239.08 L624.487 1203.06 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M648.445 1241.48 Q646.639 1246.11 644.926 1247.52 Q643.213 1248.94 640.343 1248.94 L636.94 1248.94 L636.94 1245.37 L639.44 1245.37 Q641.199 1245.37 642.172 1244.54 Q643.144 1243.7 644.324 1240.6 L645.088 1238.66 L634.602 1213.15 L639.116 1213.15 L647.218 1233.43 L655.32 1213.15 L659.834 1213.15 L648.445 1241.48 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M669.926 1205.79 L669.926 1213.15 L678.699 1213.15 L678.699 1216.46 L669.926 1216.46 L669.926 1230.53 Q669.926 1233.7 670.783 1234.61 Q671.662 1235.51 674.324 1235.51 L678.699 1235.51 L678.699 1239.08 L674.324 1239.08 Q669.394 1239.08 667.519 1237.25 Q665.644 1235.39 665.644 1230.53 L665.644 1216.46 L662.519 1216.46 L662.519 1213.15 L665.644 1213.15 L665.644 1205.79 L669.926 1205.79 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M684.301 1213.15 L688.56 1213.15 L688.56 1239.08 L684.301 1239.08 L684.301 1213.15 M684.301 1203.06 L688.56 1203.06 L688.56 1208.45 L684.301 1208.45 L684.301 1203.06 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M716.13 1214.14 L716.13 1218.13 Q714.324 1217.13 712.495 1216.64 Q710.69 1216.14 708.838 1216.14 Q704.695 1216.14 702.403 1218.77 Q700.111 1221.39 700.111 1226.14 Q700.111 1230.88 702.403 1233.52 Q704.695 1236.14 708.838 1236.14 Q710.69 1236.14 712.495 1235.65 Q714.324 1235.14 716.13 1234.14 L716.13 1238.08 Q714.347 1238.91 712.426 1239.33 Q710.528 1239.75 708.375 1239.75 Q702.519 1239.75 699.07 1236.07 Q695.62 1232.39 695.62 1226.14 Q695.62 1219.79 699.093 1216.16 Q702.588 1212.52 708.653 1212.52 Q710.62 1212.52 712.495 1212.94 Q714.37 1213.33 716.13 1214.14 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M735.319 1226.04 Q730.157 1226.04 728.167 1227.22 Q726.176 1228.4 726.176 1231.25 Q726.176 1233.52 727.657 1234.86 Q729.162 1236.18 731.731 1236.18 Q735.273 1236.18 737.403 1233.68 Q739.555 1231.16 739.555 1226.99 L739.555 1226.04 L735.319 1226.04 M743.815 1224.28 L743.815 1239.08 L739.555 1239.08 L739.555 1235.14 Q738.097 1237.5 735.921 1238.64 Q733.745 1239.75 730.597 1239.75 Q726.616 1239.75 724.255 1237.52 Q721.917 1235.28 721.917 1231.53 Q721.917 1227.15 724.833 1224.93 Q727.773 1222.71 733.583 1222.71 L739.555 1222.71 L739.555 1222.29 Q739.555 1219.35 737.611 1217.76 Q735.69 1216.14 732.194 1216.14 Q729.972 1216.14 727.866 1216.67 Q725.759 1217.2 723.815 1218.27 L723.815 1214.33 Q726.153 1213.43 728.352 1212.99 Q730.551 1212.52 732.634 1212.52 Q738.259 1212.52 741.037 1215.44 Q743.815 1218.36 743.815 1224.28 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M752.588 1203.06 L756.847 1203.06 L756.847 1239.08 L752.588 1239.08 L752.588 1203.06 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M338.171 1250.88 L338.171 1296.39 L383.682 1296.39 L383.682 1250.88 L338.171 1250.88 Z" fill="#ac8d18" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="4.55111"/>
<path clip-path="url(#clip720)" d="M456.918 1256.36 L461.594 1256.36 L461.594 1290.92 L456.918 1290.92 L456.918 1256.36 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M490.899 1269.97 Q492.497 1267.1 494.719 1265.73 Q496.941 1264.36 499.95 1264.36 Q504.001 1264.36 506.2 1267.21 Q508.399 1270.04 508.399 1275.27 L508.399 1290.92 L504.117 1290.92 L504.117 1275.41 Q504.117 1271.68 502.797 1269.87 Q501.478 1268.07 498.77 1268.07 Q495.46 1268.07 493.538 1270.27 Q491.617 1272.47 491.617 1276.26 L491.617 1290.92 L487.335 1290.92 L487.335 1275.41 Q487.335 1271.66 486.015 1269.87 Q484.696 1268.07 481.941 1268.07 Q478.677 1268.07 476.756 1270.29 Q474.835 1272.49 474.835 1276.26 L474.835 1290.92 L470.552 1290.92 L470.552 1264.99 L474.835 1264.99 L474.835 1269.02 Q476.293 1266.63 478.33 1265.5 Q480.367 1264.36 483.168 1264.36 Q485.992 1264.36 487.96 1265.8 Q489.95 1267.23 490.899 1269.97 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M514.742 1276.03 L527.219 1276.03 L527.219 1279.83 L514.742 1279.83 L514.742 1276.03 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M534.186 1256.36 L554.047 1256.36 L554.047 1260.29 L538.862 1260.29 L538.862 1270.48 L552.566 1270.48 L552.566 1274.41 L538.862 1274.41 L538.862 1290.92 L534.186 1290.92 L534.186 1256.36 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M558.237 1280.68 L558.237 1264.99 L562.496 1264.99 L562.496 1280.52 Q562.496 1284.2 563.931 1286.05 Q565.367 1287.88 568.237 1287.88 Q571.686 1287.88 573.677 1285.68 Q575.691 1283.48 575.691 1279.69 L575.691 1264.99 L579.95 1264.99 L579.95 1290.92 L575.691 1290.92 L575.691 1286.93 Q574.14 1289.29 572.079 1290.45 Q570.042 1291.59 567.334 1291.59 Q562.867 1291.59 560.552 1288.81 Q558.237 1286.03 558.237 1280.68 M568.954 1264.36 L568.954 1264.36 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M588.723 1254.9 L592.982 1254.9 L592.982 1290.92 L588.723 1290.92 L588.723 1254.9 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M601.894 1254.9 L606.153 1254.9 L606.153 1290.92 L601.894 1290.92 L601.894 1254.9 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><circle clip-path="url(#clip720)" cx="360.926" cy="1325.48" r="20.48" fill="#00a9ad" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="4.55111"/>
<path clip-path="url(#clip720)" d="M456.918 1308.2 L461.594 1308.2 L461.594 1342.76 L456.918 1342.76 L456.918 1308.2 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M490.899 1321.81 Q492.497 1318.94 494.719 1317.57 Q496.941 1316.2 499.95 1316.2 Q504.001 1316.2 506.2 1319.05 Q508.399 1321.88 508.399 1327.11 L508.399 1342.76 L504.117 1342.76 L504.117 1327.25 Q504.117 1323.52 502.797 1321.71 Q501.478 1319.91 498.77 1319.91 Q495.46 1319.91 493.538 1322.11 Q491.617 1324.31 491.617 1328.1 L491.617 1342.76 L487.335 1342.76 L487.335 1327.25 Q487.335 1323.5 486.015 1321.71 Q484.696 1319.91 481.941 1319.91 Q478.677 1319.91 476.756 1322.13 Q474.835 1324.33 474.835 1328.1 L474.835 1342.76 L470.552 1342.76 L470.552 1316.83 L474.835 1316.83 L474.835 1320.86 Q476.293 1318.47 478.33 1317.34 Q480.367 1316.2 483.168 1316.2 Q485.992 1316.2 487.96 1317.64 Q489.95 1319.07 490.899 1321.81 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M514.742 1327.87 L527.219 1327.87 L527.219 1331.67 L514.742 1331.67 L514.742 1327.87 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M534.186 1308.2 L554.047 1308.2 L554.047 1312.13 L538.862 1312.13 L538.862 1322.32 L552.566 1322.32 L552.566 1326.25 L538.862 1326.25 L538.862 1342.76 L534.186 1342.76 L534.186 1308.2 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M561.455 1308.2 L568.422 1308.2 L577.241 1331.71 L586.107 1308.2 L593.075 1308.2 L593.075 1342.76 L588.515 1342.76 L588.515 1312.41 L579.603 1336.11 L574.904 1336.11 L565.992 1312.41 L565.992 1342.76 L561.455 1342.76 L561.455 1308.2 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip720)" d="M602.357 1308.2 L609.325 1308.2 L618.144 1331.71 L627.01 1308.2 L633.977 1308.2 L633.977 1342.76 L629.417 1342.76 L629.417 1312.41 L620.505 1336.11 L615.806 1336.11 L606.894 1312.41 L606.894 1342.76 L602.357 1342.76 L602.357 1308.2 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /></svg>

```

The surface pressures can also be plotted. Note that it is only possible to plot linear meshes, meaing that we must remove the quadratic parts.

````julia
data_mesh,data_viz = create_vizualization_data(mesh,p_fmm)
fig, ax, hm = viz(data_mesh;showfacets=true, color=real.(data_viz/P₀))
wgl.Colorbar(fig[1,2],label="Re(p/p₀)");
````

![](3d_cube_wave_anechoic.png)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

