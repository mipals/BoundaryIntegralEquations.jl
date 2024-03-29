```@meta
EditURL = "../../../examples/3d_oscilating_sphere.jl"
```

# Oscilating sphere (Interior)
# Importing related packages

````julia
using BoundaryIntegralEquations # For BIEs
using IterativeSolvers          # For gmres
using SpecialFunctions          # For Bessel functions
using LinearAlgebra             # For Diagonal
using Meshes                    # For 3d mesh plots
using Plots                     # For 2d plots
import GLMakie as wgl
wgl.set_theme!(resolution=(1600, 1600))
````

# Setting up constants

````julia
frequency = 100.0;      # Frequency                [Hz]
c  = 343.0;             # Speed up sound           [m/s]
ρ₀ = 1.21;              # Mean density             [kg/m^3]
Z₀ = ρ₀*c;              # Characteristic impedance [Rayl]
vz = 1.0;               # Speed in the z-direction [m/s]
a  = 1.0;               # Radius of sphere_1m      [m]
k  = 2π*frequency/c;    # Wavenumber
````

# Loading Mesh
Loading and visualizing the triangular (spherical) mesh

````julia
mesh_file = joinpath(dirname(pathof(BoundaryIntegralEquations)),"..","examples","meshes","sphere_1m_coarse");
mesh = load3dTriangularComsolMesh(mesh_file)
````

````
Number of elements: 	 464
Number of unkowns:  	 930
Geometry defined by:	 BoundaryIntegralEquations.TriangularQuadratic{Float64}
Physics defined by: 	 BoundaryIntegralEquations.TriangularQuadratic{Float64}

````

# Analytical Solution
The analytical description of the interior pressure of an z-oscilating sphere is given by
```math
 p_\text{analytical}(x,y,z) = -3v_z z \mathrm{i}kZ_0\frac{j_1(kr)}{kr(j_0(ka) - 2j_2(ka))},
```
where ``j_i`` is the spherical bessel function of the first kind, ``r`` is the distance from origo to the point ``(x,y,z)`` and ``v_z`` is the oscilating velocity in the ``z``-direction. For the points in the interior where the pressure will evaluated ``z`` will be equal to ``r`` as the remaining two coordinates are set to 0.
From this the analytical expression is computed

````julia
z_ana = collect(0.0:0.01:1)
p_analytical = -3*vz*z_ana*im*k*Z₀.*(sphericalbesselj.(1,k*z_ana))./(k*z_ana*(sphericalbesselj.(0,k*a) - 2*sphericalbesselj.(2,k*a)));
````

# Solution using the BEM
We start by solving the BEM system using dense matrices. For this we need to first assemble the matrices

````julia
F,G,C = assemble_parallel!(mesh,k,mesh.sources,n=3,m=3,progress=false);
H  = Diagonal(1.0 .- C) - F; # Adding the integral free term
````

We now setup the right-hand side

````julia
vs = im*Z₀*k*mesh.normals[3,:]*vz;
b = G*vs; # Computing right-hand side
````

Using this we can solve for the surface pressure

````julia
p_bem = gmres(H,b;verbose=true);
````

````
=== gmres ===
rest	iter	resnorm
  1	  1	1.23e+01
  1	  2	1.07e+00
  1	  3	5.22e-01
  1	  4	5.08e-02
  1	  5	1.88e-02
  1	  6	5.39e-03
  1	  7	1.62e-04
  1	  8	2.77e-06


````

Similarly we can compute the BEM solution using the Fast Multipole Operators

````julia
Gf = FMMGOperator(mesh,k);
Ff = FMMFOperator(mesh,k);
Hf = Diagonal(1.0 .- C) - Ff;
bf = Gf*vs;                        # Computing right-hand side
p_fmm = gmres(Hf,bf;verbose=true); # Solving the linear system
````

````
=== gmres ===
rest	iter	resnorm
  1	  1	1.29e+01
  1	  2	1.13e+00
  1	  3	5.38e-01
  1	  4	5.32e-02
  1	  5	2.15e-02
  1	  6	5.42e-03
  1	  7	1.52e-04
  1	  8	2.55e-06


````

# Field evaluations
In order to compare the results with the analytical solution we must use the found surface pressures to compute pressure at the interior field points. We therefore start by creating a matrix of ``N`` points (as columns) in 3-dimensional space

````julia
N = 20;
x = zeros(N);
y = zeros(N);
z = collect(0.1:(0.9-0.1)/(N-1):0.9);
X_fieldpoints = [x'; y'; z'];
````

Using the described field-points we can assemble the (dense) matrices for the field point evaluation

````julia
Fp,Gp,_ = assemble_parallel!(mesh,k,X_fieldpoints,n=2,m=2,progress=false);
````

Using the matrices we can evalute the pressure at the field points using the previously found surface pressure

````julia
p_field  = Fp*p_bem + Gp*vs;
pf_field = Fp*p_fmm + Gp*vs;
````

Plotting the field point pressure it can be seen that all 3 methods find the correct pressures

````julia
plot(z_ana,abs.(p_analytical./Z₀),  label="Analytical")
scatter!(z,abs.(p_field./Z₀),       label="Full", markershape=:rect)
scatter!(z,abs.(pf_field./Z₀),      label="FMM")
ylabel!("p/Z₀"); xlabel!("r/a")
````

```@raw html
<?xml version="1.0" encoding="utf-8"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="600" height="400" viewBox="0 0 2400 1600">
<defs>
  <clipPath id="clip600">
    <rect x="0" y="0" width="2400" height="1600"/>
  </clipPath>
</defs>
<path clip-path="url(#clip600)" d="M0 1600 L2400 1600 L2400 0 L0 0  Z" fill="#ffffff" fill-rule="evenodd" fill-opacity="1"/>
<defs>
  <clipPath id="clip601">
    <rect x="480" y="0" width="1681" height="1600"/>
  </clipPath>
</defs>
<path clip-path="url(#clip600)" d="M175.445 1423.18 L2352.76 1423.18 L2352.76 47.2441 L175.445 47.2441  Z" fill="#ffffff" fill-rule="evenodd" fill-opacity="1"/>
<defs>
  <clipPath id="clip602">
    <rect x="175" y="47" width="2178" height="1377"/>
  </clipPath>
</defs>
<polyline clip-path="url(#clip602)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="237.067,1423.18 237.067,47.2441 "/>
<polyline clip-path="url(#clip602)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="750.584,1423.18 750.584,47.2441 "/>
<polyline clip-path="url(#clip602)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="1264.1,1423.18 1264.1,47.2441 "/>
<polyline clip-path="url(#clip602)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="1777.62,1423.18 1777.62,47.2441 "/>
<polyline clip-path="url(#clip602)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="2291.13,1423.18 2291.13,47.2441 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="175.445,1423.18 2352.76,1423.18 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="237.067,1423.18 237.067,1404.28 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="750.584,1423.18 750.584,1404.28 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="1264.1,1423.18 1264.1,1404.28 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="1777.62,1423.18 1777.62,1404.28 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="2291.13,1423.18 2291.13,1404.28 "/>
<path clip-path="url(#clip600)" d="M199.371 1454.1 Q195.76 1454.1 193.931 1457.66 Q192.125 1461.2 192.125 1468.33 Q192.125 1475.44 193.931 1479.01 Q195.76 1482.55 199.371 1482.55 Q203.005 1482.55 204.81 1479.01 Q206.639 1475.44 206.639 1468.33 Q206.639 1461.2 204.81 1457.66 Q203.005 1454.1 199.371 1454.1 M199.371 1450.39 Q205.181 1450.39 208.236 1455 Q211.315 1459.58 211.315 1468.33 Q211.315 1477.06 208.236 1481.67 Q205.181 1486.25 199.371 1486.25 Q193.56 1486.25 190.482 1481.67 Q187.426 1477.06 187.426 1468.33 Q187.426 1459.58 190.482 1455 Q193.56 1450.39 199.371 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M219.533 1479.7 L224.417 1479.7 L224.417 1485.58 L219.533 1485.58 L219.533 1479.7 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M244.602 1454.1 Q240.991 1454.1 239.162 1457.66 Q237.357 1461.2 237.357 1468.33 Q237.357 1475.44 239.162 1479.01 Q240.991 1482.55 244.602 1482.55 Q248.236 1482.55 250.042 1479.01 Q251.87 1475.44 251.87 1468.33 Q251.87 1461.2 250.042 1457.66 Q248.236 1454.1 244.602 1454.1 M244.602 1450.39 Q250.412 1450.39 253.468 1455 Q256.546 1459.58 256.546 1468.33 Q256.546 1477.06 253.468 1481.67 Q250.412 1486.25 244.602 1486.25 Q238.792 1486.25 235.713 1481.67 Q232.658 1477.06 232.658 1468.33 Q232.658 1459.58 235.713 1455 Q238.792 1450.39 244.602 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M274.764 1454.1 Q271.153 1454.1 269.324 1457.66 Q267.518 1461.2 267.518 1468.33 Q267.518 1475.44 269.324 1479.01 Q271.153 1482.55 274.764 1482.55 Q278.398 1482.55 280.204 1479.01 Q282.032 1475.44 282.032 1468.33 Q282.032 1461.2 280.204 1457.66 Q278.398 1454.1 274.764 1454.1 M274.764 1450.39 Q280.574 1450.39 283.629 1455 Q286.708 1459.58 286.708 1468.33 Q286.708 1477.06 283.629 1481.67 Q280.574 1486.25 274.764 1486.25 Q268.954 1486.25 265.875 1481.67 Q262.819 1477.06 262.819 1468.33 Q262.819 1459.58 265.875 1455 Q268.954 1450.39 274.764 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M713.385 1454.1 Q709.774 1454.1 707.945 1457.66 Q706.14 1461.2 706.14 1468.33 Q706.14 1475.44 707.945 1479.01 Q709.774 1482.55 713.385 1482.55 Q717.019 1482.55 718.825 1479.01 Q720.653 1475.44 720.653 1468.33 Q720.653 1461.2 718.825 1457.66 Q717.019 1454.1 713.385 1454.1 M713.385 1450.39 Q719.195 1450.39 722.251 1455 Q725.329 1459.58 725.329 1468.33 Q725.329 1477.06 722.251 1481.67 Q719.195 1486.25 713.385 1486.25 Q707.575 1486.25 704.496 1481.67 Q701.441 1477.06 701.441 1468.33 Q701.441 1459.58 704.496 1455 Q707.575 1450.39 713.385 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M733.547 1479.7 L738.431 1479.7 L738.431 1485.58 L733.547 1485.58 L733.547 1479.7 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M752.644 1481.64 L768.963 1481.64 L768.963 1485.58 L747.019 1485.58 L747.019 1481.64 Q749.681 1478.89 754.264 1474.26 Q758.871 1469.61 760.051 1468.27 Q762.297 1465.74 763.176 1464.01 Q764.079 1462.25 764.079 1460.56 Q764.079 1457.8 762.135 1456.07 Q760.213 1454.33 757.112 1454.33 Q754.913 1454.33 752.459 1455.09 Q750.028 1455.86 747.251 1457.41 L747.251 1452.69 Q750.075 1451.55 752.528 1450.97 Q754.982 1450.39 757.019 1450.39 Q762.389 1450.39 765.584 1453.08 Q768.778 1455.77 768.778 1460.26 Q768.778 1462.39 767.968 1464.31 Q767.181 1466.2 765.075 1468.8 Q764.496 1469.47 761.394 1472.69 Q758.292 1475.88 752.644 1481.64 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M778.824 1451.02 L797.181 1451.02 L797.181 1454.96 L783.107 1454.96 L783.107 1463.43 Q784.125 1463.08 785.144 1462.92 Q786.162 1462.73 787.181 1462.73 Q792.968 1462.73 796.348 1465.9 Q799.727 1469.08 799.727 1474.49 Q799.727 1480.07 796.255 1483.17 Q792.783 1486.25 786.463 1486.25 Q784.287 1486.25 782.019 1485.88 Q779.774 1485.51 777.366 1484.77 L777.366 1480.07 Q779.449 1481.2 781.672 1481.76 Q783.894 1482.32 786.371 1482.32 Q790.375 1482.32 792.713 1480.21 Q795.051 1478.1 795.051 1474.49 Q795.051 1470.88 792.713 1468.77 Q790.375 1466.67 786.371 1466.67 Q784.496 1466.67 782.621 1467.08 Q780.769 1467.5 778.824 1468.38 L778.824 1451.02 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1226.4 1454.1 Q1222.79 1454.1 1220.96 1457.66 Q1219.16 1461.2 1219.16 1468.33 Q1219.16 1475.44 1220.96 1479.01 Q1222.79 1482.55 1226.4 1482.55 Q1230.04 1482.55 1231.84 1479.01 Q1233.67 1475.44 1233.67 1468.33 Q1233.67 1461.2 1231.84 1457.66 Q1230.04 1454.1 1226.4 1454.1 M1226.4 1450.39 Q1232.21 1450.39 1235.27 1455 Q1238.35 1459.58 1238.35 1468.33 Q1238.35 1477.06 1235.27 1481.67 Q1232.21 1486.25 1226.4 1486.25 Q1220.59 1486.25 1217.52 1481.67 Q1214.46 1477.06 1214.46 1468.33 Q1214.46 1459.58 1217.52 1455 Q1220.59 1450.39 1226.4 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1246.57 1479.7 L1251.45 1479.7 L1251.45 1485.58 L1246.57 1485.58 L1246.57 1479.7 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1261.68 1451.02 L1280.04 1451.02 L1280.04 1454.96 L1265.96 1454.96 L1265.96 1463.43 Q1266.98 1463.08 1268 1462.92 Q1269.02 1462.73 1270.04 1462.73 Q1275.83 1462.73 1279.2 1465.9 Q1282.58 1469.08 1282.58 1474.49 Q1282.58 1480.07 1279.11 1483.17 Q1275.64 1486.25 1269.32 1486.25 Q1267.14 1486.25 1264.88 1485.88 Q1262.63 1485.51 1260.22 1484.77 L1260.22 1480.07 Q1262.31 1481.2 1264.53 1481.76 Q1266.75 1482.32 1269.23 1482.32 Q1273.23 1482.32 1275.57 1480.21 Q1277.91 1478.1 1277.91 1474.49 Q1277.91 1470.88 1275.57 1468.77 Q1273.23 1466.67 1269.23 1466.67 Q1267.35 1466.67 1265.48 1467.08 Q1263.63 1467.5 1261.68 1468.38 L1261.68 1451.02 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1301.8 1454.1 Q1298.19 1454.1 1296.36 1457.66 Q1294.55 1461.2 1294.55 1468.33 Q1294.55 1475.44 1296.36 1479.01 Q1298.19 1482.55 1301.8 1482.55 Q1305.43 1482.55 1307.24 1479.01 Q1309.07 1475.44 1309.07 1468.33 Q1309.07 1461.2 1307.24 1457.66 Q1305.43 1454.1 1301.8 1454.1 M1301.8 1450.39 Q1307.61 1450.39 1310.66 1455 Q1313.74 1459.58 1313.74 1468.33 Q1313.74 1477.06 1310.66 1481.67 Q1307.61 1486.25 1301.8 1486.25 Q1295.99 1486.25 1292.91 1481.67 Q1289.85 1477.06 1289.85 1468.33 Q1289.85 1459.58 1292.91 1455 Q1295.99 1450.39 1301.8 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1740.42 1454.1 Q1736.81 1454.1 1734.98 1457.66 Q1733.17 1461.2 1733.17 1468.33 Q1733.17 1475.44 1734.98 1479.01 Q1736.81 1482.55 1740.42 1482.55 Q1744.05 1482.55 1745.86 1479.01 Q1747.69 1475.44 1747.69 1468.33 Q1747.69 1461.2 1745.86 1457.66 Q1744.05 1454.1 1740.42 1454.1 M1740.42 1450.39 Q1746.23 1450.39 1749.28 1455 Q1752.36 1459.58 1752.36 1468.33 Q1752.36 1477.06 1749.28 1481.67 Q1746.23 1486.25 1740.42 1486.25 Q1734.61 1486.25 1731.53 1481.67 Q1728.47 1477.06 1728.47 1468.33 Q1728.47 1459.58 1731.53 1455 Q1734.61 1450.39 1740.42 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1760.58 1479.7 L1765.46 1479.7 L1765.46 1485.58 L1760.58 1485.58 L1760.58 1479.7 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1774.47 1451.02 L1796.69 1451.02 L1796.69 1453.01 L1784.14 1485.58 L1779.26 1485.58 L1791.07 1454.96 L1774.47 1454.96 L1774.47 1451.02 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1805.86 1451.02 L1824.21 1451.02 L1824.21 1454.96 L1810.14 1454.96 L1810.14 1463.43 Q1811.16 1463.08 1812.18 1462.92 Q1813.2 1462.73 1814.21 1462.73 Q1820 1462.73 1823.38 1465.9 Q1826.76 1469.08 1826.76 1474.49 Q1826.76 1480.07 1823.29 1483.17 Q1819.82 1486.25 1813.5 1486.25 Q1811.32 1486.25 1809.05 1485.88 Q1806.81 1485.51 1804.4 1484.77 L1804.4 1480.07 Q1806.48 1481.2 1808.71 1481.76 Q1810.93 1482.32 1813.4 1482.32 Q1817.41 1482.32 1819.75 1480.21 Q1822.08 1478.1 1822.08 1474.49 Q1822.08 1470.88 1819.75 1468.77 Q1817.41 1466.67 1813.4 1466.67 Q1811.53 1466.67 1809.65 1467.08 Q1807.8 1467.5 1805.86 1468.38 L1805.86 1451.02 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2243.21 1481.64 L2250.84 1481.64 L2250.84 1455.28 L2242.53 1456.95 L2242.53 1452.69 L2250.8 1451.02 L2255.47 1451.02 L2255.47 1481.64 L2263.11 1481.64 L2263.11 1485.58 L2243.21 1485.58 L2243.21 1481.64 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2272.56 1479.7 L2277.44 1479.7 L2277.44 1485.58 L2272.56 1485.58 L2272.56 1479.7 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2297.63 1454.1 Q2294.02 1454.1 2292.19 1457.66 Q2290.38 1461.2 2290.38 1468.33 Q2290.38 1475.44 2292.19 1479.01 Q2294.02 1482.55 2297.63 1482.55 Q2301.26 1482.55 2303.07 1479.01 Q2304.9 1475.44 2304.9 1468.33 Q2304.9 1461.2 2303.07 1457.66 Q2301.26 1454.1 2297.63 1454.1 M2297.63 1450.39 Q2303.44 1450.39 2306.49 1455 Q2309.57 1459.58 2309.57 1468.33 Q2309.57 1477.06 2306.49 1481.67 Q2303.44 1486.25 2297.63 1486.25 Q2291.82 1486.25 2288.74 1481.67 Q2285.68 1477.06 2285.68 1468.33 Q2285.68 1459.58 2288.74 1455 Q2291.82 1450.39 2297.63 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2327.79 1454.1 Q2324.18 1454.1 2322.35 1457.66 Q2320.54 1461.2 2320.54 1468.33 Q2320.54 1475.44 2322.35 1479.01 Q2324.18 1482.55 2327.79 1482.55 Q2331.42 1482.55 2333.23 1479.01 Q2335.06 1475.44 2335.06 1468.33 Q2335.06 1461.2 2333.23 1457.66 Q2331.42 1454.1 2327.79 1454.1 M2327.79 1450.39 Q2333.6 1450.39 2336.65 1455 Q2339.73 1459.58 2339.73 1468.33 Q2339.73 1477.06 2336.65 1481.67 Q2333.6 1486.25 2327.79 1486.25 Q2321.98 1486.25 2318.9 1481.67 Q2315.84 1477.06 2315.84 1468.33 Q2315.84 1459.58 2318.9 1455 Q2321.98 1450.39 2327.79 1450.39 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1246.55 1537.87 Q1245.56 1537.3 1244.38 1537.04 Q1243.24 1536.76 1241.84 1536.76 Q1236.87 1536.76 1234.2 1540 Q1231.56 1543.22 1231.56 1549.27 L1231.56 1568.04 L1225.67 1568.04 L1225.67 1532.4 L1231.56 1532.4 L1231.56 1537.93 Q1233.4 1534.69 1236.36 1533.13 Q1239.32 1531.54 1243.56 1531.54 Q1244.16 1531.54 1244.89 1531.63 Q1245.62 1531.7 1246.52 1531.85 L1246.55 1537.87 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1263.1 1520.52 L1268.51 1520.52 L1251.96 1574.09 L1246.55 1574.09 L1263.1 1520.52 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1290.85 1550.12 Q1283.75 1550.12 1281.02 1551.75 Q1278.28 1553.37 1278.28 1557.29 Q1278.28 1560.4 1280.32 1562.25 Q1282.39 1564.07 1285.92 1564.07 Q1290.79 1564.07 1293.72 1560.63 Q1296.68 1557.16 1296.68 1551.43 L1296.68 1550.12 L1290.85 1550.12 M1302.53 1547.71 L1302.53 1568.04 L1296.68 1568.04 L1296.68 1562.63 Q1294.67 1565.88 1291.68 1567.44 Q1288.69 1568.97 1284.36 1568.97 Q1278.88 1568.97 1275.64 1565.91 Q1272.42 1562.82 1272.42 1557.67 Q1272.42 1551.65 1276.43 1548.6 Q1280.48 1545.54 1288.47 1545.54 L1296.68 1545.54 L1296.68 1544.97 Q1296.68 1540.93 1294 1538.73 Q1291.36 1536.5 1286.56 1536.5 Q1283.5 1536.5 1280.6 1537.23 Q1277.71 1537.97 1275.03 1539.43 L1275.03 1534.02 Q1278.25 1532.78 1281.27 1532.17 Q1284.3 1531.54 1287.16 1531.54 Q1294.89 1531.54 1298.71 1535.55 Q1302.53 1539.56 1302.53 1547.71 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><polyline clip-path="url(#clip602)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="175.445,1402.99 2352.76,1402.99 "/>
<polyline clip-path="url(#clip602)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="175.445,1039.13 2352.76,1039.13 "/>
<polyline clip-path="url(#clip602)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="175.445,675.272 2352.76,675.272 "/>
<polyline clip-path="url(#clip602)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="175.445,311.413 2352.76,311.413 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="175.445,1423.18 175.445,47.2441 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="175.445,1402.99 194.343,1402.99 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="175.445,1039.13 194.343,1039.13 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="175.445,675.272 194.343,675.272 "/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="175.445,311.413 194.343,311.413 "/>
<path clip-path="url(#clip600)" d="M127.501 1388.79 Q123.89 1388.79 122.061 1392.35 Q120.255 1395.89 120.255 1403.02 Q120.255 1410.13 122.061 1413.69 Q123.89 1417.24 127.501 1417.24 Q131.135 1417.24 132.941 1413.69 Q134.769 1410.13 134.769 1403.02 Q134.769 1395.89 132.941 1392.35 Q131.135 1388.79 127.501 1388.79 M127.501 1385.08 Q133.311 1385.08 136.367 1389.69 Q139.445 1394.27 139.445 1403.02 Q139.445 1411.75 136.367 1416.36 Q133.311 1420.94 127.501 1420.94 Q121.691 1420.94 118.612 1416.36 Q115.556 1411.75 115.556 1403.02 Q115.556 1394.27 118.612 1389.69 Q121.691 1385.08 127.501 1385.08 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M123.126 1052.48 L139.445 1052.48 L139.445 1056.41 L117.501 1056.41 L117.501 1052.48 Q120.163 1049.72 124.746 1045.09 Q129.353 1040.44 130.533 1039.1 Q132.779 1036.57 133.658 1034.84 Q134.561 1033.08 134.561 1031.39 Q134.561 1028.63 132.617 1026.9 Q130.695 1025.16 127.593 1025.16 Q125.394 1025.16 122.941 1025.92 Q120.51 1026.69 117.732 1028.24 L117.732 1023.52 Q120.556 1022.38 123.01 1021.8 Q125.464 1021.23 127.501 1021.23 Q132.871 1021.23 136.066 1023.91 Q139.26 1026.6 139.26 1031.09 Q139.26 1033.22 138.45 1035.14 Q137.663 1037.04 135.556 1039.63 Q134.978 1040.3 131.876 1043.52 Q128.774 1046.71 123.126 1052.48 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M129.862 662.066 L118.056 680.515 L129.862 680.515 L129.862 662.066 M128.635 657.992 L134.515 657.992 L134.515 680.515 L139.445 680.515 L139.445 684.403 L134.515 684.403 L134.515 692.552 L129.862 692.552 L129.862 684.403 L114.26 684.403 L114.26 679.89 L128.635 657.992 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M127.917 309.55 Q124.769 309.55 122.918 311.702 Q121.089 313.855 121.089 317.605 Q121.089 321.332 122.918 323.508 Q124.769 325.661 127.917 325.661 Q131.066 325.661 132.894 323.508 Q134.746 321.332 134.746 317.605 Q134.746 313.855 132.894 311.702 Q131.066 309.55 127.917 309.55 M137.2 294.897 L137.2 299.156 Q135.441 298.323 133.635 297.883 Q131.853 297.443 130.093 297.443 Q125.464 297.443 123.01 300.568 Q120.58 303.693 120.232 310.012 Q121.598 307.999 123.658 306.934 Q125.718 305.846 128.195 305.846 Q133.404 305.846 136.413 309.017 Q139.445 312.165 139.445 317.605 Q139.445 322.929 136.297 326.147 Q133.149 329.364 127.917 329.364 Q121.922 329.364 118.751 324.781 Q115.58 320.174 115.58 311.448 Q115.58 303.253 119.468 298.392 Q123.357 293.508 129.908 293.508 Q131.667 293.508 133.45 293.855 Q135.255 294.202 137.2 294.897 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M58.657 792.296 L77.5631 792.296 L77.5631 798.185 L28.3562 798.185 L28.3562 792.296 L33.7671 792.296 Q30.5842 790.45 29.0564 787.649 Q27.4968 784.817 27.4968 780.902 Q27.4968 774.409 32.6531 770.367 Q37.8093 766.293 46.212 766.293 Q54.6147 766.293 59.771 770.367 Q64.9272 774.409 64.9272 780.902 Q64.9272 784.817 63.3994 787.649 Q61.8398 790.45 58.657 792.296 M46.212 772.372 Q39.7508 772.372 36.0905 775.045 Q32.3984 777.687 32.3984 782.334 Q32.3984 786.981 36.0905 789.655 Q39.7508 792.296 46.212 792.296 Q52.6732 792.296 56.3653 789.655 Q60.0256 786.981 60.0256 782.334 Q60.0256 777.687 56.3653 775.045 Q52.6732 772.372 46.212 772.372 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M16.4842 746.177 L16.4842 740.766 L70.0516 757.317 L70.0516 762.728 L16.4842 746.177 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M16.4842 737.106 L16.4842 699.771 L21.3858 699.771 L58.5933 729.817 L58.5933 699.039 L64.0042 699.039 L64.0042 737.838 L59.1026 737.838 L21.895 707.792 L21.895 737.106 L16.4842 737.106 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M42.4562 677.968 Q39.719 679.56 39.719 682.774 Q39.719 685.989 42.4562 687.581 Q45.1935 689.204 50.668 689.204 Q56.1743 689.204 58.9116 687.581 Q61.6488 685.989 61.6488 682.774 Q61.6488 679.56 58.9116 677.968 Q56.1743 676.345 50.668 676.345 Q45.1935 676.345 42.4562 677.968 M36.8862 682.774 Q36.8862 677.65 40.4192 674.945 Q43.9522 672.239 50.668 672.239 Q57.4157 672.239 60.9486 674.945 Q64.4816 677.65 64.4816 682.774 Q64.4816 687.931 60.9486 690.636 Q57.4157 693.342 50.668 693.342 Q43.9522 693.342 40.4192 690.636 Q36.8862 687.931 36.8862 682.774 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><polyline clip-path="url(#clip602)" style="stroke:#009af9; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="257.608,1384.24 278.149,1365.49 298.689,1346.75 319.23,1328.02 339.771,1309.31 360.311,1290.62 380.852,1271.95 401.393,1253.3 421.933,1234.69 442.474,1216.1 463.015,1197.56 483.555,1179.06 504.096,1160.6 524.637,1142.19 545.177,1123.84 565.718,1105.54 586.259,1087.3 606.799,1069.12 627.34,1051.01 647.881,1032.97 668.421,1015.01 688.962,997.122 709.503,979.317 730.043,961.597 750.584,943.965 771.125,926.426 791.665,908.982 812.206,891.637 832.747,874.394 853.287,857.258 873.828,840.23 894.369,823.315 914.909,806.516 935.45,789.836 955.991,773.279 976.531,756.848 997.072,740.545 1017.61,724.375 1038.15,708.34 1058.69,692.444 1079.23,676.689 1099.78,661.079 1120.32,645.616 1140.86,630.304 1161.4,615.145 1181.94,600.143 1202.48,585.301 1223.02,570.621 1243.56,556.105 1264.1,541.758 1284.64,527.581 1305.18,513.577 1325.72,499.749 1346.26,486.099 1366.8,472.631 1387.34,459.346 1407.89,446.247 1428.43,433.336 1448.97,420.617 1469.51,408.09 1490.05,395.759 1510.59,383.625 1531.13,371.692 1551.67,359.96 1572.21,348.432 1592.75,337.111 1613.29,325.998 1633.83,315.095 1654.37,304.404 1674.91,293.927 1695.45,283.666 1716,273.623 1736.54,263.798 1757.08,254.195 1777.62,244.813 1798.16,235.656 1818.7,226.725 1839.24,218.02 1859.78,209.544 1880.32,201.297 1900.86,193.282 1921.4,185.499 1941.94,177.949 1962.48,170.633 1983.02,163.553 2003.56,156.71 2024.11,150.105 2044.65,143.737 2065.19,137.61 2085.73,131.722 2106.27,126.075 2126.81,120.67 2147.35,115.506 2167.89,110.586 2188.43,105.909 2208.97,101.475 2229.51,97.286 2250.05,93.341 2270.59,89.6409 2291.13,86.1857 "/>
<path clip-path="url(#clip602)" d="M426.474 1201.07 L426.474 1233.07 L458.474 1233.07 L458.474 1201.07 L426.474 1201.07 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M512.961 1123.7 L512.961 1155.7 L544.961 1155.7 L544.961 1123.7 L512.961 1123.7 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M599.448 1047.26 L599.448 1079.26 L631.448 1079.26 L631.448 1047.26 L599.448 1047.26 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M685.935 972.036 L685.935 1004.04 L717.935 1004.04 L717.935 972.036 L685.935 972.036 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M772.422 898.286 L772.422 930.286 L804.422 930.286 L804.422 898.286 L772.422 898.286 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M858.909 826.272 L858.909 858.272 L890.909 858.272 L890.909 826.272 L858.909 826.272 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M945.396 756.247 L945.396 788.247 L977.396 788.247 L977.396 756.247 L945.396 756.247 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1031.88 688.45 L1031.88 720.45 L1063.88 720.45 L1063.88 688.45 L1031.88 688.45 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1118.37 623.119 L1118.37 655.119 L1150.37 655.119 L1150.37 623.119 L1118.37 623.119 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1204.86 560.474 L1204.86 592.474 L1236.86 592.474 L1236.86 560.474 L1204.86 560.474 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1291.34 500.723 L1291.34 532.723 L1323.34 532.723 L1323.34 500.723 L1291.34 500.723 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1377.83 444.059 L1377.83 476.059 L1409.83 476.059 L1409.83 444.059 L1377.83 444.059 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1464.32 390.655 L1464.32 422.655 L1496.32 422.655 L1496.32 390.655 L1464.32 390.655 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1550.81 340.654 L1550.81 372.654 L1582.81 372.654 L1582.81 340.654 L1550.81 340.654 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1637.29 294.152 L1637.29 326.152 L1669.29 326.152 L1669.29 294.152 L1637.29 294.152 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1723.78 251.141 L1723.78 283.141 L1755.78 283.141 L1755.78 251.141 L1723.78 251.141 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1810.27 211.384 L1810.27 243.384 L1842.27 243.384 L1842.27 211.384 L1810.27 211.384 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1896.75 174.215 L1896.75 206.215 L1928.75 206.215 L1928.75 174.215 L1896.75 174.215 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M1983.24 139.134 L1983.24 171.134 L2015.24 171.134 L2015.24 139.134 L1983.24 139.134 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip602)" d="M2069.73 114.809 L2069.73 146.809 L2101.73 146.809 L2101.73 114.809 L2069.73 114.809 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="442.474" cy="1217.36" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="528.961" cy="1140.11" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="615.448" cy="1063.79" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="701.935" cy="988.677" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="788.422" cy="915.04" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="874.909" cy="843.137" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="961.396" cy="773.22" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1047.88" cy="705.527" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1134.37" cy="640.297" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1220.86" cy="577.748" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1307.34" cy="518.088" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1393.83" cy="461.51" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1480.32" cy="408.187" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1566.81" cy="358.263" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1653.29" cy="311.831" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1739.78" cy="268.886" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1826.27" cy="229.191" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1912.75" cy="192.083" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="1999.24" cy="157.065" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<circle clip-path="url(#clip602)" cx="2085.73" cy="132.791" r="14.4" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="3.2"/>
<path clip-path="url(#clip600)" d="M1834.46 1377.32 L2280.18 1377.32 L2280.18 1169.96 L1834.46 1169.96  Z" fill="#ffffff" fill-rule="evenodd" fill-opacity="1"/>
<polyline clip-path="url(#clip600)" style="stroke:#000000; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="1834.46,1377.32 2280.18,1377.32 2280.18,1169.96 1834.46,1169.96 1834.46,1377.32 "/>
<polyline clip-path="url(#clip600)" style="stroke:#009af9; stroke-linecap:round; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="1858.66,1221.8 2003.81,1221.8 "/>
<path clip-path="url(#clip600)" d="M2043.83 1209.12 L2037.49 1226.32 L2050.2 1226.32 L2043.83 1209.12 M2041.2 1204.52 L2046.5 1204.52 L2059.67 1239.08 L2054.81 1239.08 L2051.66 1230.21 L2036.08 1230.21 L2032.93 1239.08 L2028 1239.08 L2041.2 1204.52 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2086.08 1223.43 L2086.08 1239.08 L2081.82 1239.08 L2081.82 1223.57 Q2081.82 1219.89 2080.39 1218.06 Q2078.95 1216.23 2076.08 1216.23 Q2072.63 1216.23 2070.64 1218.43 Q2068.65 1220.63 2068.65 1224.42 L2068.65 1239.08 L2064.37 1239.08 L2064.37 1213.15 L2068.65 1213.15 L2068.65 1217.18 Q2070.18 1214.84 2072.24 1213.68 Q2074.32 1212.52 2077.03 1212.52 Q2081.5 1212.52 2083.79 1215.3 Q2086.08 1218.06 2086.08 1223.43 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2106.36 1226.04 Q2101.2 1226.04 2099.2 1227.22 Q2097.21 1228.4 2097.21 1231.25 Q2097.21 1233.52 2098.7 1234.86 Q2100.2 1236.18 2102.77 1236.18 Q2106.31 1236.18 2108.44 1233.68 Q2110.59 1231.16 2110.59 1226.99 L2110.59 1226.04 L2106.36 1226.04 M2114.85 1224.28 L2114.85 1239.08 L2110.59 1239.08 L2110.59 1235.14 Q2109.14 1237.5 2106.96 1238.64 Q2104.78 1239.75 2101.64 1239.75 Q2097.65 1239.75 2095.29 1237.52 Q2092.96 1235.28 2092.96 1231.53 Q2092.96 1227.15 2095.87 1224.93 Q2098.81 1222.71 2104.62 1222.71 L2110.59 1222.71 L2110.59 1222.29 Q2110.59 1219.35 2108.65 1217.76 Q2106.73 1216.14 2103.23 1216.14 Q2101.01 1216.14 2098.9 1216.67 Q2096.8 1217.2 2094.85 1218.27 L2094.85 1214.33 Q2097.19 1213.43 2099.39 1212.99 Q2101.59 1212.52 2103.67 1212.52 Q2109.3 1212.52 2112.08 1215.44 Q2114.85 1218.36 2114.85 1224.28 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2123.63 1203.06 L2127.89 1203.06 L2127.89 1239.08 L2123.63 1239.08 L2123.63 1203.06 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2147.58 1241.48 Q2145.78 1246.11 2144.07 1247.52 Q2142.35 1248.94 2139.48 1248.94 L2136.08 1248.94 L2136.08 1245.37 L2138.58 1245.37 Q2140.34 1245.37 2141.31 1244.54 Q2142.28 1243.7 2143.46 1240.6 L2144.23 1238.66 L2133.74 1213.15 L2138.26 1213.15 L2146.36 1233.43 L2154.46 1213.15 L2158.97 1213.15 L2147.58 1241.48 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2169.07 1205.79 L2169.07 1213.15 L2177.84 1213.15 L2177.84 1216.46 L2169.07 1216.46 L2169.07 1230.53 Q2169.07 1233.7 2169.92 1234.61 Q2170.8 1235.51 2173.46 1235.51 L2177.84 1235.51 L2177.84 1239.08 L2173.46 1239.08 Q2168.53 1239.08 2166.66 1237.25 Q2164.78 1235.39 2164.78 1230.53 L2164.78 1216.46 L2161.66 1216.46 L2161.66 1213.15 L2164.78 1213.15 L2164.78 1205.79 L2169.07 1205.79 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2183.44 1213.15 L2187.7 1213.15 L2187.7 1239.08 L2183.44 1239.08 L2183.44 1213.15 M2183.44 1203.06 L2187.7 1203.06 L2187.7 1208.45 L2183.44 1208.45 L2183.44 1203.06 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2215.27 1214.14 L2215.27 1218.13 Q2213.46 1217.13 2211.63 1216.64 Q2209.83 1216.14 2207.98 1216.14 Q2203.83 1216.14 2201.54 1218.77 Q2199.25 1221.39 2199.25 1226.14 Q2199.25 1230.88 2201.54 1233.52 Q2203.83 1236.14 2207.98 1236.14 Q2209.83 1236.14 2211.63 1235.65 Q2213.46 1235.14 2215.27 1234.14 L2215.27 1238.08 Q2213.49 1238.91 2211.57 1239.33 Q2209.67 1239.75 2207.51 1239.75 Q2201.66 1239.75 2198.21 1236.07 Q2194.76 1232.39 2194.76 1226.14 Q2194.76 1219.79 2198.23 1216.16 Q2201.73 1212.52 2207.79 1212.52 Q2209.76 1212.52 2211.63 1212.94 Q2213.51 1213.33 2215.27 1214.14 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2234.46 1226.04 Q2229.3 1226.04 2227.31 1227.22 Q2225.32 1228.4 2225.32 1231.25 Q2225.32 1233.52 2226.8 1234.86 Q2228.3 1236.18 2230.87 1236.18 Q2234.41 1236.18 2236.54 1233.68 Q2238.69 1231.16 2238.69 1226.99 L2238.69 1226.04 L2234.46 1226.04 M2242.95 1224.28 L2242.95 1239.08 L2238.69 1239.08 L2238.69 1235.14 Q2237.24 1237.5 2235.06 1238.64 Q2232.88 1239.75 2229.74 1239.75 Q2225.76 1239.75 2223.39 1237.52 Q2221.06 1235.28 2221.06 1231.53 Q2221.06 1227.15 2223.97 1224.93 Q2226.91 1222.71 2232.72 1222.71 L2238.69 1222.71 L2238.69 1222.29 Q2238.69 1219.35 2236.75 1217.76 Q2234.83 1216.14 2231.33 1216.14 Q2229.11 1216.14 2227.01 1216.67 Q2224.9 1217.2 2222.95 1218.27 L2222.95 1214.33 Q2225.29 1213.43 2227.49 1212.99 Q2229.69 1212.52 2231.77 1212.52 Q2237.4 1212.52 2240.18 1215.44 Q2242.95 1218.36 2242.95 1224.28 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2251.73 1203.06 L2255.99 1203.06 L2255.99 1239.08 L2251.73 1239.08 L2251.73 1203.06 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M1908.48 1250.88 L1908.48 1296.39 L1953.99 1296.39 L1953.99 1250.88 L1908.48 1250.88 Z" fill="#e26f46" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="4.55111"/>
<path clip-path="url(#clip600)" d="M2028 1256.36 L2047.86 1256.36 L2047.86 1260.29 L2032.68 1260.29 L2032.68 1270.48 L2046.38 1270.48 L2046.38 1274.41 L2032.68 1274.41 L2032.68 1290.92 L2028 1290.92 L2028 1256.36 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2052.05 1280.68 L2052.05 1264.99 L2056.31 1264.99 L2056.31 1280.52 Q2056.31 1284.2 2057.75 1286.05 Q2059.18 1287.88 2062.05 1287.88 Q2065.5 1287.88 2067.49 1285.68 Q2069.51 1283.48 2069.51 1279.69 L2069.51 1264.99 L2073.77 1264.99 L2073.77 1290.92 L2069.51 1290.92 L2069.51 1286.93 Q2067.96 1289.29 2065.89 1290.45 Q2063.86 1291.59 2061.15 1291.59 Q2056.68 1291.59 2054.37 1288.81 Q2052.05 1286.03 2052.05 1280.68 M2062.77 1264.36 L2062.77 1264.36 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2082.54 1254.9 L2086.8 1254.9 L2086.8 1290.92 L2082.54 1290.92 L2082.54 1254.9 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2095.71 1254.9 L2099.97 1254.9 L2099.97 1290.92 L2095.71 1290.92 L2095.71 1254.9 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><circle clip-path="url(#clip600)" cx="1931.23" cy="1325.48" r="20.48" fill="#3da44d" fill-rule="evenodd" fill-opacity="1" stroke="#000000" stroke-opacity="1" stroke-width="4.55111"/>
<path clip-path="url(#clip600)" d="M2028 1308.2 L2047.86 1308.2 L2047.86 1312.13 L2032.68 1312.13 L2032.68 1322.32 L2046.38 1322.32 L2046.38 1326.25 L2032.68 1326.25 L2032.68 1342.76 L2028 1342.76 L2028 1308.2 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2055.27 1308.2 L2062.24 1308.2 L2071.06 1331.71 L2079.92 1308.2 L2086.89 1308.2 L2086.89 1342.76 L2082.33 1342.76 L2082.33 1312.41 L2073.42 1336.11 L2068.72 1336.11 L2059.81 1312.41 L2059.81 1342.76 L2055.27 1342.76 L2055.27 1308.2 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /><path clip-path="url(#clip600)" d="M2096.17 1308.2 L2103.14 1308.2 L2111.96 1331.71 L2120.83 1308.2 L2127.79 1308.2 L2127.79 1342.76 L2123.23 1342.76 L2123.23 1312.41 L2114.32 1336.11 L2109.62 1336.11 L2100.71 1312.41 L2100.71 1342.76 L2096.17 1342.76 L2096.17 1308.2 Z" fill="#000000" fill-rule="nonzero" fill-opacity="1" /></svg>

```

The surface pressures can also be plotted. Note that it is only possible to plot linear meshes, meaing that we must remove the quadratic parts.

````julia
data_mesh,data_viz = create_vizualization_data(mesh,p_fmm)
fig, ax, hm = viz(data_mesh;showfacets=true, color=abs.(data_viz/Z₀))
wgl.Colorbar(fig[1,2],label="|p/Z₀|");
````

![](3d_oscilating_sphere.png)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

