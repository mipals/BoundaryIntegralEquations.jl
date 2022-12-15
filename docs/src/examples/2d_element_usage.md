```@meta
EditURL = "<unknown>/../examples/2d_element_usage.jl"
```

# 2D Element Usage

Our function $f$ on the annulus:

````julia
import WGLMakie as Mke
````

Set the default resolution to something that fits the Documenter theme

````julia
Mke.set_theme!(resolution=(800, 400))
````

notes
nodes,weights  = gausslegendre(20)
m1 = [1;0.0]
m2 = [1.5;0.5]
m3 = [1;1.0]
coordinates = [m1 m2 m3]
elementInterpolation = quadratic(nodes')
interpolation = similar(coordinates*elementInterpolation)
mul!(interpolation,coordinates,elementInterpolation)
scatter(coordinates[1,:],coordinates[2,:])
plot!(interpolation[1,:],interpolation[2,:])
scatter!(interpolation[1,:],interpolation[2,:])

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

