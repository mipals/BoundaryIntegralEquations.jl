```@meta
EditURL = "<unknown>/../examples/3d_element_usage.jl"
```

# 3D Element Usage
For 3D boundary integral equations the elements of interest are `SurfaceFunctions`.
These are defined by two local coordinates, $u$ and $v$
```math
  \mathbf{N}(u,v)
```

````julia
import WGLMakie as Mke
````

Set the default resolution to something that fits the Documenter theme

````julia
Mke.set_theme!(resolution=(800, 800))

using IntegralEquations
import IntegralEquations: set_nodal_interpolation!
linear_triangular = TriangularLinear(3)
````

````
SurfaceFunction Defined by:    	 TriangularLinear{Float64}
Number of gauss nodes:         	 3

````

Plotting Gaussian Integration Points

````julia
Mke.scatter(linear_triangular.gauss_u, linear_triangular.gauss_v)
````

```@raw html
<div data-jscall-id="2"><script data-jscall-id="3" type="text/javascript">
    function register_resize_handler(remote_origin) {
        function resize_callback(event) {
            if (event.origin !== remote_origin) {
                return;
            }
            const uuid = event.data[0];
            const width = event.data[1];
            const height = event.data[2];
            const iframe = document.getElementById('9fc62682-3d9b-41f7-b76e-ff1df2c7e2cc');
            if (iframe) {
                iframe.style.width = width + "px";
                iframe.style.height = height + "px";
            }
        }
        if (window.addEventListener) {
            window.addEventListener("message", resize_callback, false);
        } else if (window.attachEvent) {
            window.attachEvent("onmessage", resize_callback);
        }
    }
    register_resize_handler('http://127.0.0.1:9284')

</script><iframe scrolling="no" id="9fc62682-3d9b-41f7-b76e-ff1df2c7e2cc" data-jscall-id="1" src="http://127.0.0.1:9284/9fc62682-3d9b-41f7-b76e-ff1df2c7e2cc" style="position: relative; display: block; width: 100%; height: 100%; padding: 0; overflow: hidden; border: none"></iframe></div>

```

We can instead

````julia
reference_triangle = TriangularLinear(3)
set_nodal_interpolation!(reference_triangle)
Mke.scatter!(reference_triangle.gauss_u, reference_triangle.gauss_v)
Mke.lines!([reference_triangle.gauss_u; 0.0],[reference_triangle.gauss_v; 0.0])
Mke.current_figure()
````

```@raw html
<div data-jscall-id="5"><script data-jscall-id="6" type="text/javascript">
    function register_resize_handler(remote_origin) {
        function resize_callback(event) {
            if (event.origin !== remote_origin) {
                return;
            }
            const uuid = event.data[0];
            const width = event.data[1];
            const height = event.data[2];
            const iframe = document.getElementById('c3f98684-5018-41ff-896b-735b4aced47e');
            if (iframe) {
                iframe.style.width = width + "px";
                iframe.style.height = height + "px";
            }
        }
        if (window.addEventListener) {
            window.addEventListener("message", resize_callback, false);
        } else if (window.attachEvent) {
            window.attachEvent("onmessage", resize_callback);
        }
    }
    register_resize_handler('http://127.0.0.1:9284')

</script><iframe scrolling="no" id="c3f98684-5018-41ff-896b-735b4aced47e" data-jscall-id="4" src="http://127.0.0.1:9284/c3f98684-5018-41ff-896b-735b4aced47e" style="position: relative; display: block; width: 100%; height: 100%; padding: 0; overflow: hidden; border: none"></iframe></div>

```

Weights should add to area of triangle

````julia
sum(linear_triangular.weights) ≈ 0.5
````

````
true
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

