using Test, Aqua, BoundaryIntegralEquations

Aqua.test_undefined_exports(BoundaryIntegralEquations)
Aqua.test_project_extras(BoundaryIntegralEquations)
Aqua.test_unbound_args(BoundaryIntegralEquations)
Aqua.test_ambiguities(BoundaryIntegralEquations)
Aqua.test_stale_deps(BoundaryIntegralEquations)
Aqua.test_piracies(BoundaryIntegralEquations)
# Aqua.test_deps_compat(BoundaryIntegralEquations) # We need FMM3D to be registered
# Aqua.test_persistent_tasks(BoundaryIntegralEquations) # We need FMM3D to be registered
