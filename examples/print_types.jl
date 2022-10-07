using AbstractTrees
using IntegralEquations

AbstractTrees.children(T::Type) = subtypes(T)

print_tree(ShapeFunction)
print_tree(CurveFunction)
print_tree(SurfaceFunction)
