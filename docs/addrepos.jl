# Add PackageNursery registry for CI
using Pkg
Pkg.Registry.add("General")
Pkg.Registry.add(RegistrySpec(url = "https://github.com/mipals/PackageRepository"))
