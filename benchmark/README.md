This directory contains benchmarks for CSetAutomorphisms. To run all the
benchmarks, launch `julia --project=benchmark` and enter:

``` julia
using PkgBenchmark
import CSetAutomorphisms

benchmarkpkg(CSetAutomorphisms)
```

To run a specific set of benchmarks, use the `script` keyword argument, for
example:

``` julia
benchmarkpkg(CSetAutomorphisms; script="benchmark/benchmarks.jl")
```
