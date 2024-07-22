# VicsekModel

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://shiraishi-mu.github.io/VicsekModel.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://shiraishi-mu.github.io/VicsekModel.jl/dev/)
[![Build Status](https://github.com/shiraishi-mu/VicsekModel.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/shiraishi-mu/VicsekModel.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package provides a simple implementation of the Vicsek model in Julia.

Algorithm is based on [the python code by Dr. Frances Turci](https://francescoturci.net/2020/06/19/minimal-vicsek-model-in-python/).

## Installation

```julia
julia> ]
pkg> add VicsekModel
```

## Usage

### Simulation with default parameters and saveing the animation

```julia
using VicsekModel

Tmax = 1000

vicsekmodel_with_anim(Tmax)
```

### Simulation with custom parameters and saveing the animation

```julia
using VicsekModel

Tmax = 1000
$\rho$ = 0.5
$\eta$ = 0.5
$\nu$ = 0.03

p = VicsekModelParameters{Float64, Int64}(η=0.3, ρ=0.4)
println(p)
var = VicsekModelVariables(p)

anim = @animate for t in 1:Tmax
    update_with_anim(t, var, p)
end
gif(anim, "animation.gif", fps = 15)
```

### Simulation without saving the animation

```julia
using VicsekModel

Tmax = 1000
$\rho$ = 0.5
$\eta$ = 0.5
$\nu$ = 0.03

p = VicsekModelParameters{Float64, Int64}(η=0.3, ρ=0.4)
var = VicsekModelVariables(p)

for t in 1:Tmax
    update(t, var, p)
end
```

### Simulation with custom initial state

```julia
using VicsekModel

Tmax = 1000
$\rho$ = 0.5
$\eta$ = 0.5
$\nu$ = 0.03

p = VicsekModelParameters{Float64, Int64}(η=0.3, ρ=0.4)
var = VicsekModelVariables(p)

# Custom initial state
var.θ = rand(Float64, p.N) * 2π
var.v = [cos.(var.θ) sin.(var.θ)]

for t in 1:Tmax
    update(t, var, p)
end
```

### Simulation with another algorithm

```julia
using VicsekModel

Tmax = 1000
$\rho$ = 0.5
$\eta$ = 0.5
$\nu$ = 0.03

p = VicsekModelParameters{Float64, Int64}(η=0.3, ρ=0.4, alg=algForLoop())
var = VicsekModelVariables(p)

for t in 1:Tmax
    update(t, var, p)
end
```

## Documentation

- [Stable](https://shiraishi-mu.github.io/VicsekModel.jl/stable/)
- [Dev](https://shiraishi-mu.github.io/VicsekModel.jl/dev/)
- [Tutorial](https://shiraishi-mu.github.io/VicsekModel.jl/dev/tutorial/)
