module VicsekModel
    using NearestNeighbors
    using Random
    using Distributions
    using Parameters
    using Plots

    abstract type VicsekModelAlgorithm end
    struct algKDTree <: VicsekModelAlgorithm end
    struct algForLoop <: VicsekModelAlgorithm end

    is_Alg(::VicsekModelAlgorithm) = ForLoop()
    is_Alg(::algKDTree) = algKDTree()
    is_Alg(::algForLoop) = algForLoop()

    @kwdef struct VicsekModelParameters{Tf <: Real, Ti <: Integer}
        L::Tf = 32.0
        ρ::Tf = 3.0
        N::Ti = Int(ρ * L^2)

        r₀::Tf = 1.0
        Δt::Tf = 1.0
        factor::Tf = 0.5
        v₀::Tf = factor * r₀ / Δt
        η::Tf = 0.15

        alg::VicsekModelAlgorithm = algKDTree()
    end 

    function Base.show(io::IO, vmpara::VicsekModelParameters)
        println("L: ", vmpara.L)
        println("ρ: ", vmpara.ρ)
        println("N: ", vmpara.N)
        println("r₀: ", vmpara.r₀)
        println("Δt: ", vmpara.Δt)
        println("factor: ", vmpara.factor)
        println("v₀: ", vmpara.v₀)
        println("η: ", vmpara.η)
        println("alg: ", typeof(vmpara.alg))
    end

    struct VicsekModelVariables{Tf <: Real}
        pos::Array{Tf, 2}
        vel::Array{Tf, 2}
        θ::Vector{Tf}
        S::Vector{Complex}
    end

    function VicsekModelVariables(p::VicsekModelParameters, seed::Integer=1)
        @unpack v₀, N, L = p
        Random.seed!(seed)
        pos = rand(Uniform(0, L), N, 2)
        θ = rand(Uniform(-π, π), N)
        vel = hcat(v₀ .* cos.(θ), v₀ .*sin.(θ))
        S = zeros(Complex, N)
        return VicsekModelVariables(pos, vel, θ, S)
    end

    function Base.show(io::IO, var::VicsekModelVariables)
        println("位置 (pos): ", var.pos)
        println("速度 (vel): ", var.vel)
        println("角度 (θ): ", var.θ)
        println("強度 (S): ", var.S)
    end

    update_with_anim(i, var, p) = update_with_anim(is_Alg(p.alg), i, var, p)

    function calculate_weights!(tree::NNTree, pos::Array{Tf, 2}, θ::Vector{Tf}, S::Vector{Complex}, r₀::Tf, ix) where Tf
        @inbounds for j in ix
            neighbors = inrange(tree, pos[j,:], r₀, true)
            @inbounds for k in neighbors
                S[j] += cis(θ[k] * im)
            end
        end
    end

    function update_position!(pos::Array{Tf, 2}, vel::Array{Tf, 2}, θ::Vector{Tf}, S::Vector{Complex}, v₀::Tf, η::Tf, ix) where Tf
        @inbounds for j in ix
            θ[j] = S[j] + η * Random.rand(Uniform(-π, π))
            vel[j, 1] = cos(θ[j])
            vel[j, 2] = sin(θ[j])
            pos[j, 1] += v₀ * vel[j, 1]
            pos[j, 2] += v₀ * vel[j, 2]
        end
    end

    function periodic_boundary!(pos::Array{Tf, 2}, L::Tf, ix::Base.OneTo) where Tf
        @inbounds for j in ix
            pos[j, 1] = ifelse(pos[j, 1] > L, pos[j, 1]-L, ifelse(pos[j, 1]<0, pos[j, 1]+L, pos[j, 1]))
            pos[j, 2] = ifelse(pos[j, 2] > L, pos[j, 2]-L, ifelse(pos[j, 2]<0, pos[j, 2]+L, pos[j, 2]))
        end

    end

    function update_with_anim(::algKDTree, i::Ti, var::VicsekModelVariables{Tf}, p::VicsekModelParameters{Tf, Ti}) where {Tf, Ti}
        @unpack pos, vel, θ, S = var
        @unpack L, ρ, N, r₀, Δt, factor, v₀, η = p
        ix = axes(pos, 1)

        tree = KDTree(transpose(pos))

        calculate_weights!(tree, pos, θ, S, r₀, ix)

        update_position!(pos, vel, θ, S, v₀, η, ix)

        periodic_boundary!(pos, L, ix)

        fig = plot(size=(600, 600), xlim=(0, L), ylim=(0, L))
        quiver!(fig, pos[:,1], pos[:,2], quiver=(vel[:,1], vel[:,2]), framestyle=:box, line_z=repeat(θ, inner=4), c=:viridis)
        return fig
    end

    function vicsekmodel_with_anim(Tmax)
        p = VicsekModelParameters{Float64, Int64}()
        println(p)
        var = VicsekModelVariables(p)

        anim = @animate for t in 1:Tmax
            update_with_anim(t, var, p)
        end
        gif(anim, "animation.gif", fps = 15)
    end
end

