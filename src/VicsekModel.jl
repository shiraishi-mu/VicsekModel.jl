module VicsekModel
    using Random
    using Distributions
    using Parameters

    @kwdef struct VicsekModelParameters{Tf, Ti}
        L::Tf = 32.0
        ρ::Tf = 3.0
        N::Ti = Int(ρ * L^2)

        r₀::Tf = 1.0
        Δt::Tf = 1.0
        factor::Tf = 0.5
        v₀::Tf = factor * r₀ / Δt
        η:Tf = 0.15
    end 

    struct VicsekModelVariables{Tf}
        pos::Array{Tf, 2}
        vel::Array{Tf, 2}
        θ::Vector{Tf}
    end

    abstract type VicsekModelAlgorithm end
    struct KDTree <: VicsekModelAlgorithm end
    struct ForLoop <: VicsekModelAlgorithm end

    is_Alg(::VicsekModelAlgorithm) = ForLoop()
    is_Alg(::KDTree) = KDTree()
    is_Alg(::ForLoop) = ForLoop()

    update_with_anim(i) = update_with_anim(VicsekModelAlgorithm(), i, var, p)

    function update_with_anim(::KDTree, i::Ti, var::VicsekModelVariables{Tf}, p::VicsekModelParameters{Tf, Ti}) where {Tf, Ti}
        @unpack pos, vel, θ = var
        @unpack L, ρ, N, r₀, Δt, factor, v₀, η = p

        tree = KDTree(transpose(pos))

        @inbounds for j in 1:N
            neighbors = inrange(tree, pos[j,:], r₀, true)
            @inbounds for k in neighbors
                S[j] += exp(θ[k] * im)
            end
        end

        @inbounds for j in 1:N
            θ[j] = S[j] + η * random(Uniform(-π, π))
            vel[j, 1] = cos(θ[j])
            vel[j, 2] = sin(θ[j])
            pos[j, 1] += v₀ * vel[j, 1]
            pos[j, 2] += v₀ * vel[j, 2]
        end

        @inbounds for j in 1:N
            pos[j, 1] = mod(pos[j, 1], L)
            pos[j, 2] = mod(pos[j, 2], L)
        end

        fig = plot(size=(600, 600), xlim=(0, L), ylim=(0, L))
        quiver!(fig, pos[:,1], pos[:,2], quiver=(vel[:,1], vel[:,2]), framestyle=:box, line_z=repeat(orient, inner=4), c=:viridis)
        return fig

        function vicsekmodel_with_anim()
            for t in 1:typemax
                update_with_anim(t)
            end
        end
    end

end

