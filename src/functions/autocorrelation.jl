# =============================================================================
# Autocorrelation function and spectrum of a wave packet.
#
# A wave packet ψ₀ evolves as ψ(t) = e^{-iHt} ψ₀ (ħ = 1). Expanding it in the
# eigenstates φₙ of H,  ψ₀ = Σₙ cₙ φₙ  with  cₙ = ⟨φₙ|ψ₀⟩, gives
#
#   autocorrelation   C(t) = ⟨ψ₀|ψ(t)⟩ = Σₙ |cₙ|² e^{-iεₙt}
#   spectrum          σ(ω) = (1/2π) ∫ C(t) e^{iωt} dt = Σₙ |cₙ|² δ(ω − εₙ)
#
# Since C(−t) = C(t)*, σ(ω) = (1/π) Re ∫₀^∞ C(t) e^{iωt} dt. A finite propagation
# time requires a damping (window) w(t), which broadens every line:
#   w(t) = e^{−γt}         →  Lorentzian lines  (1/π) γ / ((ω−εₙ)² + γ²)
#   w(t) = e^{−t²/(2τ²)}   →  Gaussian lines    τ/√(2π) e^{−(ω−εₙ)²τ²/2}
# `autocorrelation_spectrum` evaluates the broadened spectrum directly from the
# eigenpairs; `fourier_transform` computes it from sampled values of C(t) (e.g.
# from a time propagation), with the same conventions.
# =============================================================================

"""
    expansion_coefficients(ψ₀, states, dΩ)

Coefficients ``c_n = \\langle\\phi_n|\\psi_0\\rangle = \\int \\phi_n^*\\,\\psi_0\\,d\\Omega`` of
the wave packet `ψ₀` in the eigenstates `states` (e.g. `sol.states` of
[`solve_tise`](@ref)). `ψ₀` can be a Julia function of the mesh point or a finite
element function.

``\\sum_n |c_n|^2`` equals the norm of `ψ₀` only if the eigenstates span it: check it
to make sure that enough eigenstates were computed.

!!! note
    Eigenvectors of a degenerate level are linearly independent but not necessarily
    orthogonal; in that case orthonormalize them before computing the coefficients.
"""
function expansion_coefficients(ψ₀, states::AbstractVector, dΩ)
    return [sum(∫(conj(φ) * ψ₀) * dΩ) for φ in states]
end

"""
    autocorrelation(coefficients, energies, t)

Autocorrelation ``C(t) = \\sum_n |c_n|^2 e^{-i\\varepsilon_n t}`` (ħ = 1) at the times
`t` (a number or a collection), from the expansion coefficients of the wave packet
(see [`expansion_coefficients`](@ref)) and the eigenvalues.

``|C(t)|^2`` is the survival probability of the initial state.
"""
function autocorrelation(coefficients::AbstractVector, energies::AbstractVector, t::Real)
    length(coefficients) == length(energies) || throw(ArgumentError(
        "got $(length(coefficients)) coefficients and $(length(energies)) energies"))
    return sum(abs2(c) * cis(-ε * t) for (c, ε) in zip(coefficients, energies))
end
autocorrelation(coefficients, energies, times) = [autocorrelation(coefficients, energies, t) for t in times]

"""
    autocorrelation_spectrum(coefficients, energies, ω; width, lineshape = :lorentzian)

Spectrum of the wave packet, ``\\sigma(\\omega) = \\sum_n |c_n|^2 L(\\omega-\\varepsilon_n)``,
on the frequencies (energies) `ω`, with every line broadened by
- `lineshape = :lorentzian`: ``L(\\Delta) = \\frac{1}{\\pi}\\frac{w}{\\Delta^2 + w^2}``
  (half width at half maximum `w`);
- `lineshape = :gaussian`: ``L(\\Delta) = \\frac{1}{\\sqrt{2\\pi}\\,w}e^{-\\Delta^2/(2w^2)}``
  (standard deviation `w`).

These are the Fourier transforms of ``C(t)`` damped by ``e^{-wt}`` and by
``e^{-w^2t^2/2}``, respectively, so they coincide with
[`fourier_transform`](@ref) of the sampled autocorrelation with the matching window.
Each line has unit area, so ``\\int\\sigma\\,d\\omega = \\sum_n|c_n|^2``.
"""
function autocorrelation_spectrum(coefficients::AbstractVector, energies::AbstractVector, ω;
                                  width::Real, lineshape::Symbol = :lorentzian)
    width > 0 || throw(ArgumentError("the line width must be positive"))
    line = if lineshape == :lorentzian
        Δ -> width / (π * (Δ^2 + width^2))
    elseif lineshape == :gaussian
        Δ -> exp(-Δ^2 / (2width^2)) / (sqrt(2π) * width)
    else
        throw(ArgumentError("lineshape must be :lorentzian or :gaussian, got :$lineshape"))
    end
    weights = abs2.(coefficients)
    return [sum(w * line(ωᵢ - ε) for (w, ε) in zip(weights, energies)) for ωᵢ in ω]
end

"""
    fourier_transform(t, C, ω; window = :none, width = 0.0)

Spectrum ``\\sigma(\\omega) = \\frac{1}{\\pi}\\,\\mathrm{Re}\\int_0^{T} C(t)\\,w(t)\\,e^{i\\omega t}\\,dt``
of an autocorrelation function sampled at equally spaced times `t` (starting at
``t=0``), evaluated at the frequencies `ω` with the trapezoidal rule.

# Windows
- `:none`: ``w = 1``. Lines have width ``\\sim\\pi/T`` and side lobes.
- `:exponential`: ``w = e^{-\\gamma t}`` with ``\\gamma =`` `width`. Lorentzian lines
  (HWHM ``\\gamma``); requires ``\\gamma T \\gg 1``.
- `:gaussian`: ``w = e^{-\\gamma^2t^2/2}`` with ``\\gamma =`` `width`. Gaussian lines
  (standard deviation ``\\gamma``); requires ``\\gamma T \\gg 1``.

With the same `width`, the result agrees with
[`autocorrelation_spectrum`](@ref)`(c, ε, ω; width, lineshape)`. The cost is
``O(N_t N_\\omega)``.
"""
function fourier_transform(t::AbstractVector, C::AbstractVector, ω; window::Symbol = :none,
                           width::Real = 0.0)
    length(t) == length(C) || throw(ArgumentError("t and C must have the same length"))
    length(t) >= 2 || throw(ArgumentError("at least two time samples are needed"))
    Δt = t[2] - t[1]
    all(abs.(diff(t) .- Δt) .<= 1e-9 * abs(Δt)) || throw(ArgumentError("the times must be equally spaced"))
    abs(t[1]) <= 1e-12 * abs(Δt) || throw(ArgumentError("the first time must be t = 0"))
    damping = if window == :none
        ones(length(t))
    elseif window == :exponential
        width > 0 || throw(ArgumentError("the exponential window needs width > 0"))
        exp.(-width .* t)
    elseif window == :gaussian
        width > 0 || throw(ArgumentError("the Gaussian window needs width > 0"))
        exp.(-(width .* t) .^ 2 ./ 2)
    else
        throw(ArgumentError("window must be :none, :exponential or :gaussian, got :$window"))
    end
    trapezoid = fill(Δt, length(t)); trapezoid[1] = trapezoid[end] = Δt / 2
    weighted = trapezoid .* damping .* C
    return [real(sum(weighted .* cis.(ωᵢ .* t))) / π for ωᵢ in ω]
end
