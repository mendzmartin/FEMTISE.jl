#=
References:
    https://gist.github.com/Balaje/102485bb14ec6daf677f938fbd8f3ebb
    https://docs.scipy.org/doc/scipy/tutorial/arpack.html
    https://github.com/JuliaLinearAlgebra/Arpack.jl/blob/master/docs/src/eigs.md
    https://github.com/JuliaPDE/SurveyofPDEPackages#mm
=#

struct EigenOperator{H<:AbstractMatrix,E<:AbstractMatrix} #<: NonlinearOperator
    hamiltonian::H
    energy::E
end

struct EigenProblem <: FEOperator
    trial::FESpace
    test::FESpace
    op::EigenOperator
    nev::Int64
    which::Symbol #LM or SM
    explicittransform::Symbol #:shiftinvert
    tol::Float64
    maxiter::Int64
    sigma
end

export eigen_problem
"""
    eigen_problem(weakformₖ,weakformₘ,test,trial;; <keyword arguments>)

Define eigen problem as an input to solve function where we compute eigen problem by Arpack eigs function.

    ...
    # Arguments
    - `weakformₖ::Function`: forma bilineal lado izquierdo de la formulación débil
    - `weakformₘ::Function`: forma bilineal lado derecho de la formulación débil
    - `test::FESpace`: espacio de prueba, puede ser MultiFieldFESpace
    - `trial::FESpace`: espacio de solución, puede ser MultiFieldFESpace
    - `nev::Int64=10`: número de autovalores requeridos
    - `tol::Float64=10e-6`: relative tolerance for convergence of Ritz values
    - `maxiter::Integer=100`: maximum number of iterations
    - `explicittransform::Symbol=:none`: shift and invert should be explicitly invoked in julia code
                               `=:auto`:
                               `=:shiftinvert`:
    - `sigma::Float64=1.0`: the level shift used in inverse iteration.
    - `which::Symbol=:LM`: eigenvalues of largest magnitude (default)
                   `=:SM`: eigenvalues of smallest magnitude
                   `=:LR`: eigenvalues of largest real part
                   `=:SR`: eigenvalues of smallest real part
                   `=:LI`: eigenvalues of largest imaginary part (nonsymmetric or complex matrix only)
                   `=:SI`: eigenvalues of smallest imaginary part (nonsymmetric or complex matrix only)
                   `=:BE`: compute half of the eigenvalues from each end of the spectrum, biased in favor of the high end. (real symmetric matrix only)
    ...
"""
function eigen_problem(weakformₖ::Function,weakformₘ::Function,test::FESpace,trial::FESpace;
    nev::Int64=10,which::Symbol=:LM,explicittransform::Symbol=:none,tol::Float64=10^(-6),
    maxiter::Int64=100,sigma=0.0)
    # source vector (always need to be zero for eigen problems)
    F(v) = 0.0;
    opH = AffineFEOperator(weakformₖ, F, test, trial)
    opE = AffineFEOperator(weakformₘ, F, test, trial)
    op = EigenOperator(opH.op.matrix,opE.op.matrix)
    return EigenProblem(trial,test,op,nev,which,explicittransform,tol,maxiter,sigma)
end

export solve
"""
    solve(prob)

Compute eigen problem by Arpack eigs function and return eigenvalues and eigenvectors.

    ...
    # Arguments
    ...
    `prob::EigenProblem`: problem deinition
"""
function solve(prob::EigenProblem)
    ϵ,eigenvecs = eigs(prob.op.hamiltonian,prob.op.energy;
        nev=prob.nev,which=prob.which,explicittransform=prob.explicittransform,
        tol=prob.tol,maxiter=prob.maxiter,sigma=prob.sigma)
    ϕ=[FEFunction(prob.trial, eigenvecs[:,m]) for m=1:prob.nev]
    return ϵ,ϕ
end