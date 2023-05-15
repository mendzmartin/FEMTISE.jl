# Taken from:  https://gist.github.com/Balaje/102485bb14ec6daf677f938fbd8f3ebb
# Useful link: https://docs.scipy.org/doc/scipy/tutorial/arpack.html

struct EigenOperator{K<:AbstractMatrix,M<:AbstractMatrix} #<: NonlinearOperator
    stima::K
    massma::M
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

# https://github.com/JuliaLinearAlgebra/Arpack.jl/blob/master/docs/src/eigs.md
# https://github.com/JuliaPDE/SurveyofPDEPackages#mm

export EigenProblem
"""
    function EigenProblem(weakformₖ,weakformₘ,trial,test;nev,which,explicittransform,tol,maxiter,sigma)\n
    weakformₖ::Function,  (Forma bilineal lado izquierdo de la formulación débil) \n
    weakformₘ::Function,  (Forma bilineal lado derecho de la formulación débil) \n
    test::FESpace,  (Espacio de prueba, puede ser MultiFieldFESpace) \n
    trial::FESpace; (Espacio de solución, puede ser MultiFieldFESpace) \n
    nev::Int64=10,  (Número de autovalores requeridos) \n
    which::Symbol=:LM, \n
    #:LM	eigenvalues of largest magnitude (default) \n
    #:SM	eigenvalues of smallest magnitude \n
    #:LR	eigenvalues of largest real part \n
    #:SR	eigenvalues of smallest real part \n
    #:LI	eigenvalues of largest imaginary part (nonsymmetric or complex A only) \n
    #:SI	eigenvalues of smallest imaginary part (nonsymmetric or complex A only) \n
    :BE	compute half of the eigenvalues from each end of the spectrum, biased in favor of the high end. (real symmetric A only) \n
    explicittransform=:auto, \n
    :auto \n
    :none or :shiftinvert, specifying if shift and invert should be explicitly invoked in julia code \n
    tol::Float64=10^(-6), \n
    maxiter::Int64=100, \n
    sigma=0.0 nothing or a number) \n
"""
function EigenProblem(weakformₖ::Function,weakformₘ::Function,test::FESpace,trial::FESpace;
    nev::Int64=10,which::Symbol=:LM,explicittransform::Symbol=:none,tol::Float64=10^(-6),
    maxiter::Int64=100,sigma=0.0)
    L(v) = 0.0
    opK = AffineFEOperator(weakformₖ, L, test, trial)
    opM = AffineFEOperator(weakformₘ, L, test, trial)
    K = opK.op.matrix
    M = opM.op.matrix
    op = EigenOperator(K, M)
    return EigenProblem(trial,test,op,nev,which,explicittransform,tol,maxiter,sigma)
end

export solve
"""
    solve(prob::EigenProblem)
    \n Función para resolver el problema de autovalores
    \n Retorna autovalores y autovectores
"""
function solve(prob::EigenProblem)
    K = prob.op.stima
    M = prob.op.massma
    ξ,Vec = eigs(K,M;nev=prob.nev,which=prob.which,explicittransform=prob.explicittransform,
        tol=prob.tol,maxiter=prob.maxiter,sigma=prob.sigma)
    fₕs = Vector{CellField}(undef, prob.nev)
    Threads.@threads for m=1:prob.nev
        fₕ = FEFunction(prob.trial, Vec[:,m])
        fₕs[m] = fₕ
    end
    return ξ,fₕs
end