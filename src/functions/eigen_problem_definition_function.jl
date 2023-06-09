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
function eigen_problem(weakformₖ::Function,weakformₘ::Function,test::FESpace,trial::FESpace;
    nev::Int64=10,which::Symbol=:LM,explicittransform::Symbol=:none,tol::Float64=10^(-6),
    maxiter::Int64=100,sigma=0.0)
    # source vector (always need to be zero for eigen problems)
    F(v) = 0.0; 
    opH = AffineFEOperator(weakformₖ, F, test, trial)
    opE = AffineFEOperator(weakformₘ, F, test, trial)
    H = opH.op.matrix
    E = opE.op.matrix
    op = EigenOperator(H,E)
    return EigenProblem(trial,test,op,nev,which,explicittransform,tol,maxiter,sigma)
end

export solve
"""
    solve(prob::EigenProblem)
    \n Función para resolver el problema de autovalores
    \n Retorna autovalores y autovectores
"""
function solve(prob::EigenProblem)
    H = prob.op.hamiltonian
    E = prob.op.energy
    ϵ,eigenvecs = eigs( H,E;nev=prob.nev,which=prob.which,explicittransform=prob.explicittransform,
        tol=prob.tol,maxiter=prob.maxiter,sigma=prob.sigma)
    ϕ = Vector{CellField}(undef, prob.nev)
    Threads.@threads for m=1:prob.nev
        ϕₙ = FEFunction(prob.trial, eigenvecs[:,m])
        ϕ[m] = ϕₙ
    end
    return ϵ,ϕ
end