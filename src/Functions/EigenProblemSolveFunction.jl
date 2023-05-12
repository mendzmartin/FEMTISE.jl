export EigenValuesAndEigenVectors
"""
    params=(nev,tol,maxiter,explicittransform,sigma)
"""
function EigenValuesAndEigenVectors(p,q,r,dΩ,USpace,VSpace;params=(10,10e-9,100,:none,1.0))
    a,b=bilineal_forms(p,q,r,dΩ);   # Define bilinear forms and FE spaces
    # solve eigenvalue problem
    prob=EigenProblem(a,b,USpace,VSpace;nev=params[1],tol=params[2],maxiter=params[3],explicittransform=params[4],sigma=params[5]);
    ϵ,ϕ=solve(prob);
    return ϵ,ϕ;
end