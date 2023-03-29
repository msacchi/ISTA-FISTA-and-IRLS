"""
    Estimation of sparse solution x of the problem

    x = argmin ||A x- y||_2^2 + λ ||x||_1

    using the FAST Iterative Soft-Thresholding Algorithm (FISTA)

    M D Sacchi
    March 2023
"""

function FISTA(A,y,Niter,λ)

    
    Soft(x,alpha) = sign(x)*max(abs(x)-alpha, 0)

    N,M = size(A)
    e = Power_Iteration(A)
    η = 0.95/e
   
    # Start FISTA
    
    x = zeros(Float64,M)
    t = 1. 
    p = x
        
        J = zeros(Niter)
    for k = 1:Niter
            x_old = x
            x = p .-  η*A'*(A*p.-y)
            x = Soft.(x, η*λ)
            t_old=t
            t = 0.5*(1+sqrt(1+t*t^2))
            p = x +((t_old-1)/t)*(x-x_old)
            J[k] = 0.5*sum((A*x-y).^2) + λ*sum(abs.(x))
    end
        return x, J
end
        
    
