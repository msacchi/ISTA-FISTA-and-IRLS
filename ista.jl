
"""
    Estimation of sparse solution x of the problem 
	
    x = argmin ||A x- y||_2^2 + λ ||x||_1  

    using the Iterative Soft-Thresholding Algorithm (ISTA) 

    M D Sacchi 
    March 2023
"""

function ISTA(A,y,Niter,λ)

    
    Soft(x,alpha) = sign(x)*max(abs(x)-alpha, 0)

    N,M = size(A)
    e = Power_Iteration(A)
    η = 0.95/e
   
   # Start ISTA
    
    x = zeros(Float64,M)
  
    J = zeros(Niter)
    for k = 1:Niter
        u = x .-  η*A'*(A*x.-y)
        x = Soft.(u, η*λ)
        J[k] = 0.5*sum((A*x-y).^2) + λ*sum(abs.(x))
    end
        return x, J
end

    
