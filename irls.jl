

"""
    Estimation of sparse solution x of the problem

    x = argmin ||A x- y||_2^2 + λ ||x||_1

    using Iteratively Reweighed Least Squares (IRLS) 

    M D Sacchi
    March 2023
"""

    function IRLS(A,y,Niter,λ)

    x = zeros(Float64,M)
  
    J = zeros(Niter)
    G = A'*A 
    b = A'*y 
     for k in 1:Niter
        q = 1.0./(abs.(x).+0.0001)
       	Q = diagm(0 => q)
        x = vec((G + λ*Q)\b)
        J[k] = 0.5*sum((A*x-y).^2) + λ*sum(abs.(x))
     end
        return x, J
end

        

