"""
    The max eigenvalue of A'A, which is also the max eigenvalue of AA' 
    is estimated via the Power Iteration Method
"""
function Power_Iteration(A)

    
    N,M = size(A)
    
    if N > M 
        H = A'*A;
        b = randn(M,1)
    else
        H = A*A';
        b = randn(N,1)
    end
        e = 1.
        
        for k = 1:20
            tmp = H*b
            e = norm(tmp)
            tmp = tmp/e
            b = tmp
        end
        
        return e
    end


        

