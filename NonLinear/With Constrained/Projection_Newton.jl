using ForwardDiff
using LinearAlgebra


function Projection_Newton( f , x , A , b , ftol = 1e-6 , xtol = 1e-6 ,maxit = 1_000)
    xi = NaN64
    number = 1

    Vf(x) = reduce( hcat , ForwardDiff.gradient( f , x ) )
    H(x)  = inv(ForwardDiff.hessian( f , x ))

    solve = A' * ((A*A')\A)
    i  = Matrix(I , size(solve)[1] , size(solve)[2])
    pk = i - solve

    while true

        d  = -(pk * (H(x) *Vf(x)')) 
        xi = x + vec(d)

        R = b - A * ( x + pk*d )
        if     ( norm(R) != 0 )          break 
        elseif ( det(H(x))   == 0 )      break 
        elseif ( f(x) - f(xi) <= ftol )  break 
        elseif ( norm(x - xi) <= xtol )  break 
        elseif ( number == maxit )       break  
        end

  
        number += 1 ; x = xi
    
    end
    return f(xi) , xi , number
end


f(x) = x[1]^2 + x[2]^2 - x[1]*x[2] - 3*x[1] - 4*x[2] + 1


p    = [-1,1]

A = [1 2]
b = [3;]




fx , x ,  iter = Projection_Newton(f,p,A,b)

println( fx , x , iter )
