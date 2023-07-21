using Plots
using ForwardDiff
using LinearAlgebra
include("Plot_Uni_Pluri//Graphics1.jl")

function method_newton( f , x , xtol = 1e-8 , ftol = 1e-8 , maxit = 1_000 )
    xi = NaN64
    number = 1
    pontos = Vector{Float64}[x]

    Vf(x) = reduce( hcat , ForwardDiff.gradient( f , x ) )
    H(x)  = inv(ForwardDiff.hessian( f , x ))

    while true

        d  = H(x) * transpose(Vf(x))  
        xi = x - vec(d)

        boolean = [ isnan(i) for i in xi ]
        if     ( any(boolean) )          break 
        elseif ( det(H(x))   == 0 )      break 
        elseif ( f(x) - f(xi) <= ftol )  break 
        elseif ( norm(x - xi) <= xtol )  break 
        elseif ( number == maxit )       break  
        end

        push!( pontos , xi )
        number += 1 ; x = xi
    
    end
    if (length(x) == 2 ) Graph( f , pontos , x ) end

    return f(xi), xi , pontos , number
end
