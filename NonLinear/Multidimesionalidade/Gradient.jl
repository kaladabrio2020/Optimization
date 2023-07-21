using Plots
using ForwardDiff
using LinearAlgebra
include("Plot_Uni_Pluri//OneDimensional.jl")
include("Plot_Uni_Pluri//Graphics1.jl")


function MethodGradient( f , x , gtol= 1e-8 ,ftol=1e-8 , maxit=1_000)
    xi = NaN64 
    number = 1
    pontos = Vector{Float64}[x]
    g(k)   = ForwardDiff.gradient( f , k )


    while true
        d = - g(x)

        p(a)  = f(x - a * g(x)) 
        alpha = Method_Newton( p , 1 )

        xi = x + alpha * d  
        xi = round.(xi,digits=5)
        #_______________________________
        cond1 = norm(x - xi)
        cond2 = abs(f(x) - f(xi))
        if ( any(isnan.(xi)) )     break end  
        if     ( cond1 < gtol   )  break
        elseif ( cond2 <= ftol   ) break
        elseif ( number == maxit ) break
        end
        
        push!(pontos, xi)
        number += 1
        x = xi
    end
    
    if (length(x) == 2 ) Graph( f , pontos , x ) end 

    return f(xi ), xi ,  pontos , number
end

