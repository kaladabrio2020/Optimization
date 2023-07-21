using Plots
using ForwardDiff
using LinearAlgebra
include("Plot_Uni_Pluri//0_One_Dimensional.jl")
include("Plot_Uni_Pluri//1_Graphics.jl")


function Method_Gradient( f , x , gtol= 1e-8 ,ftol=1e-8 , maxit=1_000)
    xi = NaN64 
    number = 1
    pontos = Vector{Float64}[x]
    g(k) = ForwardDiff.gradient( f , k )


    while true
        d = - g(x)

        p(a)  = f(x - a * g(x)) 
        alpha = Method_Newton( p , 1 )

        xi = x + alpha * d  
        xi = round.(xi,digits=5)
        #_______________________________
        println(xi)
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
    return xi ,  pontos , number
end

f(x)   = x[1]^2 + x[2]^2 - x[1]*x[2] - 3*x[1] - 4*x[2] + 1
inicio = [2,3.2]

x , pontos , iter = Method_Gradient( f , inicio )
print(x , iter)

if (length(inicio) == 2 ) Graph(f,pontos,inicio) end 