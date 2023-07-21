using Plots
using ForwardDiff
using LinearAlgebra
include("Plot_Uni_Pluri//0_One_Dimensional.jl")
include("Plot_Uni_Pluri//1_Graphics.jl")



function GradientProjection( f , x , A , b  , ftol=1e-5 , gtol=1e-5 ,maxit=1_000 )
    d = xi = NaN64
    pontos = Vector{Float64}[x]
    number = 1

    g(x)   = reduce(hcat,ForwardDiff.gradient( f , x ))
    solve = A' * (inv(A*A') * A)
    
    i  = Matrix(I , size(solve)[1] , size(solve)[2])
    pk = i - solve
    
    
    while true
        
        d =  - 1* (pk * g(x)')
      
        h(a) = f( x + a * vec(d) )
        
        alpha = Method_Newton( h , 10 )
    
        xi = x + alpha * vec(d)
        push!(pontos,xi)
        R = b - A * ( x + pk*d )
        print(xi)
        if ( norm(R) != 0 )         break 
        elseif ( number == maxit )  break      

        end
        
        number += 1 
        x = xi
        
    end
    return  f(xi) , xi , pontos ,number
    
end





# Z
f(x) = x[1]^2 + 4*x[2]^2 - 8*x[1] - 16*x[2]

p = [5.0,0.0]
A = [3 -2;3 2]
b = [6;15]

fx , x , pontos ,iter = GradientProjection(f , p , A , b )


if (length(p)== 2) 
    Graph(f , pontos , A , b )  
end