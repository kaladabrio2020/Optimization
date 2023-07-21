using Plots
using ForwardDiff
using LinearAlgebra
include("Plot_Uni_Pluri//1_Graphics.jl")


function gradient_conjugate( f , x , ftol = 1e-6 , xtol = 1e-6 , maxit =1_000)
    xi = NaN64
    number = 1
    pontos = Vector{Float64}[x]
    

    
    g(k) = reduce(hcat,ForwardDiff.gradient( f , k )) 
    H(k)  = ForwardDiff.hessian( f , k )

    Q  = H(x) 
    d  = -1 * g(x) 

    while true
        
        p = Q * transpose(d)
        R = d * p 
        Gxd = g(x) * transpose(d)

        a = - 1 * ( Gxd[[1]][1] / R[[1]][1] ) 
        a = round( a , digits = 4 )

        xi = x + a * vec(d)
        xi = round.(xi,digits = 4)


        #Restriçao
       
        boolean = [ isnan(i) for i in xi ]
        if  ( any(boolean) )
            break 
        else   
            push!(pontos , xi)
        end

        
        if     ( norm(vec(g(x))) == 0 )                  break 
        elseif ( f(x) - f(xi) <= ftol)                   break
        #elseif ( g(x) - g(xi) <= xtol)                 break
        elseif ( vec(g(xi)) == vec(zeros(length(xi))) ) break 
        end
        

        
        R1 =  g(xi) * p   
        P1 =  d * p
        
        B = R1[[1]][1] / P1[[1]][1] 
        B = round( B, digits = 4 )
        

        di = (- 1 * vec( g(xi) ) ) + B * vec(d)
        di = round.( di , digits = 4 )     

        d  = reduce(hcat,di)
        x  =   xi 
   
        if ( number == maxit ) break end
        number +=1
        
    end
    return xi , pontos , number 
end



f(x) = x[1]^2 + x[2]^2 - x[1]*x[2] - 3*x[1] - 4*x[2] + 1

p    = [0,0]

x, pontos , iteration = gradient_conjugate( f , p )

print(x, iteration )


if (length(p) == 2 ) Graph( f , pontos , p ) end
