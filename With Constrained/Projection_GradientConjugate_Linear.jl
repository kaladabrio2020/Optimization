using ForwardDiff
using LinearAlgebra
include("Plot_Uni_Pluri//1_Graphics.jl")

function Projection_GradientConjugate(f , x , A , b , ftol = 1e-8 , xtol = 1e-8 , maxit = 1_000)
    xi = NaN64
    number = 1
    pontos = Vector{Float64}[x]
    
    g(k) = reduce(hcat,ForwardDiff.gradient( f , k )) 
    H(k)  = ForwardDiff.hessian( f , k )
    

    solve = A' * ((A*A')\A)
    i  = Matrix(I , size(solve)[1] , size(solve)[2])
    pk = i - solve
    
    Q  = H(x) 
    B  = d = 0.0 

    while true
        
        d  = (pk * ( -g(x) .+ (B * d))')'   
        p = Q * transpose(d)
        R = d * p 
        Gxd = g(x) * transpose(d)

        a = - 1 * ( Gxd[[1]][1] / R[[1]][1] ) 
        a = round( a , digits = 4 )

        xi = x + a * vec(d)
        xi = round.(xi,digits = 4)

    

        #Restriçao residuo
        
        R = b - A * ( x + pk*d' )
        if ( norm(R) != 0 ) break end

        #Demais Restriçao
        boolean = [ isnan(i) for i in xi ]
        if   ( any(boolean)   ) break 
        else push!(pontos , xi)
        end

        if     ( pk*g(x)' == 0)                         break
        elseif ( number == maxit )                      break
        elseif ( f(x) - f(xi) <= ftol)                  break
        elseif ( norm(vec(g(xi))) == 0 )                break 
        elseif ( norm(vec(g(x)) - vec(g(xi))) <= xtol)  break
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
        
        number +=1
        
    end
    return f(xi), xi , pontos , number 
end

f(x) = x[1]^2 + 4*x[2]^2 - 8*x[1] - 16*x[2]

p    = [-2,0]
A = [3 -2; 2 3]
b = [6;4]




fx , x , pontos , iter = Projection_GradientConjugate(f,p,A,b)

println( fx , x , pontos , iter )
if (length(p) == 2) Graph(f,pontos,A,b) end

