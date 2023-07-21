using Plots
using ForwardDiff
using LinearAlgebra
include("Plot_Uni_Pluri//OneDimensional.jl")
include("Plot_Uni_Pluri//Graphics1.jl")


function MethodQuasiNewton( f , x  , ftol = 1e-8 , gtol = 1e-8 , maxit = 1_000)
    xi = NaN64
    pontos = Vector{Float64}[x] 
    number = 1
    h(k) = ForwardDiff.hessian( f , k )
    g(k) = reduce(hcat,ForwardDiff.gradient( f , k ))

    H = h(x)
    
    #if (! minimum(H) >= 0) return [NaN64 for i in 1:3] end
    
    d = -H*g(x)'

    while true
        p(a)  = f(vec(x) + a * vec(d))
        alpha = Method_Newton( p , 10 )
        
        xi  = x + alpha * d
        #xi  = round.(xi,5)
      
        push!(pontos,vec(xi))
        Dxi =  alpha * d'
        

        Dgx = g(xi) - g(x)
        
        p1 = (( (Dxi * Dxi') )/ (Dxi * Dgx'))[:,1]
        
        z1 = (( H * Dgx' )' * ( H * Dgx' ))
    
        z2 = (Dgx * (H * Dgx' ))
        
        p2 = (vec(z1)./vec(z2)) 

        Hxi = H .+ ( p1 - p2 )

        H = Hxi ; x = xi

        if ( number == maxit ) break end
        number += 1

        if     ( norm(vec(g(x))) == 0 )                  break 
        elseif ( f(x) - f(xi) <= ftol)                   break
        end
    end
    
    if (length(x)==2) Graph( f , pontos , x ) end

    return f(xi),xi,number,pontos
end
