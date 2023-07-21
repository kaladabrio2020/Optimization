using LinearAlgebra
using Combinatorics

using Combinatorics

function MethodSimplex(A , Bi , Ni , c , b , a , sense = true )
    iter = 1    # Iteraçao
    B = A[:,Bi] # B
    print(B)

    
    x = []
    lista = NaN64
    #Verificando se o determinante e diferente de zero
    if ( det(B) == 0 ) return false,0,0 end 
    Binv = inv(B) 
    Binv = round.(Binv,digits=3)
        

    xb   = Binv * b #xb
    xb   = round.(xb,digits=3)

    #Verificando se solucao viavel
    if ( !all(vec(xb).>=0 )) 
        return false , 0 , 0
    end


    #Z - C = cb * inv(B) * ai - ci
    cb = c[:,Bi]        
        
    for i in Ni    
        z_c = cb * Binv * A[:,i] - c[:,i]
        append!(x, z_c[1] )
    end
        

    #Verificando se e solucao otima Z - C <= 0 para minimizacao
    if ( all(x.<=0 ) & sense )  
        for e in sort(Ni) lista = insert!(xb , e , 0) end
        return lista , Bi , Ni

    #Verificando se e solucao otima Z - C >= 0 para maximizacao
    elseif ( all(x.>= 0 ) & !sense )
        for e in sort(Ni) lista = insert!(xb , e , 0) end
        return lista , Bi , Ni

    else 
        return false , 0 , 0
    end
    return false , 0 , 0

end

## DEPOIS FACO
function Permutacao(A,b)
    for Bi in per
  
        prime = []
        boolean , Be , Ne = MethodSimplex( A , Bi , setdiff(a , Bi) , c , b , a , false )
        
        if ( boolean != false )
            push!(xi, boolean)
            push!(prime, Be) 
            push!(prime, Ne)  
            append!(bases, [prime])
    
        end 
    end
    
    
    print("====PRIMAL====")
    for i in 1:length(xi)
        f(x) = 3*x[1] + 5*x[2]
        println()
        println("Base")
        println("B =",bases[i][1] , "N = ", bases[i][2])
        println("x = ",xi[i] ,"      f(x) = ", f(xi[i]) ,"\n")
    
    end
    
end
