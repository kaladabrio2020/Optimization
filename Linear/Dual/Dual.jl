using LinearAlgebra

function nova_base( M , Bi , Ni , xb)

    index_xb  = findall( isequal( minimum( xb)), xb )[1]
    index_max = findall( isequal( maximum( M )), M  )[1]
    numero    = Ni[index_max]
    bi = Bi[index_xb]
    deleteat!(Bi , index_xb )
    insert!(Bi , index_xb , numero )

    deleteat!(Ni , index_max)
    insert!(Ni , index_max , bi)
    
    return Bi , Ni

end

function Dual_Simplex(A , Bi , Ni , c , b , a , sense = false ,MaxIter = 1) #sense = true(min) e false(max)
    lista = NaN64
    iter = 1    # Iteraçao
    B = A[:,Bi] # B

    while true
        zc   = []
        Yjn  = []
        bool = []
        if ( iter == MaxIter ) break end
        #Verificando se o determinante e diferente de zero
        if ( det(B) == 0     ) break end 

        Binv = inv(B) # Fazendo inversa da Matriz B
        Binv = round.(Binv,digits=4)

        xb   = Binv * b #xb
        xb   = round.(xb,digits=3)


        #Verificando se e primal viavel
        if ( all(vec(xb).>=0 ))  
            append!(bool,true)
            println("Primal viavel " , xb,"\n")
        

        else  
            append!(bool,false)                   
            println("Nao e primal viavel ",xb ,"\nB = ",Bi ,"N = ",Ni)
        
        end
        

        for i in Ni
            zj_cj = (c[:,Bi] * Binv) * A[:,i] - c[:,i]
            append!( zc , zj_cj[1] )
            
        end
      

        if (all(zc.>= 0) ) append!(bool,true )
        else               append!(bool,false)
        end

        if (all(bool.== true) ) println("Solucao primal e dual viavel");break end

        index = findall( isequal(minimum(xb)), xb )
        for i in Ni
            Y = Binv[index,:] * A[:,i]
            append!( Yjn , Y )
        end
        Bi , Ni = nova_base((zc./Yjn) , Bi , Ni , xb )
        B = A[:,Bi]  
        iter += 1
    
    end
    return false , false , false , false   
end


a = [ 1 ,2 ,3, 4]
c = [4 5 0 0]
A = [1 4 1 0;
     3 2 0 1]
b = [5;-7]
B = [3,4]
N = [1,2]
x = Dual_Simplex(A ,B , N , c , b , a )