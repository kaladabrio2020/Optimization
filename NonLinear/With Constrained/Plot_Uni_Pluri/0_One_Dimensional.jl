using ForwardDiff

function Method_Newton( f , xi = 0.5 , xtol = 1e-8 , ftol = 1e-8 , maxit = 1_000 )
    number = 1 
    fxi = fci = NaN64

    while true
        d_1(x) = ForwardDiff.derivative( f , x )
        d_2(x) = ForwardDiff.derivative( s -> ForwardDiff.derivative(f, s), x)

        c = xi - ( d_1(xi) / d_2(xi) )

        fxi = f(xi)  
        fci = f(c) 
        
        if ( xi  == 0)   xi  = 1 end
        if ( fxi == 0)   fxi = 1 end
            

        condicao1 = ( abs( xi - c ) / abs(xi) ) 
        condicao2 = ( abs( fxi - fci )/ abs( fxi ) )

        if ( condicao1 <= xtol ) break end
        if ( condicao2 <= ftol ) break end
        if ( number   == maxit ) break end
      
        number+=1
        xi = c  

    end
    return xi
end


