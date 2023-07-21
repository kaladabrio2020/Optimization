using ForwardDiff

function number_iteration( a , b , xtol=0.3 )
    return round( (log((b - a)/xtol))/log(2) )
end

function bissection_search( f , a , b , xtol = 1e-10 , maxit = 1_000)
    number = 1 ; df = 0.0
    xi( a , b ) = ( b + a )/2

    while true

        df = ForwardDiff.derivative( f , xi( a , b ))

        if ( number == maxit ) break
        end
        
        if ( b - a  <= xtol  ) break
        end
        
        if ( df > 0 ) b = xi( a , b ) 
        else          a = xi( a , b )
        end
        
        number += 1

    end

    return xi( a , b )

end

f(x) =2*x^2 - 2*x + 8

a = 40
b = 90
xtol = 0.3
iter = number_iteration( a , b )

x  = bissection_search(f , a  , b ,  xtol , iter)

