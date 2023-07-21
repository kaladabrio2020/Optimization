
function fibonacci(N)
    if (N <= 1) 
        return 1
    else 
        valor = 1
        for i in 1:N  valor += valor end
    end
    return valor
end

function new_numero_iteracao( xtol = 0.3 , b = 2 , eps = 1e-8)
    numero = 1
    div = xtol / b 
    po  = round( (1 + 2 * eps)/div )
    while true
        if ( fibonacci(numero + 1) >= po ) 
            return numero
        else 
            numero+=1 
        end
    end
end

function Method_Fibonacci( f , a , b , xtol = 0.3 , maxit = 1_000 )
    
    N = maxit 
    numero  = 1 
    x1 = x2 = NaN64 

    while true

        p = fibonacci(N) / fibonacci(N+1)

        x1 = a + ( p * ( b - a ))
        x2 = a + ((1-p) * ( b - a ))

        if ( b - a <= xtol )    break end
        if ( numero  == maxit)  break end
        
        if ( f(x1) < f(x2) ) a = x2 end
        if ( f(x1) > f(x2) ) b = x1 end
            
        N -= 1; numero += 1  
    end

    return x1, x2, f(x1), f(x2)  

end

f(x)=x^2+3*x+2
print(Method_Fibonacci(f,-2,2,1e-8,100))