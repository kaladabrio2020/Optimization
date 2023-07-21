using Plots

function Graph( f , pontos , A , bi)
    Plots.GRBackend()
    gr( size=( 800,600 ) )
    
    lx = [ i[1] for i in pontos ] 
    ly = [ i[2] for i in pontos ]
   
    

    uniao = union(lx, ly)
    a = -2 * abs( minimum(minimum.(uniao)) )-2
    b =  2 + abs( maximum(maximum.(uniao)) )+1
    xi = yi = range( a , stop = b , length = 50 )
    
    
    contour( xi , yi , ( xi , yi )-> f([xi;yi]), leg=false , #=c =:turbid=# )
        
    plot!( lx , ly , c=:red  )
    scatter!( lx , ly , markersize = 4 , markercolor=:white )
    
    for i=1:size(A)[1]
        p(x) = (- A[i,1]*x + (bi[i]))/A[i,2]
        display(plot!(x -> p(x),xlim=(a,b),ylim=(a,b),color=:green))
    end
    
end
