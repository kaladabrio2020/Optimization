using Plots

function Graph( f , pontos , inicio )

    Plots.GRBackend()
    gr( size=( 800,600 ) )
        
    lx = [ i[1] for i in pontos ] 
    ly = [ i[2] for i in pontos ]
    
        
    texto = "  [" * 
                string( round( lx[ length(lx) ] , digits = 3 ) ) * "  ,  "*
                string( round( ly[ length(ly )] , digits = 3 ) ) * "]  "
    
    
    uniao = union(lx, ly)
    a = -2 * abs( minimum(minimum.(uniao)) )-2
    b =  2 + abs( maximum(maximum.(uniao)) )+1
    x = y = range( a , stop = b , length = 50 )
    
    
    p = contourf( x , y , ( x , y )-> f([x;y]), leg=false , c =:turbid )
        
    plot!( lx , ly , c=:red  )
    scatter!( lx , ly , markersize = 4 , markercolor=:white )
        
    annotate!( lx[length(lx)], ly[length(ly)] , text( texto ,:right , :green , 9 ) )
    display(p)
end
