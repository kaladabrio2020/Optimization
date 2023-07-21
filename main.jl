# Pnl 
include("NonLinear/Multidimesionalidade/Gradient.jl")
include("NonLinear/Multidimesionalidade/GradientConjugate.jl")
include("NonLinear/Multidimesionalidade/Newton.jl")
include("NonLinear/Multidimesionalidade/QuasiNewton.jl")

# PL
include("Linear//Dual//Dual.jl")
include("Linear//Primal//Primal.jl")

a = [ 1 ,2 ,3, 4]
c = [-4 -5 0 0]
A = [-1 -4 1 0;
     -3 -2 0 1]

    
a = [i for i =1:size(A)[2]]

b = [4;6;18]

MethodSimplex(
    A,
    [1,2,5],
    [4,3],
    c,
    b,
    a
)