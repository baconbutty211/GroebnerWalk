## Julia
The following code is written in Julia, which only runs on linux, instructions for installing can be found here: https://julialang.org/install/
### Oscar
The code also requires Oscar, which is a julia package, instructions for installing can be found here: https://www.oscar-system.org/install/

## X-homogeneous Groebner walk algorithms
The x-homogeneous groebner walk algorithms witness, lift, flip can be found in the XHomoGWalk.jl file along with their optimised counterparts.
In GroebnerWalk.jl you can find an example (Example 5.6 from https://arxiv.org/abs/1512.02662) that tests flip, flip_optimised and therefore lift, lift_optimised, witness, witness_optimised.
Oscar does not provide a method for intially reducing a groebner basis, so the file initialReduction.jl along with some helper functions in functions.jl have been provided to implement algorithm 4.7 in https://arxiv.org/abs/1512.02662. 

## Experiments
experiment.jl is some code not relevant to the x-homogeneous groebner walk algorithms. It simply compares the lexicographical and degree reverse lexicographical groebner bases for the same ideal.
