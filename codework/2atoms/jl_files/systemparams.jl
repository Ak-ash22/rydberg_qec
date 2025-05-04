using LinearAlgebra
using DifferentialEquations
using JLD2, FileIO

const n_atoms = 2
const r = [0.0;1.0]
const g = [1.0;0.0]
#initial state of the system
const α = 1/sqrt(2)                 # Coefficient of |r> state
const β = sqrt(1-α^2)               # Coefficient of |g> state


function initialize_system()
    "
    System:
    1-----A------2
    
    Output: The initial density matrix of the system
    "    
    ψ_a = α.*r + β.*g
    
    ψ_0 = [g,ψ_a,g]
    ρ_0 = complex(reduce(kron,ψ_0) * reduce(kron,ψ_0)');
    
    return ρ_0
end

function system_constants()
    "
    Defining the system Matrices
    "
    #2x2 Matrices
    σ_x = [0 1; 1 0]
    n = [0 0; 0 1]
    Π_g = [1 0; 0 0]
    n = [0 0; 0 1]
    I = [1 0; 0 1]
    σ_minus = [0 1; 0 0]
    σ_plus = [0 0; 1 0]
    σ_z = [1 0; 0 -1]
    
    return σ_x, n, Π_g, n, I, σ_minus, σ_plus, σ_z
end

function params()
    "
    Defining the system parameters
    "
    return Dict(
        :Ω1 => 1.0,             # Rabi frequency
        :Ω2 => 1.0,             
        :γ_Decay => 0.0,        # Decay rate
        :γ_dephase => 0.0,      # Dephasing rate
        :V1_nn => -1000.0,      # rydberg interaction
        :V2_nn => -1000.0,
        :δ => 0.0,             # adiabtatic sweep rate
        :Δ1_0 => 1032.0,        # initial detuning
        :Δ2_0 => 1032.0,
        :T_optimal => 0.0     # optimal time 
    )
end

function unpack_params()
    "
    Unpack the system parameters
    "
    p = params()
    return p[:Ω1], p[:Ω2], p[:γ_Decay], p[:γ_dephase], p[:V1_nn], p[:V2_nn], p[:δ], p[:Δ1_0], p[:Δ2_0], p[:T_optimal]
end
