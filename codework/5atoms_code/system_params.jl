const n_qubits = 3
const n_ancillas = 2
const total_qubits = n_qubits + n_ancillas

const g = [1, 0]
const r = [0, 1]

const α = 0                 # Coefficient of |g> state
const β = sqrt(1-α^2)       # Coefficient of |r> state


function wavefunction(num_qubits::Int64, α, β, site::Int64)
    """
    Function to initialize the wavefunction of the N-qubit system with given parameters
    initialized at given site.

    Args:
        num_qubits:: Int64: Total number of qubits in the system
        α:: Any: Coefficient of |g> state
        β:: Any: Coefficient of |r> state
        site:: Int64: Site at which the state is initialized

    Returns:
        ψ:: Array{ComplexF64,1}: Wavefunction of the system
    """

    @assert 1<=site<=num_qubits "Site should be between 1 and $num_qubits"
    @assert 1<=num_qubits "Number of qubits should be greater than 0"
    @assert 0<=α<=1 "α should be between 0 and 1"
    @assert 0<=β<=1 "β should be between 0 and 1"

    site_states = [g for _ in 1:num_qubits]
    site_states[site] = α.*g + β.*r

    ψ = reduce(kron, site_states)
    return ψ
end


function initialize_system()
    "
    System:

    1    3    5
     a1  b3  c5   -----ancillas
      A--B--C    
     a2  b4  c6   -----ancillas
    2    4    6    

    ψ_system = tensor(A,B,C,1,2,3,4,5,6,a1,a2,b3,b4,c5,c6)

    Returns: The wavefunction of the system at t=0.
    "
    ψ_system = wavefunction(total_qubits,α,β,2)

    return ψ_system
end

function params()
    "
    Defining the system parameters
    "
    return Dict(
        :Ω => 1.0,             # Rabi frequency            
        :γ_Decay => 1e-4,        # Decay rate on the qubits
        :γ_dephase => 0.0,      # Dephasing rate on the qubits
        :V1_nn => -1000.0,      # rydberg interaction on the qubits
        # :V2_nn => -1000.0,      # rydberg interaction on the ancillas
        :δ1 => 0.267,            # adiabtatic sweep rate for qubits
        # :δ2 => 0.108,            # adiabtatic sweep rate for ancillas
        :Δ1_0 => 1032.0,        # Detuning at t=0 for qubits
        :Δ2_0 => 23.0,        # Detuning at t=0 for ancillas
        :Δac_0 => 1048.0,         # Detuning at t=0 for Correcting Atom A or Atom c
        :Δb_0 => 2063.0,          # Detuning at t=0 for Correcting Atom B
        :T1 => 240.0,    # Evolution time for step 1 -- From 2atom_optimal_decay_result.ipynb
        :T2 => 55.0,     # Evolution time for step 2
        :T3 => 242.0,     # Evolution time for Correcting Atom A or Atom C
        :T4 => 242.0     # Evolution time for Correcting Atom B
    )
end

function unpack_params()  #Need to define these as const in main.jl
    "
    Unpack the system parameters
    "
    p = params()
    return p[:Ω], p[:γ_Decay], p[:γ_dephase], p[:V1_nn], p[:δ1], p[:Δ1_0], p[:Δ2_0], p[:Δac_0], p[:Δb_0], p[:T1], p[:T2], p[:T3], p[:T4]
end
