const n_qubits = 9
const n_ancillas = 8
const total_qubits = n_qubits + n_ancillas

const a = ComplexF64[1,0,0]
const b = ComplexF64[0,1,0]
const r = ComplexF64[0,0,1]

const α = 0.0                 # Coefficient of |g> state
const β = sqrt(1-α^2)       # Coefficient of |r> state
const ϕ = 0.0


function wavefunction(num_qubits::Int64, α, β, ϕ::Float64, site::Vector)
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

    # @assert 1<=site<=num_qubits "Site should be between 1 and $num_qubits"
    @assert 1<=num_qubits "Number of qubits should be greater than 0"
    @assert 0.0<=α<=1.0 "α should be between 0 and 1"
    @assert 0.0<=β<=1.0 "β should be between 0 and 1"

    site_states = [a for _ in 1:num_qubits]
    for j in site
        site_states[j] = α .* b + (exp(1im * ϕ) * β) .* a
    end
    ψ = reduce(kron, site_states)
    return ψ
end


function initialize_system()
    "
    System:

    1    3    5
     a1   b3  c5   -----ancillas
      A-p1-B-p2-C    
     a2   b4  c6   -----ancillas
    2    4    6    

    ψ_system = tensor(A,B,C,1,2,3,4,5,6,a1,a2,b3,b4,c5,c6)

    Returns: The wavefunction of the system at t=0.
    "
    ψ_system = wavefunction(total_qubits,α,β,ϕ,[2])

    return ψ_system
end

function params()
    "
    Defining the system parameters
    "
    return Dict(
        :Ω => 1.0,             # Rabi frequency            
        :γ_Decay => 0.0,        # Decay rate on the qubits
        :γ_dephase => 0.0,      # Dephasing rate on the qubits
        :V_nn => -1000.0,      # rydberg interaction on the qubits
        :Δ_0 => 1000.0,        # Detuning at t=0 for encoding qubits and ancillas (CNOT gate)
        :T1 => (π + sqrt(2)π + π),    # Evolution time for step 1 -- Encoding A-B-C
        # :T2 => pi/2,     # Evolution time for step 2
        # :T3 => 154.0     # Evolution time for step 3
        # :T4 => 582.0     # Evolution time for step 4
        # :T5 => 582.0     # Evolution time for step 5  -- for phase flip error syndrome
        # :T6 => 582.0     # Evolution time for step 6  -- for phase flip error syndrome
    )
end

function unpack_params()  #Need to define these as const in main.jl
    "
    Unpack the system parameters
    "
    p = params()
    return p[:Ω], p[:γ_Decay], p[:γ_dephase], p[:V_nn], p[:Δ_0], p[:T1]
end
