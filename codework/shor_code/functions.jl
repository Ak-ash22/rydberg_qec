##Helper Functions 
function full_operator(gate, total_qubits, sites)
    """
    Applies an arbitrary gate on a specified site `i` in a `total_qubits`-qubit system.
    All other sites are identity operators.

    Args:
    - gate: AbstractOperator (arbitrary gate to apply on site `i`)
    - i: Array (site index to apply the gate, 1-based)
    - total_qubits: Int (total number of qubits)

    Returns:
    - operator: AbstractOperator (the full operator acting on the entire system)
    """
    # Ensure the gate is an AbstractOperator
    if !(gate isa AbstractOperator)
        throw(ArgumentError("The gate must be an AbstractOperator"))
    end

    # Create an identity operator for each qubit
    identity = transition(NLevelBasis(2), 1,1) + transition(NLevelBasis(2), 2,2)
    identity_ops = [identity for _ in 1:total_qubits]
    
    # Replace the identity operator at site `i` with the provided gate
    for j in sites
        identity_ops[j] = gate
    end

    # Return the Kronecker product of all operators
    return tensor(identity_ops...)
end

struct qubit_parameters
    Ω::Float64
    γ_Decay::Float64
    γ_dephase::Float64
    V_nn::Float64
    δ::Float64
end

struct ancilla_parameters
    Ω::Float64
    γ_Decay::Float64
    γ_dephase::Float64
    V_nn::Float64
    δ::Float64
end


function get_qubit_parameters(p::qubit_parameters,t::Float64,mode::Symbol)
    """
    Function to get the parameters of the system at time t for the qubits
   
    Modes:
        - `:T1` → 0<=t<=T1: Pulse for the qubits A-B-C
        - `:T2` → T1<=t<=T2: Pi/2 pulse for the qubits A-B-C
        - `:T3` → T2<=t<=T3: Pulse for the qubits 1-2-3-4-5-6

    Parameters: 
        - p:: qubit_parameters → Struct containing the qubit parameters
        - time t:: Float64
        - mode:: Symbol → Selects phase (`:T1` or `:T2` or `:T3`)

    Returns: 
        - parameters at time t:: Tuple
    """
    if mode == :T1
        Δt = Δ1_0 - p.δ * t
    elseif mode == :T2
        Δt = p.V_nn
    elseif mode == :T3
        Δt = Δ1_0 - p.δ * t
    else
        throw(ArgumentError("Invalid mode! Choose from `:T1`, `:T2`, `:T3`"))
    end
    return (Δt, p.Ω, p.γ_Decay,p.γ_dephase, p.V_nn)
end

function get_ancilla_parameters(p::ancilla_parameters,t::Float64)
    """
    Function to get the parameters of the system at time t for the ancillas
    such that ti<=t<=tf, where
        - ti=T3: start time of the pulse for the ancillas
        - tf=T4: end time of the pulse for the ancillas

    Parameters: 
        - time t:: Float64

    Returns: 
        - parameters at time t:: Tuple
    """
    Δt = Δ2_0 - p.δ * t
    return (Δt, p.Ω, p.γ_Decay,p.γ_dephase, p.V_nn)
end

#Hamiltonians
function hamiltonian1(p::NTuple)
    """
    Function to calculate the Hamiltonian for the driving the qubits A-B-C at time t such that
    ti<=t<=tf, where
        - ti=0: start time of the pulse for the qubits A-B-C
        - tf=T1: end time of the pulse for the qubits A-B-C

    Parameters:
        - p:: qubit_parameters
        - t:: Float64

    Returns:
        - H:: Matrix : Hamiltonian at time t
    """
    Δ_t, Ω, γ_Decay, γ_dephase, V_nn = p
    
    basis = NLevelBasis(2)
    n = transition(basis,2,2)
    σx = transition(basis,1,2) + transition(basis,2,1)

    σx_a = full_operator(σx, total_qubits, [1])
    σx_c = full_operator(σx, total_qubits, [3])

    n_a = full_operator(n,total_qubits, [1])
    n_c = full_operator(n, total_qubits, [3])

    nn_ab = full_operator(n, total_qubits, [1,2])
    nn_bc = full_operator(n, total_qubits, [2,3])


    H = sparse(Ω/2 .* (σx_a + σx_c) .+ Δ_t .* (n_a + n_c) .+ V_nn .* (nn_ab + nn_bc))

    return H
end

#Lindbald Operators
function lindbaldian_decay(γ_Decay::Float64,site::Array)
    """
    Function to calculate the Lindbaldian decay operator for the MCWF method.
    The decay operators are returned according to given respective sites.

    Args:
        γ_Decay:: Float64: Decay rate
        site:: Array: Array of sites at which the decay operator acts

    Returns:
        C:: Array{Matrix}: Array of decay operators acting on the system
    """
    basis = NLevelBasis(2)
    σ_minus = transition(basis,1,2)
    identity = transition(basis,1,1) + transition(basis,2,2)
    C = Vector{Operator}(undef, total_qubits)
    
    for i in 1:total_qubits
        if i in site
            C[i] = sqrt(γ_Decay) .* full_operator(σ_minus, total_qubits, [i])
        else
            C[i] = full_operator(identity, total_qubits, [i])
        end
    end

    return C
end

#Helper function for mcwf_dynamic
function f(t,ψ)
    p = qubit_parameters(Ω,γ_Decay,γ_dephase,V_nn,δ)
    pt = get_qubit_parameters(p,t,:T1)

    H = hamiltonian1(pt)
    C = lindbaldian_decay(pt[3],[1,2,3])
    Cdagger = [adjoint(i) for i in C]
    return H, C, Cdagger
end