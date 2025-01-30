##Helper Functions 

function full_operator(operator::Matrix,num_qubits::Int64,site::Array)
    """

    Function to calculate the matrix of the given operator in N-qubit system acting at 
    the given sites

    Args:
        operator:: Array{Float64,2}: Operator matrix
        num_qubits:: Int64: Total number of qubits in the system
        site: Int64:: Site at which the operator acts

    Returns:
        operator_matrix:: Matrix: Operator acting at the given site

    """

    @assert 1<=num_qubits "Number of qubits should be greater than 0"
    @assert size(operator) == (2,2) "Operator should be a 2x2 matrix"

    I = [1 0; 0 1]
    matrices = [I for _ in 1:num_qubits]

    for s in site
        @assert 1<=s<=num_qubits "Site should be between 1 and num_qubits"
        matrices[s] = operator
    end

    operator_matrix = reduce(kron, matrices)
    return operator_matrix
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
function hamiltonian1(p1::NTuple)
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
    Δ_t, Ω, γ_Decay, γ_dephase, V_nn = p1

    σx_a = full_operator(σ_x, 3,[1])
    σx_c = full_operator(σ_x, 3, [3])

    n_a = full_operator(n,3, [1])
    n_c = full_operator(n, 3, [3])

    nn_ab = full_operator(n, 3, [1,2])
    nn_bc = full_operator(n, 3, [2,3])


    H = Ω/2 .* (σx_a + σx_c) + Δ_t .* (n_a + n_c) + V_nn .* (nn_ab + nn_bc) 

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
    C = [I for _ in 1:length(site)]
    
    for i in site
        C[i] = sqrt(γ_Decay) * full_operator(σ_minus, 3, [i])
    end

    return C
end

function lindbaldian_dephase(γ_dephase::Float64,site::Array)
    """
    Function to calculate the Lindbaldian decay operator for the MCWF method.
    The dephase operators are returned according to given respective sites.

    Args:
        γ_Decay:: Float64: Decay rate
        site:: Array: Array of sites at which the decay operator acts

    Returns:
        C:: Array{Matrix}: Array of decay operators acting on the system
    """
    C = [I for _ in 1:length(site)]
    
    for i in site
        C[i] = sqrt(γ_dephase) * full_operator(σ_z, total_qubits, [i])
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
