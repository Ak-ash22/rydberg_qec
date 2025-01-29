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

    # @assert 1<=site<=num_qubits "Site should be between 1 and num_qubits"
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

function (p1::qubit_parameters)(t)
    """
    Function to get the parameters of the system at time t for the qubits
    such that ti<=t<=tf, where
        - ti=0: start time of the pulse for the qubits A-B-C
        - tf=T1: end time of the pulse for the qubits A-B-C

    Parameters: 
        - time t:: Float64

    Returns: 
        - parameters at time t:: Tuple
    """
    Δt = Δ1_0 - p1.δ * t
    return (Δt, p1.Ω, p1.γ_Decay,p1.γ_dephase, p1.V_nn)
end

function (p2::qubit_parameters)(t)
    """
    Function to get parameters of system at time t for the qubits
    such that ti<=t<=tf, where
        - ti=T1: start time of the pi/2 pulse for the qubits A-B-C
        - tf=T2: end time of the pi/2 pulse for the qubits A-B-C    

    Parameters:
        - time t:: Float64

    Returns:
        - parameters at time t:: Tuple
    """
    Δt = p2.V_nn
    return (Δt, p2.Ω, p2.γ_Decay,p2.γ_dephase, p2.V_nn)
end

function (p3::qubit_parameters)(t)
    """
    Function to get the parameters of the system at time t for the qubits
    such that ti<=t<=tf, where
        - ti=T2: start time of the pulse for the qubits 1-2-3-4-5-6
        - tf=T3: end time of the pulse for the qubits 1-2-3-4-5-6

    Parameters: 
        - time t:: Float64

    Returns: 
        - parameters at time t:: Tuple
    """
    Δt = Δ1_0 - p3.δ * t
    return (Δt, p3.Ω, p3.γ_Decay,p3.γ_dephase, p3.V_nn)
end

function (p4::ancilla_parameters)(t)
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
    Δt = Δ2_0 - p4.δ * t
    return (Δt, p4.Ω, p4.γ_Decay,p4.γ_dephase, p4.V_nn)
end

function hamiltonian1(p1::qubit_parameters,t)
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
    Δ_t, Ω, γ_Decay, γ_dephase, V_nn = p1(t)

    σx_a = full_operator(σ_x, 3,[1])
    σx_c = full_operator(σ_x, 3, [3])

    n_a = full_operator(n,3, [1])
    n_c = full_operator(n, 3, [3])

    nn_ab = full_operator(n, 3, [1,2])
    nn_bc = full_operator(n, 3, [2,3])


    H = Ω/2 .* (σx_a + σx_c) + Δ_t .* (n_a + n_c) + V_nn .* (nn_ab + nn_bc) 

    return H
end
