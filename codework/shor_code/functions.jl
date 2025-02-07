include("system_params.jl")
include("dependencies.jl")
const Ω, γ_Decay, γ_dephase, V_nn, δ, Δ1_0, T_optimal = unpack_params()

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
    identity = Operator(identity.basis_l, identity.basis_r, SparseMatrixCSC{ComplexF32, Int64}(identity.data))
    identity_ops = [identity for _ in 1:total_qubits]
    
    # Replace the identity operator at site `i` with the provided gate
    for j in sites
        identity_ops[j] = gate
    end

    # Return the Kronecker product of all operators
    return tensor(reverse(identity_ops)...)
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


function get_qubit_parameters(p::qubit_parameters,t::Float64)
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

    Δt = Δ1_0 - p.δ * t

    return [p.Ω, p.Ω, Δt, Δt, p.V_nn, p.V_nn]
    # return [p.Ω, p.Ω]
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
    σ_minus = Operator(σ_minus.basis_l, σ_minus.basis_r, SparseMatrixCSC{ComplexF32, Int64}(σ_minus.data))

    identity = transition(basis,1,1) + transition(basis,2,2)
    identity = Operator(identity.basis_l, identity.basis_r, SparseMatrixCSC{ComplexF32, Int64}(identity.data))
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


basis = NLevelBasis(2)
n = transition(basis,2,2)
n = Operator(n.basis_l, n.basis_r, SparseMatrixCSC{ComplexF32, Int64}(n.data))

σx = transition(basis,1,2) + transition(basis,2,1)
σx = Operator(σx.basis_l, σx.basis_r, SparseMatrixCSC{ComplexF32, Int64}(σx.data))

σx_a = full_operator(σx, total_qubits, [1])
σx_b = full_operator(σx, total_qubits, [2])
σx_c = full_operator(σx, total_qubits, [3])

n_a = full_operator(n,total_qubits, [1])
n_c = full_operator(n, total_qubits, [3])

nn_ab = full_operator(n, total_qubits, [1,2])
nn_bc = full_operator(n, total_qubits, [2,3])

p = qubit_parameters(Ω,γ_Decay,γ_dephase,V_nn,δ)

const coeff = [t->get_qubit_parameters(p,t)]
const tspan = [0.0:0.1:T_optimal;]
const H = LazySum([coeff[1](tspan[1])[i] for i ∈ 1:6],[σx_a, σx_c, n_a, n_c, nn_ab, nn_bc])
# const H = LazySum([coeff[1](tspan[1])[i] for i ∈ 1:2],[σx_a, σx_c])


function Ht(t)
    for i=1:length(H.factors)
        H.factors[i] = coeff[1](t)[i]
    end
    return H
end

#Helper function for mcwf_dynamic
const C = lindbaldian_decay(1e-3,[1,3])
const Cdagger = [adjoint(c) for c in C]

function f(t,ψ)
    H = Ht(t)
    return H, C, Cdagger
end

