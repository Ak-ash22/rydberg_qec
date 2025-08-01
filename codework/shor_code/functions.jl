include("system_params.jl")
include("dependencies.jl")
const Ω, γ_Decay, γ_dephase, V_nn, Δ_0, T1 = unpack_params()

##############################################################################################  Helper Functions 
const local_basis = NLevelBasis(3)

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
    identity = transition(NLevelBasis(3), 1,1) + transition(NLevelBasis(3), 2,2) + transition(NLevelBasis(3), 3,3)
    identity = Operator(identity.basis_l, identity.basis_r, SparseMatrixCSC{ComplexF32, Int64}(identity.data))
    identity_ops = [identity for _ in 1:total_qubits]
    
    # Replace the identity operator at site `i` with the provided gate
    for j in sites
        identity_ops[j] = gate
    end

    # Return the Kronecker product of all operators
    return tensor(reverse(identity_ops)...)

    # @assert issorted(sites) "sites must be given in ascending order"

    # cb = tensor(reverse(ntuple(_ -> local_basis, total_qubits))...)
    # ops = ntuple(_ -> gate, length(sites))
    # return LazyTensor(cb, sites, ops)
end


struct qubit_parameters
    γ_Decay::Float64
    γ_dephase::Float64
    V_nn::Float64
    Δ_0::Float64
end


function get_qubit_parameters(p::qubit_parameters,t::Float64,mode::Symbol)
    """
    Function to get the parameters of the system at time t for the qubits
   
    Modes:
        - `:T1` → 0<=t<=T1: Pulse for the qubits A-B-C
        - `:T2` → 0<=t<=T2: Pulse for the ancillas 1-2
        - `:T3` → 0<=t<=T3: Pulse for the qubits A or C
        - `:T4` → 0<=t<=T3: Pulse for the qubit B

    Parameters: 
        - p:: qubit_parameters → Struct containing the qubit parameters
        - time t:: Float64
        - mode:: Symbol → Selects phase (`:T1` or `:T2` or `:T3` or `:T4`)

    Returns: 
        - parameters at time t:: Tuple
    """

    @assert mode in [:T1, :T2, :T3, :T4a, :T4b, :T5] "Invalid mode selected"

    #### Qubit Encoding Modes
    if mode == :T1

        if t <= π
            Ω1 = 1.0
            Ω2 = 0.0
        elseif t < (π+sqrt(2)π)
            Ω1 = 0.0
            Ω2 = 1.0
        else 
            Ω1 = 1.0
            Ω2 = 0.0
        end

        return [Ω1/2, Ω2/2, Ω2/2, p.Δ_0, p.V_nn, Ω2/2, Ω2/2, p.Δ_0, p.V_nn]
    end
end

function get_full_wavefunction(psi::Array, sites::Vector)

    site_states = [a for _ in 1:total_qubits]
    for j in sites
        site_states[j] = psi[j]
    end
    ψ_system = reduce(kron, site_states)

    return ψ_system
end












################################################################################################ Lindbaldian Operators
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
    basis = NLevelBasis(3)
    σ_minus = transition(basis,1,2)
    σ_minus = Operator(σ_minus.basis_l, σ_minus.basis_r, SparseMatrixCSC{ComplexF32, Int64}(σ_minus.data))

    identity = transition(basis,1,1) + transition(basis,2,2) + transition(basis,3,3)
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

function lindbaldian_dephase(γ_dephase::Float64,total_qubits::Int64,site::Array)
    """
    Function to calculate the Lindbaldian dephase operator for the MCWF method.
    The dephase operators are returned according to given respective sites.

    Args:
        γ_dephase:: Float64: Dephase rate
        site:: Array: Array of sites at which the decay operator acts

    Returns:
        C:: Array{Matrix}: Array of dephase operators acting on the system
    """
    basis = NLevelBasis(3)
    σ_z= transition(basis,1,1) - transition(basis,2,2) 
    σ_z = Operator(σ_z.basis_l, σ_z.basis_r, SparseMatrixCSC{ComplexF32, Int64}(σ_z.data))

    identity = transition(basis,1,1) + transition(basis,2,2) + transition(basis,3,3)
    identity = Operator(identity.basis_l, identity.basis_r, SparseMatrixCSC{ComplexF32, Int64}(identity.data))
    C = Vector{AbstractOperator}(undef, total_qubits)
    
    for i in 1:total_qubits
        if i in site
            C[i] = sqrt(γ_dephase) * full_operator(σ_z, total_qubits, [i])
        else
            C[i] = full_operator(identity, total_qubits, [i])
        end
    end

    return C
end

# Matrix Constants to define the Hamiltonians and Lindbaldian operators
basis = NLevelBasis(3)
n_r = transition(basis, 3, 3)
n_r = Operator(n_r.basis_l, n_r.basis_r, SparseMatrixCSC{ComplexF32, Int64}(n_r.data))

σx_1r = transition(basis, 2, 3) + transition(basis, 3, 2)
σx_1r = Operator(σx_1r.basis_l, σx_1r.basis_r, SparseMatrixCSC{ComplexF32, Int64}(σx_1r.data))

σx_0r = transition(basis, 1, 3) + transition(basis, 3, 1)
σx_0r = Operator(σx_0r.basis_l, σx_0r.basis_r, SparseMatrixCSC{ComplexF32, Int64}(σx_0r.data))

# Parameters for the Atoms -- coefficients for the Hamiltonian
p = qubit_parameters(γ_Decay, γ_dephase, V_nn, Δ_0)

const C = lindbaldian_dephase(γ_dephase,total_qubits,[i for i in 1:total_qubits])
const Cdagger = [adjoint(c) for c in C]

################################################################################################ Step 1: Encoding A-B-C
#Timespan for driving atoms A-B-C
const tspan1 = [0.0:0.1:T1;]

#####System Hamiltonian 1 - driving atoms A-C
## Driving Atom A
σx_0r_atom1 = full_operator(σx_0r, total_qubits, [1])
σx_1r_atom1 = full_operator(σx_1r, total_qubits, [1])
const n_r_atom1 = full_operator(n_r, total_qubits, [1])
nn_r12 = full_operator(n_r, total_qubits, [1,2])

#Driving Atom B
σx_0r_atom2 = full_operator(σx_0r, total_qubits, [2])
σx_1r_atom2 = full_operator(σx_1r, total_qubits, [2])
const n_r_atom2 = full_operator(n_r, total_qubits, [2])

##Driving Atom C
σx_0r_atom3 = full_operator(σx_0r, total_qubits, [3])
σx_1r_atom3 = full_operator(σx_1r, total_qubits, [3])
const n_r_atom3 = full_operator(n_r, total_qubits, [3])
nn_r23 = full_operator(n_r, total_qubits, [2,3])

const coeff1 = [t->get_qubit_parameters(p,t,:T1)]
const H1 = LazySum([coeff1[1](tspan1[1])[i] for i ∈ 1:9],[σx_1r_atom2, σx_0r_atom1, σx_1r_atom1, n_r_atom1, nn_r12, σx_0r_atom3, σx_1r_atom3, n_r_atom3, nn_r23])

function Ht1(t)
    """
    Function to calculate the time dependent Hamiltonian for the MCWF method till time steps T1+T2.
        -- H1: Hamiltonian for driving atoms A-B-C for time 0:T1
        -- H2: Hamiltonian for applying hadamards on atoms A-B-C for time T1:T2
        -- H3: Hamiltonian for driving atoms 1-2 for time T2:T3
    
    Args:
        t:: Float64: Time
    
    Returns:
        H:: LazySum: Time dependent Hamiltonian
    """

    coeffs = coeff1[1](t)
    for i in eachindex(coeffs)
        H1.factors[i] = coeffs[i]
    end
    return H1

end


function f1(t,ψ)
"""
Function to calculate the time evolution of the system using the MCWF method.

Args:
    t:: Float64: Time

Returns:
    H:: LazySum: Time dependent Hamiltonian
    C:: Array{Operator}: Array of decay operators acting on the system
    Cdagger:: Array{Operator}: Array of adjoint decay operators acting on the system
"""

    H = Ht1(t)

    return H, C, Cdagger
end