include("system_params.jl")
include("dependencies.jl")
const Ω, γ_Decay, γ_dephase, V_nn, δ, Δ1_0, Δ2_0, Δac_0, Δb_0, T1, T2, T3, T4 = unpack_params()

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

    @assert mode in [:T1, :T2, :T3, :T4] "Invalid mode selected"

    #Qubit Encoding Modes
    if mode == :T1
        Δt = Δ1_0 - p.δ * t
        return [p.Ω/2, p.Ω/2, Δt, Δt, p.V_nn, p.V_nn]

    elseif mode == :T2
        Δt = Δ2_0 - p.δ * t
        return [p.Ω/2, p.Ω/2, Δt, Δt, p.V_nn/(2^6), p.V_nn/(2^6), p.V_nn/(2^6), p.V_nn/(2^6)]
    
    #Error Correction Modes
    elseif mode == :T3
        Δt = Δac_0 - p.δ * t
        return [p.Ω/2, Δt, p.V_nn, p.V_nn/(2^6)]

    elseif mode == :T4
        Δt = Δb_0 - p.δ * t
        return [p.Ω/2, Δt, p.V_nn, p.V_nn, p.V_nn/(2^6), p.V_nn/(2^6)]

    end

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

# Matrix Constants to define the Hamiltonians and Lindbaldian operators
basis = NLevelBasis(2)
n = transition(basis,2,2)
n = Operator(n.basis_l, n.basis_r, SparseMatrixCSC{ComplexF32, Int64}(n.data))
σx = transition(basis,1,2) + transition(basis,2,1)
σx = Operator(σx.basis_l, σx.basis_r, SparseMatrixCSC{ComplexF32, Int64}(σx.data))

# Parameters for the Atoms -- coefficients for the Hamiltonian
p = qubit_parameters(Ω,γ_Decay,γ_dephase,V_nn,δ)

#Timespan for driving atoms A-B-C and 1-2-3
const tspan = [0.0:0.1:(T1+T2);]

# System Hamiltonian 1 - driving atoms A-B-C
σx_a = full_operator(σx, total_qubits, [1])
σx_b = full_operator(σx, total_qubits, [2])
σx_c = full_operator(σx, total_qubits, [3])
const n_a = full_operator(n, total_qubits, [1])
const n_b = full_operator(n, total_qubits, [2])
const n_c = full_operator(n, total_qubits, [3])
nn_ab = full_operator(n, total_qubits, [1,2])
nn_bc = full_operator(n, total_qubits, [2,3])
const coeff1 = [t->get_qubit_parameters(p,t,:T1)]
const H1 = LazySum([coeff1[1](tspan[1])[i] for i ∈ 1:6],[σx_a, σx_c, n_a, n_c, nn_ab, nn_bc])

#System Hamiltonian 2 - driving atoms 1-2
σx_1 = full_operator(σx, total_qubits, [4])
σx_2 = full_operator(σx, total_qubits, [5])
const n_1 = full_operator(n, total_qubits, [4])
const n_2 = full_operator(n, total_qubits, [5])
nn_a1 = full_operator(n, total_qubits, [1,4])
nn_b1 = full_operator(n, total_qubits, [2,4])
nn_b2 = full_operator(n, total_qubits, [2,5])
nn_c2 = full_operator(n, total_qubits, [3,5])
const coeff2 = [t->get_qubit_parameters(p,t,:T2)]
const H2 = LazySum([coeff2[1](tspan[1])[i] for i ∈ 1:8],[σx_1, σx_2, n_1, n_2, nn_a1, nn_b1, nn_b2, nn_c2])


function Ht(t)
"""
Function to calculate the time dependent Hamiltonian for the MCWF method till time steps T1+T2.
    -- H1: Hamiltonian for driving atoms A-B-C for time 0:T1
    -- H2: Hamiltonian for driving atoms 1-2-3 for time T1:T1+T2

Args:
    t:: Float64: Time

Returns:
    H:: LazySum: Time dependent Hamiltonian
"""

    if t<T1 || t==T1
        coeffs = coeff1[1](t)
        for i in eachindex(coeffs)
            H1.factors[i] = coeffs[i]
        end
        return H1

    elseif t<(T1+T2) || t==(T1+T2)
        coeffs = coeff2[1](t-T1)
        for i in eachindex(coeffs)
            H2.factors[i] = coeffs[i]
        end
        return H2
    end
end


#Helper function for mcwf_dynamic
const C = lindbaldian_decay(1e-3,[1,2,3,4,5])    
const Cdagger = [adjoint(c) for c in C]


function Ct(t)
"""
Function to calculate the time dependent Lindbaldian decay operators for the MCWF method till time steps T1+T2.

Args:
    t:: Float64: Time

Returns:
    C:: Array{Operator}: Array of decay operators acting on the system
    Cdagger:: Array{Operator}: Array of adjoint decay operators acting on the system
"""

    return C, Cdagger
end

function f(t,ψ)
"""
Function to calculate the time evolution of the system using the MCWF method.

Args:
    t:: Float64: Time

Returns:
    H:: LazySum: Time dependent Hamiltonian
    C:: Array{Operator}: Array of decay operators acting on the system
    Cdagger:: Array{Operator}: Array of adjoint decay operators acting on the system
"""

    H = Ht(t)
    return H, Ct(t)...
end


## Error Correction
const tspan2 = [0.0:0.1:T3;]  #Time span for error correction of atom A or C
const tspan3 = [0.0:0.1:T4;]  #Time span for error correction of atom B

#Required Matrix Constants
n_abc = full_operator(n, total_qubits, [1,2,3])
n_abc = Operator(n_abc.basis_l, n_abc.basis_r, SparseMatrixCSC{ComplexF32, Int64}(n_abc.data))

#Hamiltonian for Error Correction of Atom A
const coeff3 = [t->get_qubit_parameters(p,t,:T3)]
const H_correct_a = LazySum([coeff3[1](tspan2[1])[i] for i ∈ 1:4],[σx_a, n_a, nn_ab, nn_a1])

#Hamiltonian for Error Correction of Atom B
const coeff4 = [t->get_qubit_parameters(p,t,:T4)]
const H_correct_b = LazySum([coeff4[1](tspan3[1])[i] for i ∈ 1:6],[σx_b, n_b, nn_ab, nn_bc, nn_b1, nn_b2])

#Hamiltonian for Error Correction of Atom C
const H_correct_c = LazySum([coeff3[1](tspan2[1])[i] for i ∈ 1:4],[σx_c, n_c, nn_bc, nn_c2])

function Ht_correct(t,site)
"""
Function to calculate the time dependent Hamiltonian for the MCWF method for error correction of atoms A-B-C.
    - H_correct_a: Hamiltonian for error correction of atom A
    - H_correct_b: Hamiltonian for error correction of atom B
    - H_correct_c: Hamiltonian for error correction of atom C

Args:
    t:: Float64: Time
    site:: Int: Site of the atom to be corrected

Returns:
    H:: LazySum: Time dependent Hamiltonian
"""

    if site == 1
        coeffs = coeff3[1](t)
        for i in eachindex(coeffs)
            H_correct_a.factors[i] = coeffs[i]
        end
        return H_correct_a

    elseif site == 2
        coeffs = coeff4[1](t)
        for i in eachindex(coeffs)
            H_correct_b.factors[i] = coeffs[i]
        end
        return H_correct_b

    elseif site == 3
        coeffs = coeff3[1](t)
        for i in eachindex(coeffs)
            H_correct_c.factors[i] = coeffs[i]
        end
        return H_correct_c
    end
end

function Ct_correct(t)
"""
Function to calculate the time dependent Lindbaldian decay operators for the MCWF method for error correction of atoms A-B-C.

Args:
    t:: Float64: Time

Returns:
    C:: Array{Operator}: Array of decay operators acting on the system
    Cdagger:: Array{Operator}: Array of adjoint decay operators acting on the system
"""

    return C, Cdagger
end


function f_correct(t,ψ,site)
"""
Function to calculate the time evolution of the system using the MCWF method for error correction.
    
Args:
    t:: Float64: Time
    ψ:: Array: State vector of the system
    site:: Int: Site of the atom to be corrected

Returns:
    H:: LazySum: Time dependent Hamiltonian
    C:: Array{Operator}: Array of decay operators acting on the system
    Cdagger:: Array{Operator}: Array of adjoint decay operators acting on the system
"""

    H = Ht_correct(t,site)
    return H, Ct(t)...
end

function f_correct_factory(site)
    """
    Function to create a closure for the f_correct function with a fixed site.
    
    Args:
    - site: Int (site of the atom to be corrected)
    
    Returns:
    - f_correct_site: Function (f_correct with the site argument fixed)
    """
    return (t,ψ) -> f_correct(t, ψ, site)
end

