include("system_params.jl")
include("dependencies.jl")
const Ω1, Ω2, γ_Decay, γ_dephase, V_nn, Δ_0, Δb_0, T1, T2, T3, T4, T5 = unpack_params()

##Helper Functions 
function full_operator(gate, qubits, sites)
    """
    Applies an arbitrary gate on a specified site `i` in a `qubits`-qubit system.
    All other sites are identity operators.

    Args:
    - gate: AbstractOperator (arbitrary gate to apply on site `i`)
    - i: Array (site index to apply the gate, 1-based)
    - qubits: Int (total number of qubits)

    Returns:
    - operator: AbstractOperator (the full operator acting on the entire system)
    """
    # Ensure the gate is an AbstractOperator
    if !(gate isa AbstractOperator)
        throw(ArgumentError("The gate must be an AbstractOperator"))
    end

    # Create an identity operator for each qubit
    identity = transition(NLevelBasis(3), 1, 1) + transition(NLevelBasis(3), 2, 2) + transition(NLevelBasis(3), 3, 3)
    identity = Operator(identity.basis_l, identity.basis_r, SparseMatrixCSC{ComplexF32, Int64}(identity.data))
    identity_ops = [identity for _ in 1:qubits]
    
    # Replace the identity operator at site `i` with the provided gate
    for j in sites
        identity_ops[j] = gate
    end

    # Return the Kronecker product of all operators
    return tensor(reverse(identity_ops)...)
end


struct qubit_parameters
    γ_Decay::Float64
    γ_dephase::Float64
    V_nn::Float64
    Δ_0::Float64
    Δb_0::Float64
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

    #### Hadamard Mode 1    
    elseif mode == :T2

        Ω = 1.0
        Δ = 100
        return [Ω, Ω, Ω, Ω, Ω, Ω, Δ, Δ, Δ]
        
    
    #### Ancilla driving mode
    elseif mode == :T3
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
        return [Ω1/2, Ω1/2, Ω1/2, Ω2/2, Ω2/2, Ω2/2, Ω2/2, p.Δ_0, p.Δ_0, p.V_nn, p.V_nn, p.V_nn, p.V_nn]
    
    #### Correction mode for Atom A or C   
    elseif mode == :T4a
        if t <= π
            Ω1 = 1.0
            Ω2 = 0.0
            Ω3 = 0.0
        elseif t < (π+sqrt(2)π)
            Ω1 = 0.0
            Ω2 = 1.0
            Ω3 = 0.0
        else 
            Ω1 = 0.0
            Ω2 = 0.0
            Ω3 = 1.0
        end
        return [Ω1/2, Ω3/2, Ω2/2, Ω2/2, p.Δ_0, p.V_nn]

    #### Correction mode for Atom B
    elseif mode == :T4b
        if t <= π
            Ω1 = 1.0
            Ω2 = 0.0
            Ω3 = 0.0
        elseif t < (π+sqrt(2)π)
            Ω1 = 0.0
            Ω2 = 1.0
            Ω3 = 0.0
        else 
            Ω1 = 0.0
            Ω2 = 0.0
            Ω3 = 1.0
        end
        return [Ω1/2, Ω1/2, Ω3/2, Ω3/2, Ω2/2, Ω2/2, 2*p.Δ_0, p.V_nn, p.V_nn]

    #### Hadamard Mode 2
    elseif mode == :T5

        Ω = 1.0
        Δ = 100
        return [Ω, Ω, Ω, Ω, Ω, Ω, Δ, Δ, Δ]

    end

end

#Lindbald Operators
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
    C = Vector{Operator}(undef, total_qubits)
    
    for i in 1:total_qubits
        if i in site
            C[i] = sqrt(γ_dephase) .* full_operator(σ_z, total_qubits, [i])
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
p = qubit_parameters(γ_Decay, γ_dephase, V_nn, Δ_0, Δb_0)



######################################################################################################## Encoding Atoms - Step 1
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


######################################################################################################## Applying hadamards - Step 2

virtual_z_single_gate = exp(-1im * π/4) * transition(basis, 1, 1) + exp(1im * π/4) * transition(basis, 2, 2) + transition(basis, 3, 3)
virtual_z_single_gate = Operator(virtual_z_single_gate.basis_l, virtual_z_single_gate.basis_r, SparseMatrixCSC{ComplexF32, Int64}(virtual_z_single_gate.data))
const virtual_z_full = full_operator(virtual_z_single_gate, total_qubits, [1,2,3])

const tspan2 = [0.0: 0.1: T2;]
const coeff2 = [t->get_qubit_parameters(p,t,:T2)]    
const H2 = LazySum([coeff2[1](tspan2[1])[i] for i ∈ 1:9], [σx_0r_atom1, σx_1r_atom1, σx_0r_atom2, σx_1r_atom2, σx_0r_atom3, σx_1r_atom3, n_r_atom1, n_r_atom2, n_r_atom3])

function Ht2(t)
    """
    Function to calculate the time dependent Hamiltonian for the MCWF method from time steps T1 to T2.
        -- H2: Hamiltonian for applying hadamards on atoms A-B-C for time T1:T2
    Args:
        t:: Float64: Time
    Returns:
        H:: LazySum: Time dependent Hamiltonian
    """

    coeffs = coeff2[1](t)
    for i in eachindex(coeffs)
        H2.factors[i] = coeffs[i]
    end
    return H2

end


const C_encoding = lindbaldian_dephase(γ_dephase,total_qubits,[i for i in 1:total_qubits])    
const Cdagger_encoding = [adjoint(c) for c in C_encoding]

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

    return H, C_encoding, Cdagger_encoding
end

function f2(t,ψ)
"""
Function to calculate the time evolution of the system using the MCWF method.
Args:
    t:: Float64: Time
Returns:
    H:: LazySum: Time dependent Hamiltonian
    C:: Array{Operator}: Array of decay operators acting on the system
    Cdagger:: Array{Operator}: Array of adjoint decay operators acting on the system
"""

    H = Ht2(t)

    return H, C_encoding, Cdagger_encoding 
end

# ######################################################################################################## Driving atoms 1-2 - Step 3
const tspan3 = [0.0:0.1:T3;]  #Time span for driving atoms 1-2

### Driving operators for the ancillas
σx_0r_ancilla1 = full_operator(σx_0r, total_qubits, [4])
σx_1r_ancilla1 = full_operator(σx_1r, total_qubits, [4])
σx_0r_ancilla2 = full_operator(σx_0r, total_qubits, [5])
σx_1r_ancilla2 = full_operator(σx_1r, total_qubits, [5])

### Detuning operators for the ancillas
n_r_ancilla1 = full_operator(n_r, total_qubits, [4])
n_r_ancilla2 = full_operator(n_r, total_qubits, [5])

### Rydberg Interaction operators
nn_r14 = full_operator(n_r, total_qubits, [1,4])
nn_r24 = full_operator(n_r, total_qubits, [2,4])
nn_r25 = full_operator(n_r, total_qubits, [2,5])
nn_r35 = full_operator(n_r, total_qubits, [3,5])


const coeff3 = [t->get_qubit_parameters(p,t,:T3)]
const H3 = LazySum([coeff3[1](tspan3[1])[i] for i ∈ 1:13],[σx_1r_atom1, σx_1r_atom2, σx_1r_atom3, σx_0r_ancilla1, σx_1r_ancilla1, σx_0r_ancilla2,
                    σx_1r_ancilla2, n_r_ancilla1, n_r_ancilla2, nn_r14, nn_r24, nn_r25, nn_r35])

const C = lindbaldian_dephase(γ_dephase,total_qubits,[i for i in 1:total_qubits])
const Cdagger = [adjoint(c) for c in C]

function Ht3(t)
"""
Function to calculate the time dependent Hamiltonian for the MCWF method for driving atoms 1-2.
    - H3: Hamiltonian for driving atoms 1-2 for time T2:T3
Args:
    t:: Float64: Time
Returns:
    H:: LazySum: Time dependent Hamiltonian
"""

    coeffs = coeff3[1](t)
    for i in eachindex(coeffs)
        H3.factors[i] = coeffs[i]
    end
    return H3

end

function f3(t,ψ)
"""
Function to calculate the time evolution of the system using the MCWF method for driving atoms 1-2.
Args:
    t:: Float64: Time
Returns:
    H:: LazySum: Time dependent Hamiltonian
    C:: Array{Operator}: Array of decay operators acting on the system
    Cdagger:: Array{Operator}: Array of adjoint decay operators acting on the system
"""

    H = Ht3(t)
    return H, C, Cdagger
end



# ######################################################################################################## Error Correction - Step 4
const tspan4 = [0.0:0.1:T4;]  #Time span for error correction of atom A or C or B


#Hamiltonian for Error Correction of Atom A
const coeff4 = [t->get_qubit_parameters(p,t,:T4a)]
const H_correct_a = LazySum([coeff4[1](tspan4[1])[i] for i ∈ 1:6],[σx_1r_ancilla1, σx_0r_ancilla1, σx_0r_atom1, σx_1r_atom1, n_r_atom1, nn_r14])

#Hamiltonian for Error Correction of Atom B
const coeff5 = [t->get_qubit_parameters(p,t,:T4b)]
const H_correct_b = LazySum([coeff5[1](tspan4[1])[i] for i ∈ 1:9],[σx_1r_ancilla1, σx_1r_ancilla2, σx_0r_ancilla1, σx_0r_ancilla2, σx_0r_atom2, σx_1r_atom2, n_r_atom2, nn_r24, nn_r25])

#Hamiltonian for Error Correction of Atom C
const H_correct_c = LazySum([coeff4[1](tspan4[1])[i] for i ∈ 1:6],[σx_1r_ancilla2, σx_0r_ancilla2, σx_0r_atom3, σx_1r_atom3, n_r_atom3, nn_r35])

#Hamiltonian for No Correction -- Zero Hamiltonian
const H_no_correct = LazySum([0.0],[σx_0r_atom1])

function Ht_correct(t,site)
"""
Function to calculate the time dependent Hamiltonian for the MCWF method for error correction of atoms A-B-C.
    - H_correct_a: Hamiltonian for error correction of atom A
    - H_correct_b: Hamiltonian for error correction of atom B
    - H_correct_c: Hamiltonian for error correction of atom C
    - H_no_correct: Zero Hamiltonian -- for no correction

Args:
    t:: Float64: Time
    site:: Int: Site of the atom to be corrected

Returns:
    H:: LazySum: Time dependent Hamiltonian
"""

    if site == 1
        coeffs = coeff4[1](t)
        for i in eachindex(coeffs)
            H_correct_a.factors[i] = coeffs[i]
        end
        return H_correct_a

    elseif site == 2
        coeffs = coeff5[1](t)
        for i in eachindex(coeffs)
            H_correct_b.factors[i] = coeffs[i]
        end
        return H_correct_b

    elseif site == 3
        coeffs = coeff4[1](t)
        for i in eachindex(coeffs)
            H_correct_c.factors[i] = coeffs[i]
        end
        return H_correct_c

    elseif site == 0
        return H_no_correct
    end
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
    return H, C, Cdagger
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



# ###############################################################################################Timespan for applying hadamards - Step 5
const tspan5 = [0.0: 0.1: T5;]

# const coeff5 = [t->get_qubit_parameters(p,t,:T5)]  
const H5 = LazySum([coeff2[1](tspan5[1])[i] for i ∈ 1:9], [σx_0r_atom1, σx_1r_atom1, σx_0r_atom2, σx_1r_atom2, σx_0r_atom3, σx_1r_atom3, n_r_atom1, n_r_atom2, n_r_atom3])

function Ht5(t)
    """
    Function to calculate the time dependent Hamiltonian for the MCWF method from time steps T1 to T2.
        -- H2: Hamiltonian for applying hadamards on atoms A-B-C for time T1:T2
    Args:
        t:: Float64: Time
    Returns:
        H:: LazySum: Time dependent Hamiltonian
    """

    coeffs = coeff2[1](t)
    for i in eachindex(coeffs)
        H5.factors[i] = coeffs[i]
    end
    return H5

end


function f5(t,ψ)
    """
    Function to calculate the time evolution of the system using the MCWF method for driving atoms 1-2.
    Args:
        t:: Float64: Time
    Returns:
        H:: LazySum: Time dependent Hamiltonian
        C:: Array{Operator}: Array of decay operators acting on the system
        Cdagger:: Array{Operator}: Array of adjoint decay operators acting on the system
    """
    
        H = Ht5(t)
        return H, C, Cdagger
    end
    

# ############################################################################################################### Functions to compute the dynamical phase

# struct PrecomputedOps{T}
#     σminus_sys::Vector{Matrix{T}}
#     σplus_sys::Vector{Matrix{T}}
#     σx_sys::Vector{Matrix{T}}
#     σz_sys::Vector{Matrix{T}}
#     n_sys::Vector{Matrix{T}}
#     nn_sys::Vector{Matrix{T}}
# end

# function PrecomputedOps(n_atoms::Int=3)
#     σ_x     = [0 1; 1 0]
#     σ_minus = [0 1; 0 0]
#     σ_plus  = [0 0; 1 0]
#     σ_z     = [1 0; 0 -1]
#     n       = [0 0; 0 1]
#     I₂      = Matrix{Float64}(I,2,2)

#     make(op) = [reduce(kron, [i==k ? op : I₂ for i=1:n_atoms]) for k=1:n_atoms]

#     return PrecomputedOps(
#         make(σ_minus), make(σ_plus), make(σ_x),
#         make(σ_z),     make(n),      [kron(n,n,I₂), kron(I₂,n,n)]
#     )
# end

# const OPS = PrecomputedOps()


# p_tuple = (Ω,γ_Decay,γ_dephase,V_nn,δ1,Δ1_0,Δ2_0,OPS)


# function energy_level_spaghetti(time::Array,p)
#     """
#     Function to calculate the energy levels of the system at different time steps.
#     Args:
#         time:: Array: Array of time steps
#         p:: Tuple: Parameters of the system
#     Returns:
#         - e:: Array: Array of energy levels
#         - e_vec:: Array: Array of eigenvectors
#         - Delta:: Array: Array of detuning values
#     """
#     # Unpack the solution object
#     t_vals = time  
#     e = []
#     e_vec = []
#     Delta = []

#     for t in t_vals

#         #parameters
#         Ω, γ_Decay, γ_dephase, V_nn, δ, Δ1_0, Δ2_0, ops = p
#         Δ_t = Δ1_0 - δ*t

#         # unpack the pre‑built operators ↓
#         σx_sys, σminus_sys, σplus_sys,
#         σz_sys, n_sys, nn_sys = ops.σx_sys, ops.σminus_sys, ops.σplus_sys,
#                                 ops.σz_sys, ops.n_sys, ops.nn_sys

#         #Hamiltonian
#         H = Ω/2 .* σx_sys[1] + Ω/2 .* σx_sys[3] + V_nn .* nn_sys[1] + V_nn.* nn_sys[2] +
#             Δ_t .* n_sys[1] + Δ_t .* n_sys[3]

#         eigenvals = eigen(H,sortby=nothing).values
#         eigvecs = eigen(H,sortby=nothing).vectors
#         push!(e_vec, eigvecs)
#         push!(e, eigenvals)
#         push!(Delta, Δ_t)
#     end

#     return e, e_vec, Delta
# end


# function trapezoidal_integrate(t,E1,E2)
#     """
#     Function to compute the dynamical phase using the trapezoidal rule for numerical integration.
#     Args:
#         t:: Array: Array of time steps
#         E1:: Array: Array of energy levels for the first state
#         E2:: Array: Array of energy levels for the second state
#     Returns:
#         - Δ0_dyn:: Float64: Dynamical phase
#         - ΔEs:: Array: Array of energy differences
#     """
#         # Trapezoidal rule for numerical integration
#         ΔEs = E2 .- E1
#         Δ0_dyn = -2*π* sum(diff(t) .* (ΔEs[1:end-1] .+ ΔEs[2:end]) / 2)
#         return Δ0_dyn, ΔEs
# end


# function compute_dynamical_phase(t::Array)
#     """
#     Function to compute the dynamical phase of the system.
#     Args:
#         t:: Array: Array of time steps
#         E1:: Array: Array of energy levels for the first state
#         E2:: Array: Array of energy levels for the second state

#     Returns:
#         - Δ0_dyn:: Float64: Dynamical phase
#     """

#     # Calculate the energy levels
#     energy_levels, _, delta = energy_level_spaghetti(t,p_tuple)

#     #ggg energy level
#     e_ggg = [subarray[4] for subarray in energy_levels]

#     #rrr energy level
#     e_rrr = [subarray[1] for subarray in energy_levels]

#     Δ0_dyn, _ = trapezoidal_integrate(t,e_ggg,e_rrr)
#     Δ0_dyn = mod(-Δ0_dyn,2π)

#     return Δ0_dyn,e_ggg,e_rrr, delta
# end

# function apply_dynamical_phase(Δ0_dyn::Float64, ψ_obtained::Ket)
#     """
#     Function to apply the dynamical phase correction to the system.
#     Args:
#         - Δ0_dyn:: Float64: Dynamical phase

#     Returns:
#         - ψ_target:: Array: Target state after applying the dynamical phase correction
#     """

#     # Rz_minus = Operator(NLevelBasis(2), [exp(-im * Δ0_dyn) 0; 0 exp(im * Δ0_dyn)])
#     # Rz_correction = full_operator(Rz_minus, num_qubits, [1,3])
#     # ψ_corrected = Rz_correction * (ψ_obtained/norm(ψ_obtained))

#     # ψ_obtained /= norm(ψ_obtained)
#     # ψ_obtained.data[end] = exp(im * Δ0_dyn) * ψ_obtained.data[end]
#     # ψ_corrected = ψ_obtained/norm(ψ_obtained)

#     ggg = reduce(kron, [g, g, g])
#     rrr = reduce(kron, [r, r, r])
#     amp_ggg = ψ_obtained.data' * ggg
#     amp_rrr = ψ_obtained.data' * rrr
#     ϕ_actual = angle(amp_rrr) - angle(amp_ggg)
#     ϕ_actual = mod(ϕ_actual, 2π)

#     ψ_obtained.data[end] = exp(im * ϕ_actual) * ψ_obtained.data[end]
#     ψ_corrected = ψ_obtained/norm(ψ_obtained)
    
#     println("Actual dynamical phase: ", ϕ_actual)
#     # println("Calculated dynamical phase: ", Δ0_dyn)
#     return ψ_corrected
# end


########################################################################################################## Metrics required for Plots visualization
n0 = transition(NLevelBasis(3),1,1)
n1 = transition(NLevelBasis(3),2,2)

const n0_atom1 = full_operator(n0,total_qubits,[1])
const n0_atom2 = full_operator(n0,total_qubits,[2])
const n0_atom3 = full_operator(n0,total_qubits,[3])
const n0_ancilla1 = full_operator(n0,total_qubits,[4])
const n0_ancilla2 = full_operator(n0,total_qubits,[5])

const n1_atom1 = full_operator(n1,total_qubits,[1])
const n1_atom2 = full_operator(n1,total_qubits,[2])
const n1_atom3 = full_operator(n1,total_qubits,[3])
const n1_ancilla1 = full_operator(n1,total_qubits,[4])
const n1_ancilla2 = full_operator(n1,total_qubits,[5])

#Stabilizer Generators
σx = transition(NLevelBasis(3),1,2) + transition(NLevelBasis(3),2,1)
const S1 = full_operator(σx,total_qubits,[1,2])
const S2 = full_operator(σx,total_qubits,[2,3])