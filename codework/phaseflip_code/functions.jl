include("system_params.jl")
include("dependencies.jl")
const Ω, γ_Decay, γ_dephase, V_nn, δ1, δ2, Δ1_0, Δ2_0, Δac_0, Δb_0, T1, T2_y, T2_z, T3, T4, T5 = unpack_params()

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
    δ1::Float64
    δ2::Float64
end

# struct ancilla_parameters
#     Ω::Float64
#     γ_Decay::Float64
#     γ_dephase::Float64
#     V_nn::Float64
#     δ::Float64
# end


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
        Δt = Δ1_0 - p.δ1 * t
        return [p.Ω/2, p.Ω/2, Δt, Δt, p.V_nn, p.V_nn]

    #### Hadamard Mode 1    
    elseif mode == :T2
        return [p.Ω/2, p.Ω/2, p.Ω/2]
    
    #### Ancilla driving mode
    elseif mode == :T3
        Δt = Δ2_0 - p.δ2 * t
        return [p.Ω/2, p.Ω/2, Δt, Δt, p.V_nn/(2^6), p.V_nn/(2^6), p.V_nn/(2^6), p.V_nn/(2^6)]
    
    #### Correction mode for Atom A or C   
    elseif mode == :T4a
        Δt = Δac_0 - p.δ2 * t
        return [p.Ω/2, Δt, p.V_nn, p.V_nn/(2^6)]

    #### Correction mode for Atom B
    elseif mode == :T4b
        Δt = Δb_0 - p.δ2 * t
        return [p.Ω/2, Δt, p.V_nn, p.V_nn, p.V_nn/(2^6), p.V_nn/(2^6)]

    #### Hadamard Mode 2
    elseif mode == :T5
        return [p.Ω/2, p.Ω/2, p.Ω/2]

    end

end

#Lindbald Operators
function lindbaldian_dephase(γ_dephase::Float64,site::Array)
    """
    Function to calculate the Lindbaldian dephase operator for the MCWF method.
    The dephase operators are returned according to given respective sites.

    Args:
        γ_dephase:: Float64: Dephase rate
        site:: Array: Array of sites at which the decay operator acts

    Returns:
        C:: Array{Matrix}: Array of dephase operators acting on the system
    """
    basis = NLevelBasis(2)
    σ_z= transition(basis,1,1) - transition(basis,2,2)
    σ_z = Operator(σ_z.basis_l, σ_z.basis_r, SparseMatrixCSC{ComplexF32, Int64}(σ_z.data))

    identity = transition(basis,1,1) + transition(basis,2,2)
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
basis = NLevelBasis(2)
n = transition(basis,2,2)
n = Operator(n.basis_l, n.basis_r, SparseMatrixCSC{ComplexF32, Int64}(n.data))
σx = transition(basis,1,2) + transition(basis,2,1)
σx = Operator(σx.basis_l, σx.basis_r, SparseMatrixCSC{ComplexF32, Int64}(σx.data))

# Parameters for the Atoms -- coefficients for the Hamiltonian
p = qubit_parameters(Ω,γ_Decay,γ_dephase,V_nn,δ1,δ2)



######################################################################################################## Encoding Atoms - Step 1
#Timespan for driving atoms A-B-C
const tspan1 = [0.0:0.1:T1;]

# System Hamiltonian 1 - driving atoms A-C
σx_a = full_operator(σx, total_qubits, [1])
σx_b = full_operator(σx, total_qubits, [2])
σx_c = full_operator(σx, total_qubits, [3])
const n_a = full_operator(n, total_qubits, [1])
const n_b = full_operator(n, total_qubits, [2])
const n_c = full_operator(n, total_qubits, [3])
const nn_ab = full_operator(n, total_qubits, [1,2])
const nn_bc = full_operator(n, total_qubits, [2,3])
const coeff1 = [t->get_qubit_parameters(p,t,:T1)]
const H1 = LazySum([coeff1[1](tspan1[1])[i] for i ∈ 1:6],[σx_a, σx_c, n_a, n_c, nn_ab, nn_bc])



######################################################################################################## Applying hadamards - Step 2
const tspan2_y = [0.0:0.1:T2_y;]  #Time span for applying hadamards on atoms A-B-C
const tspan2_z = [0.0:0.1:T2_z;]  #Time span for applying hadamards on atoms A-B-C
σy = -im * transition(basis,1,2) + im * transition(basis,2,1)
σy_a = full_operator(σy, total_qubits, [1])
σy_b = full_operator(σy, total_qubits, [2])
σy_c = full_operator(σy, total_qubits, [3])

σz = transition(basis,1,1) - transition(basis,2,2)
σz_a = full_operator(σz, total_qubits, [1])
σz_b = full_operator(σz, total_qubits, [2])
σz_c = full_operator(σz, total_qubits, [3])
const coeff2 = [t->get_qubit_parameters(p,t,:T2)]
const H2_y = LazySum([coeff2[1](tspan2_y[1])[i] for i ∈ 1:3],[σy_a, σy_b, σy_c])
const H2_z = LazySum([coeff2[1](tspan2_z[1])[i] for i ∈ 1:3],[σz_a, σz_b, σz_c])


######################################################################################################## Driving atoms 1-2 - Step 3
const tspan3 = [0.0:0.1:T3;]  #Time span for driving atoms 1-2

σx_1 = full_operator(σx, total_qubits, [4])
σx_2 = full_operator(σx, total_qubits, [5])
const n_1 = full_operator(n, total_qubits, [4])
const n_2 = full_operator(n, total_qubits, [5])
nn_a1 = full_operator(n, total_qubits, [1,4])
nn_b1 = full_operator(n, total_qubits, [2,4])
nn_b2 = full_operator(n, total_qubits, [2,5])
nn_c2 = full_operator(n, total_qubits, [3,5])
const coeff3 = [t->get_qubit_parameters(p,t,:T3)]
const H3 = LazySum([coeff3[1](tspan3[1])[i] for i ∈ 1:8],[σx_1, σx_2, n_1, n_2, nn_a1, nn_b1, nn_b2, nn_c2])


const tspan = [0.0:0.1:(T1+T2_z+T2_y);]

function Ht(t)
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

    if t<T1 || t==T1
        coeffs = coeff1[1](t)
        for i in eachindex(coeffs)
            H1.factors[i] = coeffs[i]
        end
        return H1

    elseif t<(T1+T2_z) || t==(T1+T2_z)
        coeffs = coeff2[1](t)
        for i in eachindex(coeffs)
            H2_z.factors[i] = coeffs[i]
        end
        return H2_z
    
    elseif t<(T1+T2_y+T2_z) || t==(T1+T2_y+T2_z)
        coeffs = coeff2[1](t)  # Update to subtract (T1 + T2_z)
        for i in eachindex(coeffs)
            H2_y.factors[i] = coeffs[i]
        end
        return H2_y

    # elseif t<(T1+T2+T3) || t==(T1+T2+T3)
    #     coeffs = coeff3[1](t-(T1+T2))
    #     for i in eachindex(coeffs)
    #         H3.factors[i] = coeffs[i]
    #     end
        # return H3
    end
end


#Helper function for mcwf_dynamic
const C = lindbaldian_dephase(γ_dephase,[1,2,3,4,5])    
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


######################################################################################################## Error Correction - Step 4
const tspan4 = [0.0:0.1:T4;]  #Time span for error correction of atom A or C or B

#Required Matrix Constants
n_abc = full_operator(n, total_qubits, [1,2,3])
n_abc = Operator(n_abc.basis_l, n_abc.basis_r, SparseMatrixCSC{ComplexF32, Int64}(n_abc.data))

#Hamiltonian for Error Correction of Atom A
const coeff4 = [t->get_qubit_parameters(p,t,:T4a)]
const H_correct_a = LazySum([coeff4[1](tspan4[1])[i] for i ∈ 1:4],[σx_a, n_a, nn_ab, nn_a1])

#Hamiltonian for Error Correction of Atom B
const coeff5 = [t->get_qubit_parameters(p,t,:T4b)]
const H_correct_b = LazySum([coeff5[1](tspan4[1])[i] for i ∈ 1:6],[σx_b, n_b, nn_ab, nn_bc, nn_b1, nn_b2])

#Hamiltonian for Error Correction of Atom C
const H_correct_c = LazySum([coeff4[1](tspan4[1])[i] for i ∈ 1:4],[σx_c, n_c, nn_bc, nn_c2])

#Hamiltonian for No Correction -- Zero Hamiltonian
const H_no_correct = LazySum([0.0],[σx_a])

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



###############################################################################################Timespan for applying hadamards - Step 5
const tspan6 = [0.0:0.1:T5;]  #Time span for driving atoms 1-2

σy = -im * transition(basis,1,2) + im * transition(basis,2,1)
σy_a = full_operator(σy, total_qubits, [1])
σy_b = full_operator(σy, total_qubits, [2])
σy_c = full_operator(σy, total_qubits, [3])
const coeff6 = [t->get_qubit_parameters(p,t,:T5)]
const H_end = LazySum([coeff6[1](tspan6[1])[i] for i ∈ 1:3],[σy_a, σy_b, σy_c])


function Ht_end(t)
"""
Function to calculate the time dependent Hamiltonian for the MCWF method for applying hadamrd back to atoms A-B-C.
    - H_end: Hamiltonian for applying hadamards on atoms A-B-C for time T4:T5

Args:
    t:: Float64: Time
Returns:
    H:: LazySum: Time dependent Hamiltonian
"""
    
    coeffs = coeff6[1](t)
    for i in eachindex(coeffs)
        H_end.factors[i] = coeffs[i]
    end
    return H_end
end

function f_end(t,ψ)
    """
    Function to calculate the time evolution of the system using the MCWF method.
    
    Args:
        t:: Float64: Time
    
    Returns:
        H:: LazySum: Time dependent Hamiltonian
        C:: Array{Operator}: Array of decay operators acting on the system
        Cdagger:: Array{Operator}: Array of adjoint decay operators acting on the system
    """
    
    H = Ht_end(t)
    return H, Ct(t)...
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


# p_tuple = (Ω,γ_Decay,γ_dephase,V_nn,δ,Δ1_0,Δ2_0,OPS)


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
#         Δ_t = Δ1_0 + δ*t

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
#     energy_levels, _, _ = energy_level_spaghetti(t,p_tuple)

#     #ggg energy level
#     e_ggg = [subarray[4] for subarray in energy_levels]

#     #rrr energy level
#     e_rrr = [subarray[1] for subarray in energy_levels]

#     Δ0_dyn, _ = trapezoidal_integrate(t,e_ggg,e_rrr)
#     Δ0_dyn = mod(-Δ0_dyn,2π)

#     return Δ0_dyn
# end

# function apply_dynamical_phase(Δ0_dyn::Float64)
#     """
#     Function to apply the dynamical phase correction to the system.
#     Args:
#         - Δ0_dyn:: Float64: Dynamical phase

#     Returns:
#         - ψ_target:: Array: Target state after applying the dynamical phase correction
#     """

#     ψ_target = α * reduce(kron,[g,g,g,g,g]) + exp(im * Δ0_dyn) * β * reduce(kron,[r,r,r,g,g])
#     ψ_target = ψ_target / norm(ψ_target)
#     return ψ_target
# end
