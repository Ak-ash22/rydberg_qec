include("systemparams.jl")

const σ_x, n, Π_g, n, I, σ_minus, σ_plus, σ_z = system_constants()
const Ω1, Ω2, γ_Decay, γ_dephase, V1_nn, V2_nn, δ, Δ1_0, Δ2_0, T_optimal = unpack_params()


#Defining the struct for the parameters
struct Parameters
    Ω1::Float64
    Ω2::Float64
    γ_Decay::Float64
    γ_dephase::Float64
    V1_nn::Float64
    V2_nn::Float64
    δ::Float64
end

function (p::Parameters)(t)
    """
    Function to get the parameters of the system at time t for both the atoms

    Parameters: 
        - time t: Float64

    Returns: 
        - parameters at time t: Tuple
    """
    tf = 0
    if t<tf
        return (Δ1_0-p.δ *t, Δ2_0, p.Ω1, p.Ω2, p.γ_Decay, p.γ_dephase, p.V1_nn, p.V2_nn)
    else
        return (Δ1_0-p.δ *t, Δ2_0-p.δ *(t-tf), p.Ω1, p.Ω2, p.γ_Decay, p.γ_dephase, p.V1_nn, p.V2_nn)
    end
end


function master_eqn(dρ,ρ,p,t)
    """ 
    Function to solve the master equation for the system of 1+2 atoms and update the time derivative of density matrix at time t

    Parameters:
        - dρ: Array, ρ: Array, p: Parameters(function), t: Float64
    Returns: 
        --
    """
    #parameters
    Δ1_t, Δ2_t, Ω1, Ω2, γ_Decay, γ_dephase, V1_nn, V2_nn = p(t)

    #Hamiltonian elements
    nn_sys = [kron(n,n,I), kron(I,n,n), kron(n,I,n)]
    σminus_sys = [reduce(kron,[I,σ_minus,I]),reduce(kron,[I,I,σ_minus])]
    σplus_sys = [reduce(kron,[I,σ_plus,I]), reduce(kron,[I,I,σ_plus])]
    σz_sys = [reduce(kron,[I,σ_z,I]), reduce(kron,[I,I,σ_z])]
    n_sys = [reduce(kron,[I,n,I]), reduce(kron,[I,I,n])]
    σx_sys = [reduce(kron,[I,σ_x,I]), reduce(kron,[I,I,σ_x])]

    #Hamiltonian
    H = Ω1/2 .* σx_sys[1] + Ω2/2 .* σx_sys[2] + V1_nn .* nn_sys[1] + V2_nn .* nn_sys[3] +
        Δ1_t .* n_sys[1] + Δ2_t .*n_sys[2]

    #lindbaldian terms
    l_decay = zeros(2^(n_atoms+1),2^(n_atoms+1))
    l_dephase = zeros(2^(n_atoms+1),2^(n_atoms+1))

    for i in 1:2
        l_decay += γ_Decay.*(σminus_sys[i] * ρ * σplus_sys[i] - 0.5.*((σplus_sys[i]*σminus_sys[i]*ρ)
                         + (ρ*σplus_sys[i]*σminus_sys[i])))
        l_dephase += γ_dephase.*(σz_sys[i]* ρ *σz_sys[i] - ρ)
    end

    dρ .= -1im .* (H*ρ - ρ*H) + l_decay + l_dephase
end    


function rydberg_populations(sol)
    """
    Function to calculate the populations of the Rydberg states -
        rydberg1: Population of the first atom in Rydberg state
        rydberg2: Population of the second atom in Rydberg state
        rydberg12: Population of the both atoms in Rydberg states

    Parameters: 
        - sol: ODESolution object

    Returns:
        - lists containing real values of the populations
    """
    a = reduce(kron,[I,n,I])
    rydberg1 = [real(tr(a*ρ)) for ρ in sol.u]

    b = reduce(kron,[I,I,n])
    rydberg2 = [real(tr(b*ρ)) for ρ in sol.u]

    c = reduce(kron,[I,n,n])
    rydberg12 = [real(tr(c*ρ)) for ρ in sol.u]

    return rydberg1, rydberg2, rydberg12
end


function solve_master_eqn(p::Parameters, tspan)
    """
    Function to solve the master equation using a numerical ODE solver to get the density matrix at the final time

    Parameters: 
        -p: Parameters(tuple), tspan: Tuple

    Returns: 
        -ODESolution object
    """
    ρ_0 = initialize_system()
    eqn = ODEProblem(master_eqn, ρ_0, tspan, p)
    sol = solve(eqn, Rodas3(autodiff=false), saveat = 1)
    return sol
end