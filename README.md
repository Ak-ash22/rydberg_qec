# Quantum Error Correction Codes using Rydberg Facilitation

[![Julia](https://img.shields.io/badge/Julia-1.11.4-9558B2?logo=julia)](https://julialang.org/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Overview

This repository contains the complete implementation and analysis of quantum error correction (QEC) codes using Rydberg atom facilitation mechanisms. The work was completed as part of a Master's thesis in the Erasmus Mundus Quantum Technology and Engineering (QuanTEEM) program at RPTU, Kaiserslautern, Germany.

The project explores primitive quantum error correcting codes through the implementation on Rydberg atom systems, leveraging the unique properties of Rydberg interactions for error detection and correction.

## Key Features

- **Multiple QEC Implementations**: Bit-flip, phase-flip, and Shor's 9-qubit code
- **Rydberg Atom Simulations**: Full quantum dynamics simulation using master equation solvers
- **Monte Carlo Wave Function Method**: Stochastic trajectory simulations for open quantum systems
- **Error Analysis**: Comprehensive error detection and correction protocols with storage
- **Scalability Studies**: Analysis from single atoms to multi-qubit systems
- **Parallel Computing Support**: High-performance computing capabilities for large-scale simulations

## Table of Contents

- [Installation](#installation)
- [Project Structure](#project-structure)
- [Quantum Error Correction Codes](#quantum-error-correction-codes)
- [Usage](#usage)
- [System Parameters](#system-parameters)
- [Results and Analysis](#results-and-analysis)
- [Dependencies](#dependencies)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Installation

### Prerequisites

- Julia 1.11.4 or higher
- Recommended: 16GB+ RAM for large-scale simulations
- Optional: HPC cluster access for parallel computations

### Setup

1. Clone the repository:
```bash
git clone https://github.com/Ak-ash22/rydberg_qec.git
cd rydberg_qec
```

2. Start Julia and activate the project environment:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

3. Verify installation:
```julia
using QuantumOptics
using Plots
```

## Project Structure

```
rydberg_qec/
├── codework/
│   ├── 1atom/              # Single atom dynamics and adiabatic evolution
│   ├── 2atoms/             # Two-atom system studies
│   ├── analytical_work/    # Analytical calculations and parameter optimization
│   ├── bitflip_code/       # 3-qubit bit-flip error correction code
│   │   ├── error_correction_with_storage.jl
│   │   ├── error_detection.jl
│   │   ├── functions.jl
│   │   ├── system_params.jl
│   │   └── cluster_data_analysis/
│   ├── phaseflip_code/     # 3-qubit phase-flip error correction code
│   │   ├── error_detection_and_correction.jl
│   │   ├── functions.jl
│   │   ├── system_params.jl
│   │   └── cluster_data_analysis/
│   └── shor_code/          # 9-qubit Shor code implementation
│       ├── error_correction.jl
│       ├── functions.jl
│       ├── system_params.jl
│       ├── parallel_main.jl
│       └── cluster_data_analysis/
├── images/                 # Generated plots and figures
├── Project.toml           # Julia project dependencies
└── README.md
```

## Quantum Error Correction Codes

### 1. Bit-Flip Code (3-Qubit)

The bit-flip code protects against bit-flip errors (X errors) by encoding a logical qubit into three physical qubits:
- **Logical States**: |0⟩_L = |000⟩, |1⟩_L = |111⟩
- **Encoding**: Uses Rydberg facilitation for CNOT gate implementation
- **Error Detection**: Two ancilla qubits measure error syndromes
- **Error Correction**: Majority voting with conditional corrections

**Key Files**:
- `codework/bitflip_code/error_correction_with_storage.jl`: Full protocol with storage time
- `codework/bitflip_code/functions.jl`: Hamiltonian construction and operators

### 2. Phase-Flip Code (3-Qubit)

The phase-flip code protects against phase-flip errors (Z errors):
- **Logical States**: |+⟩_L = |+++⟩, |-⟩_L = |---⟩
- **Encoding**: Hadamard gates + bit-flip encoding in dual basis
- **Error Detection**: Three-level atom system for enhanced detection
- **Dynamical Phase Correction**: Handles accumulated dynamical phases

**Key Files**:
- `codework/phaseflip_code/error_detection_and_correction.jl`: Complete protocol
- `codework/phaseflip_code/functions.jl`: Advanced three-level Hamiltonians

### 3. Shor's Code (9-Qubit)

Shor's code provides protection against both bit-flip and phase-flip errors:
- **Encoding**: Combines bit-flip and phase-flip encoding
- **Error Detection**: Eight ancilla qubits for syndrome measurement
- **Universal Protection**: Corrects arbitrary single-qubit errors
- **Scalability**: Demonstrates feasibility of larger QEC codes

**Key Files**:
- `codework/shor_code/error_correction.jl`: Error correction protocol
- `codework/shor_code/parallel_main.jl`: Parallelized simulations

## Usage

### Running Bit-Flip Error Correction

```julia
# Navigate to bitflip code directory
cd("codework/bitflip_code")

# Include the main simulation file
include("error_correction_with_storage.jl")

# Run simulation with 1000 trajectories and storage parameter s=1
main(1000, 1)
```

### Running Phase-Flip Error Correction

```julia
# Navigate to phaseflip code directory
cd("codework/phaseflip_code")

# Include the main simulation file
include("error_detection_and_correction.jl")

# Run simulation with specified parameters
main(1000, 1)  # 1000 trajectories, storage parameter s=1
```

### Running Shor's Code Simulation

```julia
# Navigate to Shor code directory
cd("codework/shor_code")

# Include the parallel implementation
include("parallel_main.jl")

# Run with parallel workers
# This automatically distributes computation across available cores
```
It is advised to run Shor's code simulation only in the HPC cluster together with optimizing few parameters, as the runtime could be really long.

### Analyzing Results

All implementations save results as `.jld2` files. To analyze:

```julia
using JLD2, FileIO

# Load saved data
data = load("results_data/simulation_results.jld2")

# Access fidelity, population, and error statistics
fidelity = data["fidelity"]
populations = data["populations"]
```

Jupyter notebooks in each subdirectory provide detailed analysis and visualization templates.

## System Parameters

### Physical Parameters

The simulations use realistic Rydberg atom parameters:

- **Rabi Frequency (Ω)**: 1.0 MHz (normalized units)
- **Rydberg Interaction (V_nn)**: -1000 MHz (blockade regime)
- **Detuning (Δ)**: 1000-2000 MHz (various protocols)
- **Decay Rate (γ_decay)**: 0.0 - 0.001 MHz (various noise levels)
- **Dephasing Rate (γ_dephase)**: 0.0 - 1.0e-5 MHz

### Customization

Modify parameters in `system_params.jl` files:

```julia
function params()
    return Dict(
        :Ω => 1.0,              # Rabi frequency
        :γ_Decay => 0.001,      # Decay rate
        :γ_dephase => 0.0,      # Dephasing rate
        :V_nn => -1000.0,       # Rydberg interaction
        :Δ_0 => 1000.0,         # Detuning
        # ... additional parameters
    )
end
```

## Results and Analysis

### Performance Metrics

The implementations track multiple performance indicators:

1. **Fidelity**: Overlap with ideal error-corrected state
2. **Population Dynamics**: Time evolution of qubit states
3. **Error Detection Rate**: Probability of correct syndrome measurement
4. **Error Correction Success**: Post-correction fidelity improvement
5. **Storage Time Dependence**: Coherence loss during idle periods

### Data Analysis Tools

- **Jupyter Notebooks**: Interactive analysis in each code directory
  - `result_analysis.ipynb`: Standard analysis workflow
  - `storage_result_analysis.ipynb`: Storage time studies
  - `cluster_result.ipynb`: HPC simulation analysis

- **Cluster Computing**: SLURM scripts for large-scale parameter sweeps
  - `slurm_script.sh`: Job submission templates
  - Parallel trajectory averaging for statistical accuracy

### Visualization

Results are visualized using Plots.jl:
- Fidelity vs. time
- Population dynamics (individual qubits and logical state)
- Error syndrome statistics
- Noise resilience curves

## Dependencies

Core Julia packages (specified in `Project.toml`):

```toml
QuantumOptics = ">=0.7"          # Quantum mechanics framework
Plots = ">=1.20"                  # Visualization
DifferentialEquations = ">=7.0"   # ODE/SDE solvers
JLD2 = ">=0.4"                    # Data storage
SymPy = ">=3.0"                   # Symbolic mathematics
Latexify = ">=0.16"               # LaTeX export
SparseArrays = ">=1.0"            # Efficient sparse matrices
Random = ">=1.0"                  # Random number generation
Distributed = ">=1.0"             # Parallel computing
FileIO = ">=1.0"                  # File I/O operations
```

### Installing Dependencies

All dependencies are automatically installed via:
```julia
Pkg.instantiate()
```

## Theoretical Background

### Rydberg Facilitation Mechanism

Rydberg atoms exhibit strong dipole-dipole interactions, creating a "blockade" effect where nearby atoms cannot simultaneously occupy the Rydberg state. This natural blockade enables:

1. **Two-Qubit Gates**: Controlled interactions for CNOT gates
2. **Error Detection**: Selective coupling to ancilla qubits
3. **State-Selective Operations**: Conditional dynamics based on neighbor states

### Error Correction Protocol

The general protocol follows these steps:

1. **Initialization**: Prepare logical qubit state
2. **Encoding**: Map logical state to physical qubits using facilitation
3. **Storage**: Allow evolution under noise (optional)
4. **Syndrome Measurement**: Couple ancillas to detect errors
5. **Error Correction**: Apply conditional corrections based on syndromes
6. **Decoding**: Map back to logical state (if needed)
7. **Fidelity Measurement**: Compare with ideal target state

### Monte Carlo Wave Function Method

Open quantum system dynamics are simulated using the MCWF approach:
- Stochastic quantum trajectories
- Jump operators for dissipation (decay, dephasing)
- Ensemble averaging over many trajectories
- Faithful representation of measurement backaction

## Citation

If you use this code in your research, please cite:

```bibtex
@mastersthesis{rydberg_qec_2024,
  author = {Akash},
  title = {Quantum Error Correction Codes using Rydberg Facilitation},
  school = {RPTU Kaiserslautern-Landau},
  year = {2024},
  type = {Master's Thesis},
  program = {Erasmus Mundus QuanTEEM}
}
```

## Future Work

Potential extensions of this work:

- [ ] Implementation of surface codes with Rydberg arrays
- [ ] Fault-tolerant gate operations with error correction
- [ ] Realistic noise models from experimental data
- [ ] Integration with quantum algorithms (VQE, QAOA)
- [ ] Hardware-efficient encoding schemes
- [ ] Real-time adaptive error correction

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- **RPTU Kaiserslautern-Landau**: Host institution
- **Erasmus Mundus QuanTEEM Program**: Funding and support
- **QuantumOptics.jl Contributors**: Excellent quantum simulation framework
- **Quantum Error Correction Community**: Theoretical foundations

## Contact

**Author**: Akash Malemath  
**Email**: 14akash2000@gmail.com  
**GitHub**: [@Ak-ash22](https://github.com/Ak-ash22)  
**Repository**: [rydberg_qec](https://github.com/Ak-ash22/rydberg_qec)

---

*For questions, issues, or collaboration opportunities, please open an issue on GitHub or contact the author directly.* 
