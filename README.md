
<p align="center">
  
</p>

# <img src="Ulysses-logo.png" alt="ULYSSES logo" width="90"> ULYSSES: Universal LeptogeneSiS Equation Solver

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![ReadTheDocs](https://readthedocs.org/projects/ulysses-universal-leptogenesis-equation-solver/badge/?version=latest)](https://ulysses-universal-leptogenesis-equation-solver.readthedocs.io/en/latest/)



## Introduction

**ULYSSES** is an open-source Python package designed to compute the baryon asymmetry of the Universe through **leptogenesis**, within the **type-I seesaw** framework. It numerically solves the momentum-averaged **Boltzmann equations (BEs)** or **quantum density matrix equations (DMEs)** that govern the evolution of particle asymmetries in the early Universe.

The package features:
- Fast and modular solver infrastructure
- Flexible physics model input via plugin system
- Full compatibility with user-defined scenarios

It is documented across two key publications:
- [Version 1 – arXiv:2007.09150](https://arxiv.org/abs/2007.09150)
- [Version 2 – arXiv:2301.05722](https://arxiv.org/abs/2301.05722)

---

## Scientific Background

The baryon asymmetry is measured as the baryon-to-photon ratio:  

$$\eta_B \approx 6 \times 10^{-10}$$

**Leptogenesis** explains this asymmetry via CP-violating decays of heavy right-handed neutrinos (RHNs) in the early Universe. These decays produce a net lepton number, which electroweak sphalerons convert to baryon number.

ULYSSES implements this through coupled differential equations:

$$\frac{dN_N}{dz} = -D (N_N - N_N^{\text{eq}}), \quad \frac{dN_{B-L}}{dz} = \varepsilon D (N_N - N_N^{\text{eq}}) - W N_{B-L}$$

with temperature parameter $z = M_1 / T$.

---

## Core Features

### Version 1 (v1) – arXiv:2007.09150
- Momentum-averaged BEs
- Flavored/unflavored solvers
- Resonant leptogenesis
- Optional scatterings and spectator processes
- Plugin-based model extensions

### Version 2 (v2) – arXiv:2301.05722
- **ARS** (low-scale) leptogenesis with 2 RHNs (`BEARS`, `BEARS_INTERP`)
- **Primordial black hole** (PBH) leptogenesis
- **Complete BEs** with full quantum statistics (Cases 2–4)
- Preconfigured 2D scans for parameter space exploration

### Version 3 (v3) – This release
- **Case S2**: Full phase-space Boltzmann equations including $\Delta L = 1$ scattering ($s$- and $t$-channel, $Nl \to qt$) via precomputed AB grids (arXiv:0907.0205)
- **ARS with 3 RHNs** (`BEARS_3RHN`): extension of the ARS mechanism to the three right-handed neutrino case
- **Extended model interface** (`--extended` flag): allows models to declare extra parameters beyond the standard PMNS+CI set, opening the solver to coupled multi-sector scenarios; `1BE1F_DM_FreezeIn` serve as worked toy-model examples demonstrating freeze-in and RHN dark matter production
---

## Installation

### PyPI (Recommended)
```bash
pip install ulysses --user
```

### From Source (for developers)

If you intend to modify the source code or contribute to the project, install from a local clone of the GitHub repository:

1. Clone the repository:
```bash
git clone https://github.com/earlyuniverse/ulysses.git
```
2. Navigate into the directory:
```bash
cd ulysses
```
3. Install in editable mode:
```bash
pip install -e . --user
```

### macOS notes

ULYSSES builds a small C library (`NumbaQuadpack`) at install time using CMake. On macOS you need:

1. **Xcode Command Line Tools** (provides the C compiler and linker):
   ```bash
   xcode-select --install
   ```
2. **CMake** — installed automatically by pip from `pyproject.toml`, but you can also install it via Homebrew:
   ```bash
   brew install cmake
   ```

Apple Silicon (M1/M2/M3) is fully supported. If you see an `ImportError: Could not find libcquadpack.dylib` after a successful build, it usually means an older version of ULYSSES was installed that compiled the library with a wrong suffix — a clean reinstall fixes it:
```bash
pip uninstall ulysses -y
pip install .
```

### Environment Setup (optional)

Add to your `~/.bashrc` or `~/.zshrc`:

```bash
export ULYSSES=/path/to/ulysses
export PYTHONPATH=$PYTHONPATH:$ULYSSES
export PATH=$PATH:$ULYSSES/bin
```

---

## Quick Start

### Step 1: Create a parameter file

```ini
# Lightest neutrino mass (log10, eV)
m      -2
# RHN masses (log10, GeV)
M1     11
M2     11.6
M3     12
# Casas-Ibarra angles (deg)
delta  213
a21     81
a31    476
x1      90
x2      87
x3     180
y1    -120
y2       0
y3    -120
# PMNS angles — NuFit 6.1 best fit, NO
t12    33.76
t13     8.62
t23    43.27
```

### Step 2: Run a simulation

```bash
uls-calc -m 1BE1F examples/1N1F.dat -o evolution.pdf
```

The terminal will display the ULYSSES banner, the model name, all input parameters, and the computed $\eta_B$, $Y_B$, and $\Omega_b h^2$.

---

## Advanced Usage

ULYSSES includes a suite of command-line tools in `bin/`:

| Command       | Description                                                                 |
|---------------|-----------------------------------------------------------------------------|
| `uls-calc`    | Single-point evolution of the asymmetry from a parameter card.             |
| `uls-scan`    | 1D scan over a user-defined parameter range.                                |
| `uls-scan2D`  | 2D grid scan over two parameters.                                           |
| `uls-nest`    | Nested sampling scan for Bayesian inference or model selection.             |
| `uls-models`  | List all available pre-defined physics models.                              |

### Extended mode (model-specific parameters)

Models that require parameters beyond the standard PMNS+CI set (e.g. dark matter modules) use the `--extended` flag:

```bash
uls-calc -m 1BE1F_DM_FreezeIn --extended examples/1N1F_dm.dat -o dm_evolution.pdf
```

---

## Available Physics Models

### Standard Boltzmann equation models

| Model      | Example file  | Description                                      |
|------------|---------------|--------------------------------------------------|
| `1BE1F`    | `1N1F.dat`    | One-flavour BE, 1 RHN                            |
| `1BE2F`    | `1N2F.dat`    | Two-flavour BE, 1 RHN                            |
| `1BE3F`    | `1N3F.dat`    | Three-flavour BE, 1 RHN                          |
| `2BE1F`    | `2N1F.dat`    | One-flavour BE, 2 RHNs                           |
| `2BE2F`    | `2N2F.dat`    | Two-flavour BE, 2 RHNs                           |
| `2BE3F`    | `2N3F.dat`    | Three-flavour BE, 2 RHNs                         |
| `1BE1Fsf`  | `1N1F.dat`    | `1BE1F` evolved in scale factor                  |

### Density matrix equation models

| Model      | Example file  | Description                                      |
|------------|---------------|--------------------------------------------------|
| `1DME`     | `1N3F.dat`    | Density matrix equations, 1 RHN                  |
| `2DME`     | `2N3F.dat`    | Density matrix equations, 2 RHNs                 |
| `3DME`     | `3N3F.dat`    | Density matrix equations, 3 RHNs                 |
| `3DMEsct`  | `3N3F.dat`    | 3DME with scattering effects                     |

### Resonant leptogenesis

| Model      | Example file  | Description                                      |
|------------|---------------|--------------------------------------------------|
| `2RES`     | `Res.dat`     | 2BE3F in the resonant regime (experimental)      |
| `2RESmix`  | `Res.dat`     | Resonant leptogenesis with flavour mixing        |
| `2RESsp`   | `Res.dat`     | Resonant leptogenesis with spectator processes   |

### Complete Boltzmann equations (full quantum statistics)

| Model            | Example file | Description                                                |
|------------------|--------------|------------------------------------------------------------|
| `1BE1F_Case2`    | `1N1F.dat`   | Full quantum stats, kinetic equilibrium assumed            |
| `1BE1F_Case3`    | `1N1F.dat`   | Full quantum stats, no kinetic equilibrium                 |
| `1BE1F_Case4`    | `1N1F.dat`   | Full quantum stats, no kinetic equilibrium, FD/BE          |
| `1BE1F_CaseS2`   | `1N1F.dat`   | Case 4 + $\Delta L=1$ scattering via precomputed AB grids  |

### ARS (low-scale) leptogenesis

| Model           | Example file          | Description                                   |
|-----------------|-----------------------|-----------------------------------------------|
| `BEARS`         | `2RHNosc.dat`         | Temperature-independent ARS, 2 RHNs           |
| `BEARS_INTERP`  | `2RHNosc.dat`         | Temperature-dependent ARS, 2 RHNs             |
| `BEARS_3RHN`    | `Bears_3RHN_alt.dat`  | ARS with 3 RHNs                               |

### Primordial black holes

| Model        | Example file | Description                                      |
|--------------|--------------|--------------------------------------------------|
| `1BE1F_PBH`  | `PBH.dat`    | Leptogenesis from PBH evaporation                |

### Dark matter (Toy model for extended flag interface)

| Model                  | Example file   | Description                                              |
|------------------------|----------------|----------------------------------------------------------|
| `1BE1F_DM_FreezeIn`    | `1N1F_dm.dat`  | Freeze-in DM from RHN decays via feeble dark coupling    |

---

## Citation

Please cite the relevant paper(s) for the features you use:

### ULYSSES v1
```bibtex
@article{Granelli:2020pim,
    author = "Granelli, Alessandro and Moffat, Kristian and Perez-Gonzalez, Yuber F. and Schulz, Holger and Turner, Jessica",
    title = "{ULYSSES: Universal LeptogeneSiS Equation Solver}",
    eprint = "2007.09150",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "FERMILAB-PUB-20-275-T, SISSA 17/2020/FISI, IPPP/20/30",
    doi = "10.1016/j.cpc.2020.107813",
    journal = "Comput. Phys. Commun.",
    volume = "262",
    pages = "107813",
    year = "2021"
}
```

### ULYSSES v2
```bibtex
@article{Granelli:2023vcm,
    author = "Granelli, Alessandro and Leslie, Christopher and Perez-Gonzalez, Yuber F. and Schulz, Holger and Shuve, Brian and Turner, Jessica and Walker, Rosie",
    title = "{ULYSSES, universal LeptogeneSiS equation solver: Version 2}",
    eprint = "2301.05722",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "IPPP/23/02",
    doi = "10.1016/j.cpc.2023.108834",
    journal = "Comput. Phys. Commun.",
    volume = "291",
    pages = "108834",
    year = "2023"
}
```

---

## Documentation, License, and Contributing

### Documentation

- **Version 1 Manual**: [arXiv:2007.09150](https://arxiv.org/abs/2007.09150)
- **Version 2 Manual**: [arXiv:2301.05722](https://arxiv.org/abs/2301.05722)
- **API Documentation**: [Read the Docs](https://ulysses-universal-leptogenesis-equation-solver.readthedocs.io/)

### License

ULYSSES is licensed under the **MIT License**. See the [LICENSE file](https://github.com/earlyuniverse/ulysses/blob/main/LICENSE) for full details.

### Contributing

We welcome contributions from the community:

- Report bugs or feature requests via [GitHub Issues](https://github.com/earlyuniverse/ulysses/issues)
- Submit pull requests for code fixes, model implementations, or documentation
- Share new physics modules via the plugin architecture

For major changes, please open a discussion before proceeding.

---

## Useful Links

- **GitHub Repository**: [earlyuniverse/ulysses](https://github.com/earlyuniverse/ulysses)
- **arXiv:2007.09150** – ULYSSES v1: [link](https://arxiv.org/abs/2007.09150)
- **arXiv:2301.05722** – ULYSSES v2: [link](https://arxiv.org/abs/2301.05722)
- **Documentation (ReadTheDocs)**: [link](https://ulysses-universal-leptogenesis-equation-solver.readthedocs.io/)

---
