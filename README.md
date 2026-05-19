
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
- [Version 3 – arXiv:2605.16540](https://arxiv.org/abs/2605.16540)

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

### Version 3 (v3) – arXiv:2605.16540
- **ARS with 3 RHNs** (`BEARS_3RHN`): extension of the ARS mechanism to the three right-handed neutrino case, with optional mass-dependent indirect rate corrections (`--ars-indirect`)
- **Case S2**: Full phase-space Boltzmann equations including $\Delta L = 1$ scattering ($s$- and $t$-channel, $Nl \to qt$) without assuming kinetic equilibrium nor Maxwell-Boltzmann statistics
- **Extended model interface** (`--extended` flag): allows models to declare extra parameters beyond the standard PMNS+CI set, opening the solver to coupled multi-sector scenarios; `1BE1F_DM_FreezeIn` serves as a worked toy-model example demonstrating freeze-in and RHN dark matter production
- **Three Yukawa parameterisations**: Casas–Ibarra with Euler angles (default), single-imaginary CI (arXiv:2106.16226), and direct matrix entry — auto-detected from the parameter file
---

## Installation

### PyPI (Recommended)

```bash
pip install ulysses
```

Pre-built binary wheels are provided for Linux (x86\_64 and aarch64), macOS (Intel and Apple Silicon), and Windows, covering Python 3.9–3.13. **No C compiler is required** when a wheel is available for your platform.

A C compiler is only needed if pip falls back to building from source (i.e. on an unsupported platform or when installing directly from the repository). On Linux/macOS this means `gcc` or `clang`; on Windows, MSVC (Visual Studio Build Tools).

### From Source (for developers)

If you intend to modify the source code or contribute to the project, install from a local clone:

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
pip install -e .
```

This will compile the `NumbaQuadpack` C extension using your system compiler (`gcc`/`clang` on Linux/macOS, MSVC on Windows).

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

The parameter file is a plain-text card with one `key value` pair per line. Lines starting with `#` are comments.

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
# PMNS angles (deg) — optional, default to NuFit 6.1 best fit (NO)
t12    33.76
t13     8.62
t23    43.27
```

`t12`, `t13`, `t23` are optional. If omitted, ULYSSES uses the NuFit 6.1 best-fit values for the chosen mass ordering (normal ordering by default).

### Step 2: Run a simulation

```bash
uls-calc -m 1BE1F examples/1N1F.dat -o evolution.pdf
```

The terminal displays the ULYSSES banner, model name, all input parameters, and the computed $\eta_B$, $Y_B$, and $\Omega_b h^2$. The output file can be `.pdf`, `.txt`, `.dat`, or `.csv`.

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

### CLI flags

All tools share a common set of flags. The most commonly used ones for `uls-calc` are:

| Flag | Default | Description |
|------|---------|-------------|
| `-m MODEL` | `1DME` | Physics model to run (see model table below) |
| `-o FILE` | — | Output file for evolution plots or data (`.pdf`, `.txt`, `.dat`, `.csv`) |
| `--inv` | off | Use inverted neutrino mass ordering |
| `--loop` | off | Use loop-corrected Yukawa couplings |
| `--initial ` | `0` | Initial RHN abundance: `0` = vanishing, `1` = thermal |
| `--extended` | off | Allow model-specific parameters beyond the standard PMNS+CI set |
| `--lambda ` | `1e3` | ARS quasi-static cutoff scale $\Lambda$ (GeV) for `BEARS_3RHN` |
| `--ars-indirect` | off | Enable mass-dependent hindrance rate corrections in `BEARS_3RHN` |
| `--zrange Z0,Z1,N` | `0.1,30,500` | Start, end, and number of steps for the $z$ evolution variable |
| `--xrange X0,X1,N` | `1e-6,min([1, 20*131.7/M1]),500` | Start, end, and steps for the 3RHN ARS $x = T_\text{ew}/T$ variable |
| `--zcut` | `1.0` | Stitch cut value used in 2RHN ARS models |
| `-v` | off | Enable debug output |

### Extended mode (model-specific parameters)

Models that require parameters beyond the standard PMNS+CI set (e.g. dark matter modules) use the `--extended` flag:

```bash
uls-calc -m 1BE1F_DM_FreezeIn --extended examples/1N1F_dm.dat -o dm_evolution.pdf
```

The `--extended` flag tells the solver to pass unrecognised keys in the parameter file through to the model rather than rejecting them as errors. For `1BE1F_DM_FreezeIn`, the two extra keys are:

```ini
# --- standard PMNS+CI block (same as 1BE1F) ---
m      -100
M1     12
M2     15
M3     16
x1    180
y1      1.4
x2    180
y2     11.2
x3    180
y3     11
delta  270
a21      0
a31      0

# --- model-specific extensions (require --extended) ---
lam    1e-7    # feeble dark coupling λ
m_dm   1e6    # dark matter mass (GeV)
```

---

## Yukawa Parameterisations

ULYSSES supports three ways to specify the Yukawa coupling matrix. The parameterisation is **auto-detected** from the keys present in the parameter file — no explicit flag is needed.

### 1. Casas–Ibarra with Euler angles (default)

Standard parameterisation. Use when the parameter file contains `x1`, `y1`, `x2`, …

```ini
x1   90
y1  -120
x2   87
y2    0
x3  180
y3  -120
```

### 2. Single-imaginary Casas–Ibarra (arXiv:2106.16226)

Activated when the parameter file contains `xN1`. Used for ARS models with 3 RHNs.

```ini
xN1   2
xN2   3
xnu1  4
xnu2  5
x     1
y   103.13
```

### 3. Manual Yukawa matrix

Activated when the parameter file contains `Y11_mag`. Each entry is specified as a magnitude and a phase (radians).

```ini
Y11_mag  8.149e-05
Y11_phs -1.819
Y12_mag  8.149e-05
Y12_phs -0.248
# ... (all 9 entries required)
```

---

## Neutrino Mass Splittings

The default mass squared splittings come from **NuFit 6.1** (2025):

| Parameter | Default value | Description |
|-----------|---------------|-------------|
| `m2solar` | $7.537 \times 10^{-5}$ eV² | Solar splitting $\Delta m^2_{21}$ |
| `m2atm`   | $2.521 \times 10^{-3}$ eV² | Atmospheric splitting, normal ordering |
| `m2atminv`| $2.500 \times 10^{-3}$ eV² | Atmospheric splitting, inverted ordering |

To override, add the key(s) to your parameter file (values in eV²):

```ini
m2solar   7.53e-5
m2atm     2.52e-3
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

| Model           | Example file          | Description                                                      |
|-----------------|-----------------------|------------------------------------------------------------------|
| `BEARS`         | `2RHNosc.dat`         | Temperature-independent ARS, 2 RHNs                             |
| `BEARS_INTERP`  | `2RHNosc.dat`         | Temperature-dependent ARS, 2 RHNs                               |
| `BEARS_3RHN`    | `Bears_3RHN_alt.dat`  | ARS with 3 RHNs; supports `--ars-indirect` and `--lambda` flags  |

### Primordial black holes

| Model        | Example file | Description                                      |
|--------------|--------------|--------------------------------------------------|
| `1BE1F_PBH`  | `PBH.dat`    | Leptogenesis from PBH evaporation                |

### Dark matter (toy model for extended flag interface)

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

### ULYSSES v3
Paper forthcoming. In the meantime, please cite both v1 and v2 above when using v3 features.

---

## Documentation, License, and Contributing

### Documentation

- **Version 1 Manual**: [arXiv:2007.09150](https://arxiv.org/abs/2007.09150)
- **Version 2 Manual**: [arXiv:2301.05722](https://arxiv.org/abs/2301.05722)
- **Version 3 Manual**: [arXiv:2605.16540](https://arxiv.org/abs/2605.16540)
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
- **arXiv:2605.16540** – ULYSSES v3: [link](https://arxiv.org/abs/2605.16540)
- **Documentation (ReadTheDocs)**: [link](https://ulysses-universal-leptogenesis-equation-solver.readthedocs.io/)

---
