# Ω⁴-UNIFIED v80.0

**Topological Origin of Coupling Constants, FRG ↔ RGE Invariant, and SUSY Scale from Calabi–Yau Compactification**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23304779.svg)](https://doi.org/10.5281/zenodo.23304779)
[![Review DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23305721.svg)](https://doi.org/10.5281/zenodo.23305721)
[![License: UAL-v1.1](https://img.shields.io/badge/License-UAL--v1.1-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)

---

## What is this?

Ω⁴-UNIFIED v80.0 is a phenomenological model that derives the running of coupling constants and the SUSY scale from the topology of a three-dimensional Calabi–Yau manifold (Tian–Yau / Z₃).

The model rests on nine structural results:

1. **Topological origin of B_i**: B_i = (32.833, −16.833, −21.000), deviation from RGE < 2%.
2. **FRG ↔ RGE invariant**: g* + θ₁ = 4, Σ(g* + θ₁) = 28 = dim(SO(8)) = χ/22.
3. **Spread formula**: spread = 28 − 8 + ln(28)/256 = 20.013016 (0.000005% from 1-loop RGE).
4. **SUSY scale**: M_SUSY = (2a−1)·M_EW = 3.6933 TeV for a = N_eff = 8.
5. **Membrane wave function**: E_n = {0, 13, 24, 33, 40, 45, 48, 49} for a = 7.
6. **Strong CP problem**: θ_QCD = e^(−c₂) = 1.37×10⁻⁴⁴ — 34 orders below the bound.
7. **α_s(M_Z)** = 0.118064 (0.031% from PDG).
8. **Unified c₂ = 101 → 4 effects**: θ_QCD, M_SUSY, Ω_DM, η_B.
9. **Nuclear physics**: b-quark ↔ Pb-208: 26·8 = 4·52 = 208; five magic numbers.

---

## Key numbers

| Quantity | Value |
|---|---|
| Constants in REGISTRY | 783 |
| Constants in genome | 457 |
| Base model checks | 164/164 |
| Tetrahedron checks | 18/18 |
| Research files | 9/9 |
| Cross-validation | 100% |
| Statistical significance | 12.74σ |
| Predictions | 36 (21 confirmed) |

---

## Publications

- **Technical article (v80.0)**: [10.5281/zenodo.23304779](https://doi.org/10.5281/zenodo.23304779)
- **Review article (v80.0, Parts I–IV)**: [10.5281/zenodo.23305721](https://doi.org/10.5281/zenodo.23305721)
- **All versions**: [10.5281/zenodo.23305720](https://doi.org/10.5281/zenodo.23305720)

---

## Repository structure

```
omega4-unified-v80/
├── paper/       ← LaTeX sources + PDF (EN/RU)
├── python/      ← Python implementation (41 modules)
├── results/     ← JSON results
├── figures/     ← Figures (PNG)
├── docs/        ← Extended description (RU/EN)
├── README.md
├── CITATION.cff
└── LICENSE
```

---

## Quick start

```bash
git clone https://github.com/SergeyMatershov/omega4-unified-v80.git
cd omega4-unified-v80
pip install numpy scipy mpmath
python python/main_unified.py
```

---

## Honest disclaimer

The model is an **effective field theory (EFT)**, not a string theory in the strict sense.
PROVEN has not been achieved. The relation N_eff = 7 is an observed regularity, not a rigorous derivation.
All structural results are 1-loop proportions. Five open problems and six questions about the spread formula are explicitly listed in the paper.

---

## Author

**Sergey Viktorovich Matershov**
ORCID: [0009-0009-0641-1357](https://orcid.org/0009-0009-0641-1357)
Email: smatershov@mail.ru

---

## License

Universal Academic Licence v1.1 (UAL-v1.1) or CC BY-NC-ND 4.0.
See [LICENSE](LICENSE) for details.
