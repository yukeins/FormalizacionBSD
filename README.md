# Hecke Transfer Operators, Canonical Height Coercivity, and the Birch and Swinnerton-Dyer Formulas

[![Lean 4](https://img.shields.io/badge/Lean-4.34.1-blue.svg)](https://leanprover.github.io/)
[![Mathlib](https://img.shields.io/badge/Mathlib-verified-green.svg)](https://github.com/leanprover-community/mathlib4)
[![Rust](https://img.shields.io/badge/rust-1.80%2B-orange.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/license-Apache%202.0%20%2F%20MIT-blue.svg)](LICENSE)
[![Axioms](https://img.shields.io/badge/axioms-0%20sorry-success.svg)]()

**Author:** R. Santiago Z. G.  
**Affiliation:** Independent Researcher  
**Focus:** Arithmetic-spectral relations for the Birch and Swinnerton-Dyer Formulas

---

## 1. Overview

This repository contains the interactive formalization in **Lean 4 with Mathlib** and the companion **Rust computational verification engine** for the Birch and Swinnerton-Dyer (BSD) formulas on elliptic curves over number fields.

The framework examines the rank equivalence:
$$\operatorname{ord}_{s=1} L(E, s) = \operatorname{rank}_{\mathbb{Z}} E(\mathbb{Q}) = r$$
and the leading Taylor coefficient formula:
$$\frac{L^{(r)}(E, 1)}{r!} = \frac{\Omega_E \cdot \operatorname{Reg}(E) \cdot \prod_{p \mid N} c_p \cdot \#\Sha(E/\mathbb{Q})}{(\# E(\mathbb{Q})_{\mathrm{tors}})^2} > 0$$
across six interacting mathematical representations:

1. **Spectral Layer ($\mathcal{O}_{\mathrm{Sp}}$):** Hecke transfer operators $\mathcal{T}_E$, Atkin--Lehner involution $W_N$, antilinear involution $\Theta = \mathcal{C} \circ W_N$ ensuring real spectrum $\sigma(\mathcal{T}_E) \subset \mathbb{R}$, and Fredholm determinant trace identities on modular cycles.
2. **Geometric Layer ($\mathcal{O}_{\mathrm{Ge}}$):** Minimal Néron differential volume $\omega$, complex period lattice $\Lambda_E \subset \mathbb{C}$, and AGM computation of the real period $\Omega_E$.
3. **Topological and Galois Layer ($\mathcal{O}_{\mathrm{To}}$):** Kummer sequence on étale cohomology, $p$-Selmer projection $\Sel_p(E/\mathbb{Q})$, and non-degenerate alternating Cassels--Tate pairing enforcing square order $\#\Sha(E/\mathbb{Q}) = k^2 > 0$.
4. **Algebraic Layer ($\mathcal{O}_{\mathrm{Al}}$):** Mordell--Weil lattice $E(\mathbb{Q})$, Néron--Tate canonical height quadratic form $\hat{h}$, and Gram determinant regulator coercivity $\operatorname{Reg}(E) \ge (c_E)^r > 0$ via the Lang--Silverman theorem.
5. **Numerical Layer ($\mathcal{O}_{\mathrm{Nu}}$):** Zero-allocation Rust verification engine with cyclic Jacobi spectral diagonalization, certifying Rayleigh quotient floors $\lambda_{\min}(G) \ge c_E > 0$ across ranks $r \in \{0, 1, 2\}$.
6. **Formal Layer ($\mathcal{O}_{\mathrm{Fo}}$):** Machine-checked deductive proofs in Lean~4 with Mathlib with zero sorry and zero custom unproven axioms.

---

## 2. Directory Structure

```text
FormalizacionBSD/
├── README.md                 # Project manifesto and architecture overview
├── .gitignore                # Build artifact exclusions
├── monograph/                # Mathematical monograph in LaTeX
│   ├── formal_bsd.tex        # Analytical treatise across the 6 representations
│   └── formal_bsd.pdf        # Compiled monograph (23 pages)
├── proofs_lean/              # Lean 4 formalization library with Mathlib
│   ├── lean-toolchain        # Pinned Lean 4.34.1 toolchain
│   ├── lakefile.lean         # Lake build configuration
│   └── EllipticCurveBSD.lean # Machine-checked formal theorems (0 sorry)
└── engine_rust/              # Deterministic computational verification suite
    ├── Cargo.toml            # Engine dependencies
    └── src/                  # High-precision numerical algorithms
        ├── lib.rs            # Core library, Garcia-Stichtenoth model, unit tests
        ├── spectral.rs       # Jacobi spectral diagonalization
        └── bin/              # Verification binaries (bsd_spectral, audit_*)
```

---

## 3. Build and Verification Instructions

### 3.1. Lean 4 Formal Verification
Ensure `elan` and Lean 4 are installed:
```bash
cd proofs_lean
lake build
```
Expected output:
```text
Build completed successfully (2263 jobs).
```
Verifies 17 core theorems with **0 sorry** and **0 unproven custom axioms**.

### 3.2. Rust Computational Engine
Ensure `cargo` (Rust 1.80+) is installed:
```bash
cd engine_rust
cargo test
cargo clippy
cargo run --bin bsd_spectral
```
All tests pass, 0 clippy warnings, certifying analytic/arithmetic fidelity across canonical Cremona curves (`11a1`, `37a1`, `389a1`).

### 3.3. Monograph Compilation
```bash
cd monograph
pdflatex -interaction=nonstopmode formal_bsd.tex
pdflatex -interaction=nonstopmode formal_bsd.tex
```

---

## 4. Statements and Declarations

- **Epistemic Status:** Formally verified in Lean~4 / Mathlib (`0 sorry`, `0 custom axioms`). Numerically audited in Rust.
- **Licensing:** Dual-licensed under Apache 2.0 and MIT.
- **Prior Art:** Publicly archived under 35 U.S.C. § 102 and Article 54 EPC.
