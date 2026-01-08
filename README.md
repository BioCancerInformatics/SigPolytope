# SigPolytope  
### Geometric Multidimensional Representation of Multi-Omic Signatures

**SigPolytope** is a geometric framework for the **representation, exploration, and interpretation of multi-omic signatures as structured multidimensional entities**.

This repository hosts:

- 📦 the **R package (`SigPolytope`)** implementing the geometric formalism
- 🧭 the **SigPolytope Shiny application** for interactive exploration and interpretation
- 📊 reproducible code and assets supporting the manuscript figures and atlas outputs

SigPolytope operationalizes a central conceptual shift:  
**omic signatures are not vectors, scores, or gene lists — they are structured, multidimensional informational entities whose biological meaning is inherently geometric**.

---

## Figure 1 — Nomenclature and regulatory circuitry concept

<!-- INSERT FIGURE 1 HERE (manuscript Figure 1) -->
<!-- Suggested repo path (recommended): docs/assets/Figure1.png -->
<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Figure%201.png" width="1000">
</p>

---

## Scientific motivation

### Why vectorial representations fail

Most biomarker and signature analyses implicitly assume that biological meaning can be captured by:

- gene lists  
- weighted scores  
- correlation vectors  
- network edges  

These representations are **reductionist**. They collapse coordinated, cross-layer biological organization into one-dimensional summaries that cannot preserve:

- regulatory alignment or opposition  
- phenotypic directionality  
- latent heterogeneity  
- multi-axis deformation  
- concordance versus discordance between interacting programs  

As a result, signatures that are biologically distinct often appear equivalent, while signatures that are mechanistically aligned may appear unrelated.

---

## Core concept: signatures have geometry

SigPolytope is based on a rigorous reconceptualization:

> **An omic signature is a multidimensional informational entity whose identity emerges from the coordinated organization of molecular, phenotypic, immune, microenvironmental, and clinical dimensions.**

Such entities cannot be faithfully represented as vectors.  
They must be represented as **geometric objects**.

---

## Geometric framework

### Latent multidimensional embedding

Each signature is encoded as an **18-dimensional latent regulatory vector** integrating:

- correlation structure  
- tumor–normal directionality  
- survival associations (OS, DSS, DFI, PFI)  
- Cox/log-rank strength and direction  
- microenvironmental scores  
- immune context and classification  

All latent dimensions are globally scaled and embedded into a **shared coordinate system**, allowing signatures to be compared structurally rather than descriptively.

For visualization, the latent space is projected into three dimensions while preserving geometric relationships.

### From points to polytopes

A signature is not represented by a single point.

Instead:

- each latent coordinate is perturbed symmetrically  
- a local multidimensional neighborhood is generated  
- the minimal enclosing envelope is computed  

The result is a **convex polytope**.

> The polytope is not a visualization artifact — **it is the signature**, expressed as a measurable geometric object.

---

## Figure 2 — Conceptual geometric representations (simulated data)

<!-- INSERT FIGURE 2 HERE (manuscript Figure 2; composite 2x3 panels) -->
<!-- Suggested repo path (recommended): docs/assets/Figure2_Composite_2x3.png -->
<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Figure%202_Composite_2x3_600dpi.png" width="1000">
</p>

---

## Controlled vocabulary (Box 1)

To ensure interpretability and terminological consistency, SigPolytope adopts a controlled geometric vocabulary that defines the core constructs used across the manuscript, R package, and Shiny interface.

<!-- INSERT BOX 1 IMAGE HERE (Box explaining terms) -->
<!-- Suggested repo path (recommended): docs/assets/Box1_SigPolytope.png -->
<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Box1_SigPolytope_600dpi_rectangular_lines.png" width="1000">
</p>

---

## Regulatory circuitries as dual polytopes

SigPolytope focuses on **metabolic regulatory circuitries**, defined as paired entities composed of:

- **Regulatory signature (sig)**  
  Upstream regulators (e.g., TFs, miRNAs, lncRNAs, CNVs, methylation)

- **Interaction signature (int)**  
  Downstream metabolic or functional programs

Each circuitry is therefore represented as **two convex polytopes embedded in the same latent space**.

This dual representation preserves:

- alignment versus opposition  
- balance versus asymmetry  
- shared versus divergent phenotypic axes  

---

## Intrinsic geometric descriptors

From the dual-polytope construction, SigPolytope derives **intrinsic, interpretable geometric measures**:

- **Barycenter distance**  
  Multidimensional regulatory concordance vs. discordance

- **Convex hull volume**  
  Latent dimensional complexity (single-axis vs. multi-axis behavior)

- **Volume asymmetry ratio**  
  Balance versus dominance between regulatory and interacting components

Together, these define **geometric phenotypes** that cannot be inferred from gene overlap, networks, or scalar statistics.

---

## What SigPolytope enables

SigPolytope allows users to:

- Explore **large atlases of regulatory circuitries** in a unified geometric space  
- Interactively visualize **paired 3D polytopes** (signature × interaction)  
- Quantify regulatory concordance, discordance, and asymmetry  
- Identify low-, intermediate-, and high-complexity multidimensional architectures  
- Compare signatures based on **structure rather than molecular composition**  
- Detect redundancy, divergence, and mechanistic equivalence geometrically  
- Support principled stratification and prioritization of candidate biomarkers  

Geometry becomes an **analytic and interpretive foundation**, not a cosmetic layer.

---

## Components in this repository

### 📦 R package: `SigPolytope`

The R package provides the analytical core:

- latent tensor construction (18D)
- polytope generation and convex hull computation
- geometric descriptors (distance, volume, asymmetry, anisotropy)
- reproducible figure generation

> Package source lives in: `R/`, `man/`, `NAMESPACE`, `DESCRIPTION`, `inst/`

### 🧭 Shiny application: SigPolytope Atlas

The Shiny app provides interactive exploration and interpretation of the geometric atlas:

- filtering by cancer type, pathway, layers, immune states
- dual-polytope 3D visualization (Plotly/WebGL)
- structured summaries linking geometry to biological context

> App source lives in: `shiny/` (or `inst/shiny/`, depending on your chosen structure)

**Live app:** https://sigpolytope.shinyapps.io/geometricatlas/

---

## 3D polytope visualization (HTML)

SigPolytope generates **interactive Plotly/WebGL HTML figures** representing dual polytopes.

Because GitHub READMEs cannot embed interactive HTML, previews are provided as static or animated images linking to the full figure:

```md
[![3D polytope preview](docs/assets/polytope_preview.gif)](https://biocancerinformatics.github.io/SigPolytope/www/polytope_LGG_6708_LGG_5904_dual_hulls_18D.html)
