# SigPolytope  
### Geometric Multidimensional Representation of Multi-Omic Signatures

**SigPolytope** is a Shiny application developed within the **OncoMetabolismGPS – Geometry** framework for the **geometric representation, exploration, and interpretation of multi-omic regulatory signatures**.

SigPolytope operationalizes a conceptual shift:  
**omic signatures are not vectors, scores, or gene lists — they are structured, multidimensional informational entities whose biological meaning is inherently geometric**.

This application provides the first interactive implementation of that paradigm.

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

---

### From points to polytopes

A signature is not represented by a single point.

Instead:

- each latent coordinate is perturbed symmetrically  
- a local multidimensional neighborhood is generated  
- the minimal enclosing envelope is computed  

The result is a **convex polytope**.

> The polytope is not a visualization artifact — **it is the signature**, expressed as a measurable geometric object.

---

## Regulatory circuitries as dual polytopes

SigPolytope focuses on **metabolic regulatory circuitries**, defined as paired entities composed of:

- **Regulatory signature (sig)**  
  Upstream regulators (e.g. TFs, miRNAs, lncRNAs, CNVs)

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

## 3D polytope visualization

SigPolytope generates **interactive Plotly/WebGL HTML figures** representing dual polytopes.

Because GitHub READMEs cannot embed interactive HTML, previews are provided as static or animated images linking to the full figure:

```md
[![3D polytope preview](docs/assets/polytope_preview.gif)](https://biocancerinformatics.github.io/SigPolytope/www/polytope_LGG_6708_LGG_5904_dual_hulls_18D.html)
