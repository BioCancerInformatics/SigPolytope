# SigPolytope  
### Geometric Multidimensional Representation of Multi-Omic Signatures

**SigPolytope** is a geometric framework for the **representation, exploration, and interpretation of multi-omic signatures as structured multidimensional entities**.

This repository hosts:

- 🧭 the **SigPolytope Shiny Atlas** for interactive geometric exploration  
- 📊 reproducible code and assets supporting the manuscript figures and atlas outputs  
- 📁 data structures and computational pipelines implementing the geometric formalism  

SigPolytope operationalizes a central conceptual shift:  
**omic signatures are not vectors, scores, or gene lists — they are structured, multidimensional informational entities whose biological meaning is inherently geometric**.


---

## Visual entry point to the framework

The figure below is the **conceptual gateway** to SigPolytope.  
It synthesizes the analytical progression established in the manuscript, from:

- multi-omic data integration  
- regulatory circuitry construction  
- latent geometric embedding  
- multidimensional polytope representation  

This image is not included as a decorative element. It is the **entry-level abstraction of the full framework**, and it is also the basis for the redesigned Shiny entrance interface.

<!-- Recommended path for the updated entry-level figure -->
<!-- Example: SigPolytope Shiny/www/Entry_Interface_Framework.png -->

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Figure_SigPolytope_large_near_A4_300dpi.png" width="1200" alt="SigPolytope conceptual entry figure">
</p>

> In the Shiny application, this figure functions as the **full-page entrance interface**, from which users transition into the interactive SigPolytope environment.


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

**For visualization, the 18D latent space is projected into three dimensions.**  
Importantly, the geometric representation is **not intrinsically 3D**: each signature is an **18D polytope**, and the Shiny Atlas displays a **3D projection** of this high-dimensional structure for interpretability.


### From points to polytopes

A signature is not represented by a single point.

Instead:

- each latent coordinate is perturbed symmetrically  
- a local multidimensional neighborhood is generated  
- the minimal enclosing envelope is computed  

The result is a **convex polytope in the 18D latent space**.

> The polytope is not a visualization artifact — **it is the signature**, expressed as a measurable high-dimensional geometric object.  
> The interface displays a **3D projection** of this object for exploration.


---

## Nomenclature and regulatory circuitry concept

The framework formalizes regulatory circuitry as a structured geometric system in which molecular interactions are not merely catalogued, but embedded into a latent multidimensional space.

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Figure%201.png" width="1000" alt="SigPolytope nomenclature and regulatory circuitry concept">
</p>


---

## Conceptual geometric representations

These panels illustrate how latent multidimensional organization can be translated into interpretable geometric forms.

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Figure%202_Composite_2x3_600dpi.png" width="1000" alt="Conceptual geometric representations in SigPolytope">
</p>


---

## Controlled vocabulary (Box 1)

To ensure interpretability and terminological consistency, SigPolytope adopts a controlled geometric vocabulary that defines the core constructs used across the manuscript and the Shiny Atlas.

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Box1_SigPolytope_600dpi_rectangular_lines.png" width="1000" alt="Controlled vocabulary for SigPolytope">
</p>


---

## Regulatory circuitries as dual polytopes

SigPolytope focuses on **metabolic regulatory circuitries**, defined as paired entities composed of:

- **Regulatory signature (sig)**  
- **Interaction signature (int)**  

Each circuitry is therefore represented as **two convex polytopes embedded in the same 18D latent space**.

This dual representation preserves:

- alignment versus opposition  
- balance versus asymmetry  
- shared versus divergent phenotypic axes  

For interpretability, the Shiny Atlas renders **paired 3D projections of the two 18D polytopes** in a shared coordinate system.


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
- Enter the framework through a **full-page conceptual interface** aligned with the manuscript logic  
- Interactively visualize **paired 3D projections of 18D polytopes** (signature × interaction)  
- Quantify regulatory concordance, discordance, and asymmetry  
- Identify low-, intermediate-, and high-complexity multidimensional architectures  
- Compare signatures based on **structure rather than molecular composition**  
- Detect redundancy, divergence, and mechanistic equivalence geometrically  
- Support principled stratification and prioritization of candidate biomarkers  

Geometry becomes an **analytic and interpretive foundation**, not a cosmetic layer.


---

## Components in this repository

### 🧭 SigPolytope Shiny Atlas

The Shiny application provides interactive access to the geometric atlas through a staged interface:

- **Level 1**: full-page conceptual entry figure  
- **Level 3**: interactive SigPolytope exploration environment  

Main features include:

- filtering by cancer type, pathway, omic layers, immune states  
- **dual-polytope 3D projections** of **18D geometric representations** (Plotly/WebGL)  
- structured summaries linking geometry to biological context  
- exportable figures and summary tables  

> App source lives in: `SigPolytope Shiny/`

**Live app:** https://sigpolytope.shinyapps.io/geometricatlas/


### 📊 Reproducible computational pipelines

This repository also contains:

- latent tensor construction routines (18D)  
- polytope generation and convex hull computation  
- geometric descriptors (distance, volume, asymmetry, anisotropy)  
- code for generating manuscript and atlas figures  
- static assets used in the repository and Shiny entry interface  

These pipelines implement the SigPolytope formalism in a transparent, reproducible manner.
