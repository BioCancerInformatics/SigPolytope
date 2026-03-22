# SigPolytope  
### Geometric Multidimensional Representation of Multi-Omic Signatures

**SigPolytope** is a geometric framework for the **representation, exploration, and interpretation of multi-omic signatures as structured multidimensional entities**.

This repository hosts:

- 🧭 the **SigPolytope Shiny Atlas** for interactive geometric exploration  
- 📊 reproducible code and assets supporting the manuscript figures and atlas outputs  
- 📁 data structures and computational pipelines implementing the geometric formalism  

SigPolytope operationalizes a central conceptual shift:  
**omic signatures are not vectors, scores, or gene lists — they are structured, multidimensional informational entities whose biological meaning is inherently geometric**.

## Conceptual overview of the framework

The figure below provides a **schematic overview of the SigPolytope analytical workflow**, from multi-omic signatures to latent geometric representation and polytope-based interpretation.

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Figure_SigPolytope_large_near_A4_300dpi.png" width="1200">
</p>

## Scientific motivation

### Limits of vectorial representations

Multi-omic signatures are widely used in biomarker discovery, yet they are typically treated as vectors or composite scores, collapsing intrinsically multidimensional biological organization into one-dimensional summaries :contentReference[oaicite:1]{index=1}.

These representations fail to preserve:

- internal structure  
- regulatory directionality  
- multidimensional dependencies  
- latent heterogeneity  

This limitation motivates representations that preserve **structural organization rather than marginal associations**.

## Core concept: signatures have geometry

> **An omic signature is a multidimensional informational entity whose identity emerges from coordinated organization across biological layers.**

Signatures:

- are not defined by gene membership alone  
- encode structured relationships across molecular, phenotypic, immune, and clinical dimensions  
- require representations that preserve internal organization :contentReference[oaicite:2]{index=2}  


## Geometric framework

### Latent multidimensional embedding

Each signature is embedded into a **shared 18-dimensional latent space**, integrating:

- correlation structure  
- tumor–normal directionality  
- survival associations  
- immune and microenvironmental context  

The 3D representation is a **projection**, not the intrinsic space.

### From points to polytopes

Signatures are represented as **convex polytopes**, not points:

- latent coordinates are perturbed  
- a multidimensional neighborhood is generated  
- a convex hull defines the structure  

> The polytope is the signature — not a visualization artifact.

## Nomenclature and regulatory circuitry concept

The framework formalizes regulatory circuitries as structured paired entities.

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Figure%201.png" width="1000">
</p>

## Conceptual geometric representations

Distinct internal organizations generate distinct geometric structures.

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Figure%202_Composite_2x3_600dpi.png" width="1000">
</p>

## Controlled vocabulary (Box 1)

To ensure interpretability, the framework adopts a controlled geometric vocabulary.

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Box1_SigPolytope_600dpi_rectangular_lines.png" width="1000">
</p>

## Regulatory circuitries as dual polytopes

Each circuitry is represented as:

- a **regulatory polytope (sig)**  
- an **interaction polytope (int)**  

embedded in the same latent space.

## Intrinsic geometric descriptors

SigPolytope derives:

- **Barycenter distance** → concordance vs discordance  
- **Convex hull volume** → latent complexity  
- **Volume asymmetry** → regulatory balance  
- **Anisotropy** → dominant axes  

These properties are not captured by vector-based methods.

## What SigPolytope enables

- structural comparison of signatures  
- detection of redundancy and divergence  
- quantification of multidimensional behavior  
- principled biomarker prioritization  

Geometry becomes an **analytical foundation**, not a visualization layer.

## Components in this repository

### 🧭 Shiny Atlas

Interactive exploration of geometric circuitries:

- dual polytope visualization  
- filtering across biological dimensions  
- structured interpretation  

Live app: https://sigpolytope.shinyapps.io/geometricatlas/

### 📊 Computational pipelines

Includes:

- latent vector construction  
- polytope generation  
- geometric descriptors  
- manuscript figure generation  

All results are **deterministically reproducible** :contentReference[oaicite:3]{index=3}.

## Conceptual scope

SigPolytope is a **representational framework**, not a predictive model.

It provides a geometric basis for:

- comparison  
- interpretation  
- structural analysis of multi-omic systems  

## Citation

Please cite the corresponding manuscript when available.
