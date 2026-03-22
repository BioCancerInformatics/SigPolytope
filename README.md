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

The figure below provides a **schematic overview of the SigPolytope analytical framework**, synthesizing the progression from:

- multi-omic signatures  
- regulatory circuitry construction  
- latent-space embedding  
- convex polytope representation  
- geometric descriptor extraction  

This schematic corresponds to the conceptual workflow underlying the framework and complements the formal description provided in the manuscript and supplementary materials.

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Figure_SigPolytope_large_near_A4_300dpi.png" width="1200" alt="SigPolytope conceptual framework">
</p>

## Scientific motivation

### Limits of vectorial representations

Multi-omic signatures are widely used in biomarker discovery, precision oncology, and systems biology. However, they are typically represented as:

- gene lists  
- composite scores  
- correlation vectors  
- network relationships  

These representations collapse intrinsically multidimensional biological organization into one-dimensional summaries. As a consequence, they fail to preserve:

- internal structure  
- directionality across phenotypic axes  
- regulatory alignment or opposition  
- latent heterogeneity  
- multidimensional dependencies across biological layers  

As emphasized in the manuscript, this reductionist framing limits both interpretability and comparability of omic signatures :contentReference[oaicite:0]{index=0}.

## Core concept: signatures as structured entities

SigPolytope is based on a rigorous reconceptualization:

> **An omic signature is a multidimensional informational entity whose biological meaning emerges from its internal structural organization across molecular, phenotypic, immune, microenvironmental, and clinical dimensions.**

In this formulation:

- signatures are not defined by molecular membership alone  
- biological meaning arises from **coordinated structure across layers**  
- similar gene lists may encode distinct regulatory architectures  
- distinct molecular compositions may converge to similar biological behavior :contentReference[oaicite:1]{index=1}  

This motivates a representational framework that preserves structure rather than collapsing it.

## Geometric framework

### Latent multidimensional embedding

Each signature is embedded into a **shared 18-dimensional latent space** integrating:

- correlation structure  
- tumor–normal directionality  
- survival associations (OS, DSS, DFI, PFI)  
- microenvironmental context  
- immune classification  
- statistical strength of associations  

Latent dimensions are globally scaled, enabling direct comparison across signatures.

For visualization purposes, this space is projected into three dimensions.  
However:

> **the intrinsic representation remains 18-dimensional, and all geometric properties are defined in this full latent space**.

### From vectors to polytopes

Instead of representing signatures as single points:

- each latent coordinate is symmetrically perturbed  
- a local multidimensional structure is generated  
- a convex envelope is computed  

The result is a **convex polytope**, representing the signature as a geometric object.

> The polytope is not a visualization artifact — it is the mathematical representation of the signature in latent space.

## Regulatory circuitries as paired geometric entities

SigPolytope focuses on **metabolic regulatory circuitries**, defined as paired entities composed of:

- a **regulatory signature (sig)**  
- an **interaction signature (int)**  

Each circuitry is represented as **two polytopes embedded in the same latent space**, preserving:

- concordance versus discordance  
- symmetric versus asymmetric structure  
- shared versus divergent phenotypic axes  

This dual-polytope representation enables direct geometric comparison of regulatory configurations.

## Intrinsic geometric descriptors

From this representation, SigPolytope derives intrinsic geometric measurements:

- **Barycenter distance**  
  → multidimensional concordance versus discordance  

- **Convex hull volume**  
  → latent dimensional complexity  

- **Volume asymmetry ratio**  
  → balance between regulatory and interaction components  

- **Anisotropy and orientation**  
  → dominant phenotypic axes  

These descriptors quantify structural properties that are not accessible through vector-based or network-based representations.

## Conceptual geometric representations

The following panels illustrate how distinct internal organizations produce distinct geometric forms in latent space.

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Figure%202_Composite_2x3_600dpi.png" width="1000" alt="Conceptual geometric representations">
</p>

## Controlled vocabulary (Box 1)

To ensure conceptual clarity and interpretability, the framework adopts a controlled geometric vocabulary defining core constructs used across analyses.

<p align="center">
  <img src="https://github.com/BioCancerInformatics/SigPolytope/blob/main/SigPolytope%20Shiny/www/Box1_SigPolytope_600dpi_rectangular_lines.png" width="1000" alt="Controlled vocabulary">
</p>

## What SigPolytope enables

SigPolytope enables:

- structural comparison of multi-omic signatures in a shared latent space  
- identification of concordant and discordant regulatory programs  
- quantification of latent complexity and multidimensional engagement  
- detection of redundancy and functional equivalence  
- principled prioritization of candidate biomarkers  
- interpretation of regulatory systems beyond molecular overlap  

Geometry becomes an **analytical and interpretive layer**, complementing—but not replacing—vector-based approaches.

## Components in this repository

### 🧭 SigPolytope Shiny Atlas

The Shiny application provides an interactive interface for exploring geometric representations of regulatory circuitries:

- filtering by cancer type, pathway, omic layers, immune states  
- visualization of **paired 3D projections of 18D polytopes**  
- structured summaries linking geometry to biological context  
- exportable figures and summary tables  

> App source: `SigPolytope Shiny/`  
> Live app: https://sigpolytope.shinyapps.io/geometricatlas/

### 📊 Reproducible computational pipelines

This repository includes:

- latent vector construction (18D)  
- convex polytope generation  
- geometric descriptor computation  
- benchmarking and validation routines  
- scripts used to generate manuscript figures  

All geometric representations and metrics are **deterministically derived** from the latent tensors, ensuring reproducibility :contentReference[oaicite:2]{index=2}.

## Conceptual scope

This framework is methodological and representational in nature.

It does not replace existing statistical or predictive models, but instead provides:

- a **structural representation of multi-omic signatures**  
- a **geometric basis for comparison and interpretation**  
- a **foundation for downstream analytical and translational applications**  

By treating signatures as geometric objects, SigPolytope reframes biological inference as a problem of **multidimensional structural organization** rather than scalar association.

## Citation

If you use SigPolytope, please cite the corresponding manuscript and repository resources.
