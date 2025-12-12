# SigPolytope

**SigPolytope** is a Shiny dashboard from **OncoMetabolismGPS – Geometry** for exploring **multi-omic regulatory circuitries as geometric objects**. Instead of treating signatures as isolated gene lists, each circuitry is represented as a **paired 3D polytope** (signature × interaction) embedded in a latent multidimensional space.

---

## 3D polytope preview

Below is an illustrative example of the kind of **3D polytope** rendered by the app (Plotly/WebGL):

![Example 3D polytope rendered by SigPolytope](docs/assets/polytope_example.png)

> **Tip:** After you run the app locally once, you can export an actual figure from your own dataset using the “Download 3D Polytope (.html)” button, or generate a static PNG as shown in the “Generate the preview image” section below.

---

## What you can do in SigPolytope

- **Filter regulatory circuitries** by cancer type and metabolic context (superfamily, pathways, metabolic cell death, regimes).
- **Visualize a selected circuitry in 3D** as a paired polytope (signature × interaction).
- Inspect **discordance** and **complexity** measures (e.g., barycenter distance, volume ratio, geometric regimes).
- Explore **network-based views** of meaningful interactions and interpretive summaries.
- Export interactive outputs (HTML) for sharing and reproducibility.

---

## App structure (high level)

- **Overview & Concept**: framework introduction and quick-start.
- **Regulatory Circuitries (3D)**: filters → selection → **3D polytope** visualization + summary panel.
- **Regulatory Circuitries (network)**: interactive **visNetwork** view + textual summary.
- **About & Citation**: citation details and team information.

---

## Requirements

- R (>= 4.2 recommended)
- Packages used include: `shiny`, `shinydashboard`, `plotly`, `geometry`, `dplyr`, `DT`, `visNetwork`, `ggraph`, `tidygraph`, `htmlwidgets`, and others depending on enabled modules.

---

## Data

This app expects data files under the `data/` directory.

Typical inputs:
- `data/All_data.rds` (pre-packed datasets)
- `data/Regulatory_circuitries_geometric.rds` (geometric circuitries input)
- `data/TCGA_Cancer_types.xlsx` (optional; recommended to convert to `.rds` for faster startup)

> For deployment on shinyapps.io (especially on free plans), **precomputing heavy objects** (tensor/embedding) and loading them from RDS is strongly recommended.

---

## Run locally

From the repository root:

```r
# In R:
shiny::runApp()
