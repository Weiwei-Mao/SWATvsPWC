---
marp: true
paginate: true
size: 16:9
style: |
  section {
    font-size: 20px;
  }
  table {
    font-size: 18px;
  }
  .col-container {
    display: flex;
    gap: 40px;
  }
  .col {
    flex: 1;
  }
  .col-title {
    font-weight: bold;
    font-size: 22px;
  }
---

<!-- _paginate: false -->

# SWAT+ vs PWC
## A Comparative Analysis of Pesticide Transport Models

<br><br>



---

<!-- _footer: "" -->

# Overview

| Feature | SWAT+ | PWC |
|---------|-------|-----|
| **Developer** | USDA-ARS + Texas A&M | US EPA |
| **Open Source** | Fully open | Partial (from EPA) |
| **Spatial Scale** | Watershed/HRU | Field/Single Point |
| **Time Step** | Daily | Sub-daily |
| **Purpose** | Watershed modeling | Risk assessment |

---

<div class="col-container">

<div class="col">

## PWC

- **Chemicals**: Maximum 3 (1 parent + 2 metabolite generations)
- **Metabolite transport**: Full (volatilization, sorption, runoff, leaching, erosion, plant uptake)

</div>

<div class="col">

## SWAT+

- **Chemicals**: Multiple pesticides simultaneously
- **Metabolite transport**: Full (volatilization, sorption, runoff, leaching, erosion, plant uptake)
- **Evolution**: Classic SWAT (parent only) → SWAT+ (metabolites ~2020) → 2025 plant uptake

</div>

</div>

---

<div class="col-container">

<div class="col">

<span class="col-title">PWC: Field & Water Body</span>

**Field (PRZM/TPEZ):**
- Sorption (Linear + Freundlich)
- Degradation (3-phase: aqueous/sorbed/gaseous)
- Volatilization (complete model)
- Transport: Runoff, erosion, leaching, plant uptake

**Water Body (VVWM):**
- Degradation, volatilization
- Settling, resuspension
- Benthic processes

</div>

<div class="col">

<span class="col-title">SWAT+: HRU & Channel</span>

**HRU:**
- Sorption (Linear only)
- Degradation (single-phase)
- Foliar/soil decay, washoff
- Transport: Surface runoff, lateral flow, tile drain, percolation, sediment

**Channel:**
- Reaction, metabolism (parent→daughter)
- Volatilization, settling/resuspension
- Diffusion, benthic reaction/burial

</div>

</div>

---

<div class="col-container">

<div class="col">

<span class="col-title">PWC: Application Methods</span>

- **Foliar application**
- **Soil surface** (4cm layer)
- **Uniform incorporation** (any depth)
- **Specific depth** (user-defined)
- **T-Band** (2cm band)
- **Linear decrease** with depth
- **Linear increase** with depth
- **Custom** distribution

</div>

<div class="col">

<span class="col-title">SWAT+: Application Methods</span>

- **Foliar** (LAI-based interception)
- **Soil surface**
- **Soil incorporation** (uniform mixing)

</div>

</div>

### PWC offers more refined application methods (8 vs 3)

---

<div class="col-container">

<div class="col">

<span class="col-title">PWC: Algorithm Detail</span>

- **Volatilization**: Complete canopy resistance model (Henry's Law + boundary layer)
- **Sorption**: Linear + Freundlich + non-equilibrium (predictor-corrector)
- **Numerical**: Predictor-corrector, sub-daily time steps, tridiagonal matrix solver

</div>

<div class="col">

<span class="col-title">SWAT+: Algorithm Detail</span>

- **Volatilization**: Simplified aquatic flux only (no canopy model)
- **Sorption**: Linear only (equilibrium assumption)
- **Numerical**: Explicit Euler, daily time steps, 1st-order decay simplification

</div>

</div>

---

<div class="col-container">

<div class="col">

<span class="col-title">PWC Output</span>

- **Field scale**: Daily fluxes for runoff, erosion, volatilization, degradation, washoff, plant uptake
- **Water body**: Estimated Environmental Concentration (EEC)
- **Purpose**: Regulatory risk assessment

</div>

<div class="col">

<span class="col-title">SWAT+ Output</span>

- **HRU scale**: Mass balance for plant, soil, sediment, runoff, lateral flow, tile drain, percolation
- **Channel scale**: Dissolved/sorbed transport, reaction, metabolism, volatilization, settling
- **Purpose**: Watershed management analysis

</div>

</div>

---

<!-- _footer: "" -->

# Scenario Comparison

| Scenario | SWAT+ | PWC | Recommended |
|----------|-------|-----|-------------|
| Watershed load | ✓ | ✗ | **SWAT+** |
| Risk assessment | △ | ✓ | **PWC** |
| Management eval | ✓ | △ | **SWAT+** |
| Regulatory compliance | ✗ | ✓ | **PWC** |
| Short-term events | △ | ✓ | **PWC** |
| Long-term trends | ✓ | △ | **SWAT+** |
| Spatial distribution | ✓ | ✗ | **SWAT+** |
| Multiple pollutants | ✓ | ✗ | **SWAT+** |
| Volatile pesticides | △ | ✓ | **PWC** |
| Nonlinear sorption | ✗ | ✓ | **PWC** |

Note: ✓ = Excellent, △ = Limited, ✗ = Not supported

---

<div class="col-container">

<div class="col">

<span class="col-title">PWC</span>

- **Chemicals**: 1 parent + 2 metabolites
- **Scale**: Single point (field)
- **Volatilization**: Complete canopy resistance model
- **Sorption**: Linear + nonlinear (Freundlich)
- **Degradation**: Three-phase (aqueous/sorbed/gaseous)
- **Application**: 8 methods
- **Time step**: Sub-daily
- **Regulatory**: ✓ EPA accepted

</div>

<div class="col">

<span class="col-title">SWAT+</span>

- **Chemicals**: Multiple pesticides
- **Scale**: Distributed (watershed/HRU)
- **Volatilization**: Simplified aquatic flux only
- **Sorption**: Linear only
- **Degradation**: Single-phase
- **Application**: 3 methods
- **Time step**: Daily
- **Regulatory**: ✗ Not EPA accepted

</div>

</div>

---

<!-- _paginate: false -->
<!-- _footer: "" -->

# Thank You!

<br>

### Questions & Discussion
