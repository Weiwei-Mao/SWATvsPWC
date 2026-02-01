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
- **Evolution**: Classic SWAT (parent only) → SWAT+ metabolites (Rathjens et al., 2022) → plant uptake (Rathjens et al., 2025)

</div>

</div>

---

<div class="col-container">

<div class="col">

<span class="col-title">PWC: Field & Water Body</span>

**Field (PRZM/TPEZ):**
- Sorption (Linear + Freundlich)
- Degradation (3-phase: aqueous/sorbed/gaseous)
- Volatilization (canopy resistance model)
- Transport: Runoff, erosion, leaching, plant uptake

**Water Body (VVWM):**
- Degradation, volatilization
- Benthic processes (no resuspension)

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
- **Soil incorporation** (uniform mixing to specified depth)

</div>

</div>

### PWC offers more refined application methods (8 vs 3)

---

<div class="col-container">

<div class="col">

<span class="col-title">PWC: Algorithm Detail</span>

- **Volatilization**: Canopy resistance model (Henry's Law + boundary layer)
- **Sorption**: Linear + Freundlich; optional non-equilibrium
- **Numerical**: Sub-daily time steps; tridiagonal solver

</div>

<div class="col">

<span class="col-title">SWAT+: Algorithm Detail</span>

- **Volatilization**: Simplified aquatic flux
- **Sorption**: Linear only (equilibrium assumption)
- **Numerical**: Daily time steps

</div>

</div>

---

<div class="col-container">

<div class="col">

<span class="col-title">PWC Output</span>

**Field (PRZM/TPEZ):**
- Daily fluxes: runoff, erosion, volatilization, degradation, drift, washoff
- Bottom outflow

**Water Body (VVWM):**
- Water concentration
- Benthic concentration

</div>

<div class="col">

<span class="col-title">SWAT+ Output</span>

**HRU:**
- Mass balance: plant, soil, sediment
- Transport: runoff, lateral flow, tile drain, percolation

**Channel:**
- Dissolved/sorbed transport
- Reaction, metabolism, volatilization, settling

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
- **Volatilization**: Canopy resistance model
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
- **Volatilization**: Simplified aquatic flux
- **Sorption**: Linear only
- **Degradation**: Single-phase
- **Application**: 3 methods
- **Time step**: Daily
- **Regulatory**: ✗ Not EPA accepted

</div>

</div>

---

# References

1. Rathjens, H., Kiesel, J., Miguez, M.B., Winchell, M., Arnold, J.G., Sur, R. (2022). Simulation of Pesticide and Metabolite Concentrations Using SWAT+ Landscape Routing and Conditional Management Applications. Water 14(9):1332. doi:10.3390/w14091332
2. Rathjens, H., Kiesel, J., Arnold, J., Reinken, G., Sur, R. (2025). Technical note: Extending the SWAT2012 and SWAT+ models to simulate pesticide plant uptake processes. Hydrol. Earth Syst. Sci. 29:6703–6714. doi:10.5194/hess-29-6703-2025

<!-- _paginate: false -->
<!-- _footer: "" -->

# Thank You!

<br>
