# SWAT+ vs PWC Pesticide Model Comparison (PPT Outline)

---

## Slide 1: Title

**SWAT+ vs PWC: A Comparative Analysis of Pesticide Transport Models**
- Based on Source Code Analysis
- Physical Processes and Algorithm Differences
- Model Selection Recommendations

---

## Slide 2: Model Overview

| Feature | SWAT+ | PWC |
|---------|-------|-----|
| **Developer** | USDA-ARS + Texas A&M | US EPA |
| **Open Source** | Fully open | Partial (software from EPA, incomplete source) |
| **Spatial Scale** | Watershed/Subbasin/HRU (Distributed) | Field/Single Point |
| **Time Step** | Daily | Sub-daily supported |
| **Purpose** | Integrated watershed modeling | Pesticide risk assessment |

---

## Slide 3: Pesticides per Simulation

### PWC
- **Maximum 3 chemicals**
- Structure: 1 parent + 2 generations of metabolites
  - Parent (mother compound)
  - Daughter (1st generation metabolite)
  - Granddaughter (2nd generation metabolite)

### SWAT+
- **Multiple pesticides simultaneously**
- Each pesticide can have multiple metabolites

**Conclusion**: SWAT+ can simulate more pesticide types simultaneously

---

## Slide 4: Physical Processes - PWC

### Field Processes (PRZM/TPEZ)

```
Pesticide Application
    │
    ▼
┌────────────────────────────────────────────────────────┐
│                    Soil Profile                          │
├────────────────────────────────────────────────────────┤
│  Foliar Washoff → Sorption → Degradation → Volatilization │
│       │              │            │             │       │
│       ▼              ▼            ▼             ▼       │
│    WOFLUX        Sorbed      DKFLUX       PVFLUX           │
│   (washoff)    (adsorbed)  (3-phase)     (volatilization)   │
│                                                        │
│       └────────────────┴──────────────┴───────────────┘   │
│                          │                               │
│                          ▼                               │
│              ┌─────────────────────┐                      │
│              │   Transport         │                      │
│              ├─────────────────────┤                      │
│              │ • Runoff   ROFLUX   │                      │
│              │ • Erosion  ERFLUX   │                      │
│              │ • Leaching DCOFLUX  │                      │
│              │ • Uptake   UPFLUX   │                      │
│              └─────────────────────┘                      │
└────────────────────────────────────────────────────────┘
```

### Three-Phase Degradation

| Phase | Rate | Description |
|-------|------|-------------|
| Aqueous | DWRATE | Dissolved phase degradation |
| Sorbed | DSRATE | Soil particle adsorbed phase degradation |
| Gaseous | DGRATE | Soil air phase degradation |

---

## Slide 5: Physical Processes - SWAT+

### HRU-Scale Processes

```
Pesticide Application
    │
    ▼
┌────────────────────────────────────────────────────────┐
│                    HRU (Field Unit)                       │
├────────────────────────────────────────────────────────┤
│  Foliar decay → Soil decay → Washoff → Plant uptake      │
│                                                        │
│              ┌─────────────────────┐                     │
│              │  Transport Output   │                     │
│              ├─────────────────────┤                     │
│              │ • Surface runoff surq│                     │
│              │ • Lateral flow  latq│                     │
│              │ • Tile drain   tileq│                     │
│              │ • Percolation  perc │                     │
│              │ • Sediment     sed  │                     │
│              └─────────────────────┘                     │
└────────────────────────────────────────────────────────┘
```

### Channel Processes

```
┌────────────────────────────────────────────────────────┐
│                      Channel                             │
├────────────────────────────────────────────────────────┤
│  • Reaction      (aqueous degradation)                 │
│  • Metabolism     (parent → daughter)                  │
│  • Volatilization                                      │
│  • Settling                                           │
│  • Resuspension                                       │
│  • Diffusion       (benthic-water exchange)            │
│  • Benthic reaction                                   │
│  • Burial                                            │
└────────────────────────────────────────────────────────┘
```

---

## Slide 6: Application Methods

### PWC: 8 Application Methods

| CAM | Method | Description |
|-----|--------|-------------|
| 1 | Soil surface | Linear decrease to 4cm |
| 2 | Foliar | Canopy interception + ground decrease |
| 3 | Uniform incorporation | Uniform distribution to specified depth |
| 4 | Specific depth | All at specified depth |
| 5 | T-Band | Top 2cm concentrated +下层uniform |
| 6 | Linear decrease | Linear decrease with depth |
| 7 | Linear increase | Linear increase with depth |
| 8 | Custom | User-defined mixing depth |

### SWAT+: 3 Basic Methods
- Foliar application (LAI-based distribution)
- Soil surface application
- Soil incorporation

**Conclusion**: PWC has more refined application methods

---

## Slide 7: Volatilization

### PWC: Complete Canopy Resistance Model

```
Soil Surface → Henry's Law: C_air = C_water × H
                        │
                        ▼
                   Air diffusion
                   (porosity correction)
                        │
                        ▼
                Boundary layer resistance
                        │
                        ▼
            Canopy resistance (with vegetation)
                        │
                        ▼
              Flux = -CONDUC × C × H
```

### SWAT+: Simplified Volatilization

```
Aquatic volatilization
Flux = -aq_volat × C
```

**Conclusion**: PWC has more complete volatilization model

---

## Slide 8: Sorption

| Feature | PWC | SWAT+ |
|---------|-----|-------|
| Linear sorption | S = Kd × C | S = Kd × C |
| Freundlich nonlinear | S = Kf × C^N | ❌ Not supported |
| Nonequilibrium | Two-domain model | ❌ Not supported |
| Numerical treatment | Predictor-corrector | Linearized |

**Conclusion**: PWC has more complex sorption model

---

## Slide 9: Numerical Methods

| Aspect | SWAT+ | PWC |
|--------|-------|-----|
| **Time integration** | Explicit Euler | Predictor-corrector |
| **Sub-daily step** | ❌ Daily only | ✅ Supported |
| **Advection-dispersion** | Simplified 1st-order | Tridiagonal matrix |
| **Nonlinear treatment** | Linearization | Iterative solution |

---

## Slide 10: Output Comparison

### PWC Output

```
Fluxes (g/cm²/day):
  ROFLUX  - Runoff
  ERFLUX  - Erosion
  PVFLUX  - Volatilization
  DKFLUX  - Degradation
  WOFLUX  - Washoff
  UPFLUX  - Plant uptake
  DCOFLUX - Bottom outflow

Concentration:
  EEC - Estimated Environmental Concentration (µg/L)
```

### SWAT+ Output

```
Mass balance:
  plant, soil, sed, surq, latq, tileq, perc
  apply_s, apply_f, decay_s, decay_f, wash
  metab_s, metab_f, pl_uptake

Channel:
  Dissolved/Sorbed outflow
  Reaction, Metabolism, Volatilization, Settling, Resuspension, Diffusion
```

---

## Slide 11: Scenario Comparison Matrix

| Scenario | SWAT+ | PWC | Recommended |
|----------|-------|-----|-------------|
| Watershed pesticide load | ✅ | ❌ | **SWAT+** |
| Pesticide risk assessment | ⚠️ | ✅ | **PWC** |
| Management practice evaluation | ✅ | ⚠️ | **SWAT+** |
| Regulatory compliance | ❌ | ✅ | **PWC** |
| Short-term events | ⚠️ | ✅ | **PWC** |
| Long-term trends | ✅ | ⚠️ | **SWAT+** |
| Spatial distribution | ✅ | ❌ | **SWAT+** |
| Multiple pollutants | ✅ | ❌ | **SWAT+** |
| Volatile pesticides | ⚠️ | ✅ | **PWC** |
| Nonlinear sorption | ❌ | ✅ | **PWC** |

---

## Slide 12: Key Differences Summary

| Aspect | SWAT+ | PWC |
|--------|-------|-----|
| **# of pesticides** | Multiple | 1+2 metabolites |
| **Spatial scale** | Distributed | Single point |
| **Volatilization** | Simplified | Complete model |
| **Sorption** | Linear | Linear + nonlinear |
| **Degradation** | Single-phase | Three-phase |
| **Application methods** | 3 types | 8 types |
| **Time step** | Daily | Sub-daily supported |
| **Regulatory recognition** | ❌ | EPA ✅ |

---

## Slide 13: Q&A

**Q1: Can both models be used together?**
A1: Yes, nested application - SWAT+ for watershed analysis, PWC for detailed point assessment

**Q2: Which is more accurate?**
A2: Depends on assessment objectives. PWC is more accurate for field concentration, SWAT+ is more comprehensive for watershed loading

**Q3: How to choose?**
A3: Consider research objectives, data availability, and regulatory acceptance requirements

---

*End of document*
