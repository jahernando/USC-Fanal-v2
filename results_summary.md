# FANAL — Comparative Results Summary

Results from the solution notebooks for each collaboration.

## Experimental Setup

| Parameter | Alpha | Beta | Gamma | Delta | Epsilon |
|-----------|------:|-----:|------:|------:|--------:|
| Exposure (kg·y) | 500 | 1000 | 1000 | 3000 | 3000 |
| σ\_E (keV) | 5.25 | 5.25 | 7.93 | 7.91 | 10.9 |
| Selection | 1 track, E\_blob2 > 0.4 MeV, E ∈ (2.4, 2.7) MeV | = | = | = | = |
| RoI (MeV) | (2.43, 2.48) | (2.43, 2.48) | (2.43, 2.48) | (2.43, 2.48) | (2.43, 2.48) |

## Expected Background in the RoI

| Source | Alpha | Beta | Gamma | Delta | Epsilon |
|--------|------:|-----:|------:|------:|--------:|
| n(Bi) expected | 2.31 | 21.73 | 9.97 | 66.68 | 57.14 |
| n(Tl) expected | 0.22 | 4.01 | 0.91 | 12.18 | 11.96 |
| **Total bkg expected** | **2.53** | **25.74** | **10.88** | **78.85** | **69.10** |

## Observed Data and Signal Extraction

| Quantity | Alpha | Beta | Gamma | Delta | Epsilon |
|----------|------:|-----:|------:|------:|--------:|
| **n observed in RoI** | **6** | **36** | **26** | **117** | **93** |
| n signal (count) = n\_obs − n\_bkg | 3.47 | 10.26 | 15.12 | 38.15 | 23.90 |
| n(ββ) from fit | 4.27 | 11.43 | 13.85 | 29.16 | 22.50 |

## Confidence Intervals on Signal

| Quantity | Alpha | Beta | Gamma | Delta | Epsilon |
|----------|------:|-----:|------:|------:|--------:|
| FC 68% CI (count) | (1.40, 6.72) | (4.65, 17.05) | (9.94, 20.86) | (27.69, 49.05) | (14.56, 33.97) |
| FC 95% CI (count) | (0.10, 10.13) | (1.03, 23.50) | (6.11, 26.86) | (18.20, 60.92) | (7.63, 43.68) |
| Profile 68% CI (fit) | (2.50, 6.88) | (7.43, 15.99) | (8.53, 19.77) | (19.69, 39.52) | (13.01, 31.98) |

## Statistical Significance

| Quantity | Alpha | Beta | Gamma | Delta | Epsilon |
|----------|------:|-----:|------:|------:|--------:|
| **z₀ (σ)** | **3.46** | **3.76** | **2.97** | **3.35** | **2.12** |
| p-value | 2.70×10⁻⁴ | 8.59×10⁻⁵ | 1.47×10⁻³ | 4.01×10⁻⁴ | 1.68×10⁻² |

## Half-life Measurements

| Method | Alpha | Beta | Gamma | Delta | Epsilon |
|--------|------:|-----:|------:|------:|--------:|
| t₁/₂ counting (×10²⁵ y) | 20.7 | 14.0 | 9.47 | 11.3 | 17.5 |
| t₁/₂ 68% CI counting (×10²⁵ y) | (10.7, 51.1) | (8.42, 30.9) | (6.87, 14.4) | (8.76, 15.5) | (12.3, 28.8) |
| t₁/₂ fit (×10²⁵ y) | 16.8 | 12.5 | 10.4 | 14.8 | 19.1 |
| t₁/₂ 68% CI fit (×10²⁵ y) | (10.4, 28.7) | (8.97, 19.3) | (7.26, 16.8) | (10.9, 21.9) | (13.5, 33.1) |

## Combined Result (Counting Experiment)

The five collaborations are combined using a joint Poisson likelihood with a common decay-rate parameter $\mu = 1/\mathcal{T}_{1/2}$:

$$\mathcal{L}(\mu) = \prod_c \mathrm{Poisson}(n_c \mid b_c + \mu \, s_c)$$

| Quantity | Combined | Range (individual) |
|----------|------:|------:|
| n observed (total) | 278 | 6 – 117 |
| n background (total) | 187.10 | 2.53 – 78.85 |
| n signal (total) | 90.90 | 3.47 – 38.15 |
| **t₁/₂ (×10²⁵ y)** | **13.28** | 9.47 – 20.7 |
| t₁/₂ 68% CI (×10²⁵ y) | (11.2, 16.1) | — |
| t₁/₂ 95% CI (×10²⁵ y) | (9.7, 20.0) | — |
| **z₀ (σ)** | **6.59** | 2.12 – 3.76 |
| p-value | 2.26×10⁻¹¹ | 1.68×10⁻² – 8.59×10⁻⁵ |

**Compatibility test**: χ² = 2.23, ndof = 4, p = 0.694 — all five measurements are consistent.

## Key Observations

- **Discovery**: No individual collaboration exceeds 5σ, but the combination reaches **6.59σ**, crossing the discovery threshold.
- **Best individual significance**: Beta (3.76σ), followed by Alpha (3.46σ) and Delta (3.35σ).
- **Worst significance**: Epsilon (2.12σ) — the poorest energy resolution (10.9 keV) dilutes the signal despite 3000 kg·y exposure.
- **Energy resolution matters**: Beta (5.25 keV, 1000 kg·y) outperforms Delta (7.91 keV, 3000 kg·y) in significance, demonstrating that better resolution is more effective than simply increasing exposure when backgrounds are high.
- **Signal excess**: All collaborations observe more events than the expected background, consistent with a ββ0ν signal.
- **Counting vs. fit**: The fit-based t₁/₂ is generally consistent with the counting estimate, with deviations depending on how well the 3-component fit separates signal from background.
- **Compatibility**: The five individual measurements are fully compatible (p = 0.69), validating the combination.
