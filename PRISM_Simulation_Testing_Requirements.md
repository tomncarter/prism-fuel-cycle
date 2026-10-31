# PRISM Simulation Testing Requirements

---

## T1 — Criticality and Reactivity

**T1.1 — BOC k-eff within design band**
For each core configuration, beginning-of-cycle (BOEC) k-eff must fall within ±200 pcm of the target value (nominally k-eff ≈ 1.005 after control rod withdrawal allowance). Applies to startup U-Zr core, EG23 U-Pu-Zr equilibrium, and EG24 U-TRU-Zr equilibrium.

**T1.2 — EOC k-eff ≥ 1.0**
End-of-cycle (EOEC) k-eff must remain ≥ 1.0 (with statistical uncertainty band) to confirm the 12-month cycle length is achievable without excess burnup swing. Applies to all four core configurations.

**T1.3 — Cycle length reproducibility**
At equilibrium, the BOEC→EOEC reactivity swing must be stable (< 50 pcm variation) across three consecutive cycles, confirming true equilibrium has been reached.

**T1.4 — Conversion ratio agreement**
Computed conversion ratio (CR = fissile production / fissile destruction per cycle) must match the Triplett 2012 Table II values within ±5%: CR = 0.72 (UNF Recycle), 1.06 (Breakeven), 1.23 (Breeder), 1.00 (Weapons Pu). This is the primary benchmark against published data.

---

## T2 — Burnup and Fuel Performance

**T2.1 — Average discharge burnup**
Computed average discharge burnup must agree with Triplett 2012 Table II within ±5 MWd/kg HM: target 87.51 MWd/kg HM (UNF Recycle), 106.3 (Breakeven/Breeder/Weapons Pu). This cross-checks the depletion timestep and batch management logic.

**T2.2 — Peak discharge burnup**
Peak discharge burnup (hottest assembly in the outgoing batch) must not exceed the Triplett 2012 Table II values by more than 5%: 122.7 MWd/kg HM (UNF Recycle), ~147 MWd/kg HM (Breakeven/Breeder). Exceedance would indicate a shuffling pattern error producing unrealistic peaking.

**T2.3 — TRU weight fraction at BOC**
Fresh fuel TRU weight fraction (inner/outer zone) must match Triplett 2012 Table II within ±1 wt%: inner 18.9%, outer 22.7% (UNF Recycle); ~21% (Breakeven/Breeder). This verifies the material transfer and fresh fuel fabrication logic.

**T2.4 — In-vessel cooling period enforcement**
No discharged assembly may be transferred to the NFRC in fewer than one full irradiation cycle (12 months) from the cycle in which it left the active core. This is a hard logic requirement on the Python material-transfer script.

---

## T3 — Safety Parameters

**T3.1 — Doppler coefficient sign and magnitude**
The Doppler coefficient α_D must be negative throughout the transition (all phases, all core configurations). Computed value should fall in the range −0.003 to −0.001 Δk/k per °C (consistent with published IFR/PRISM values). A positive Doppler at any transition phase is a disqualifying result.

**T3.2 — Total power reactivity coefficient**
The total power reactivity coefficient (including Doppler, coolant expansion, and structural feedback) must be negative over the full operating power range for each transition phase. This is the fundamental passive safety requirement from NUREG-1368.

**T3.3 — Sodium void worth**
Sodium void worth (full-core voiding) must remain below the prompt criticality threshold (i.e., Δk_void < β_eff). The value is expected to be positive for this core geometry; the test is that it does not approach prompt criticality. Report the computed value with uncertainty for each core configuration.

**T3.4 — Delayed neutron fraction β_eff**
β_eff must be computed for each transition phase and must remain > 0.003 (consistent with Pu-heavy fast spectrum cores). A value below this threshold would call into question the feasibility of normal control rod operation.

**T3.5 — Prompt neutron lifetime Λ**
Λ must be computed and must fall in the range 3–6 × 10⁻⁷ s, consistent with fast reactor benchmarks. Values outside this range indicate a geometry or cross-section library error.

---

## T4 — Isotopic Evolution

**T4.1 — Pu vector at EG23 equilibrium**
At EG23 equilibrium (UNF recycle, cycles ≥ ~8), the discharge Pu vector must show ²³⁹Pu fraction stabilization. Compare against the Hoang et al. (2022) S-PRISM B&B upper-bound reference (74.5% ²³⁹Pu, 72.9% ²⁴¹Pu) noting that PRISM's lower burnup (87.5 vs. 290 MWd/kg) should yield a higher ²³⁹Pu fraction. Any result showing ²³⁹Pu fraction below the deep-burnup reference is flagged for review.

**T4.2 — Minor actinide buildup**
MA-to-TRU ratio (wt%) must be tracked each cycle and must reach a stable equilibrium value in EG24 mode. The equilibrium MA/TRU ratio must be lower in EG24 than in EG23 (since fast-spectrum fission of Am and Cm is the mechanism of MA destruction — the defining benefit of the EG24 cycle).

**T4.3 — Mass balance closure**
At every cycle, the total actinide mass in: (fresh fuel loaded) + (in-vessel inventory) + (cooling assemblies) + (NFRC input) must balance the previous cycle inventory corrected for fission and neutron capture to within 0.1% by heavy metal mass. This is both a physics and a code logic check.

**T4.4 — Fission product inventory**
Dominant neutron-absorbing fission products (Xe-135, Sm-149, Nd-148 as burnup monitor) must show physically plausible cycle-averaged concentrations. Xe-135 and Sm-149 must reach saturation equilibrium within the first depletion timestep of each cycle.

---

## T5 — Numerical and Statistical Convergence

**T5.1 — k-eff statistical uncertainty**
Each transport calculation must achieve k-eff 1σ uncertainty ≤ 15 pcm, as specified in the paper's methodology. Runs failing this threshold must be rerun with higher particle count before results are used.

**T5.2 — Reaction rate uncertainty**
Tally-based one-group reaction rates (fission, capture) for all actinide isotopes tracked in the depletion must have 1σ relative uncertainty < 0.5%. Isotopes failing this threshold in the dominant actinide set (U-238, Pu-239, Pu-240, Pu-241, Am-241) are flagged for increased tally resolution.

**T5.3 — Timestep sensitivity**
The 30-day depletion timestep must be validated by running one full equilibrium cycle at 15-day resolution. The resulting change in discharge burnup, CR, and Pu vector must be < 1% (relative). If larger, the 30-day timestep is not adequate and must be refined.

**T5.4 — Particle count sensitivity**
One benchmark case (UNF recycle equilibrium cycle) must be run at 2×, 5×, and 10× the baseline particle count (5 × 10⁶). k-eff must be stable to within 2σ across all runs, confirming that the baseline count is not source-converged to a local mode.

---

## T6 — Code Verification and Benchmarking (Appendix B)

**T6.1 — EBR-II benchmark**
Reproduce the published EBR-II Run 138B k-eff to within ±300 pcm using the ENDF/B-VIII.0 library. This is the standard fast reactor criticality benchmark for OpenMC validation in the sodium-cooled metallic fuel context.

**T6.2 — OECD/NEA SFR benchmark**
Reproduce the OECD/NEA sodium fast reactor criticality benchmark (e.g., the Phenix EOL benchmark or the BN-600 LEU core) to within ±500 pcm. Demonstrates that the transport geometry and cross-section treatment are correct for fast-spectrum hexagonal lattice cores.

**T6.3 — Independent depletion check**
For the first three burnup steps of the UNF recycle equilibrium cycle, reproduce isotopic concentrations using an independent depletion code (e.g., ORIGEN or Serpent) to within ±5% for all actinides contributing > 1% of total TRU mass. Flags any systematic error in the Bateman solver or one-group cross-section collapsing.

**T6.4 — Conversion ratio analytical check**
For the breakeven configuration (CR = 1.06), confirm analytically that the ratio of U-238 capture rate to total fissile destruction rate equals the computed CR to within 2%. This provides an independent check that the tally normalization and material update logic are self-consistent.

---

## T7 — Transition Logic

**T7.1 — EG23 equilibrium establishment**
The simulation must demonstrate that a stable EG23 equilibrium exists (k-eff, CR, and discharge isotopics all cycle-to-cycle stable) before the EG24 transition is initiated. "Stable" is defined as < 100 pcm variation in BOEC k-eff and < 1% relative variation in CR over three consecutive cycles.

**T7.2 — Monotonic TRU inventory evolution**
During the EG23→EG24 transition, the total TRU inventory in the active core must evolve monotonically in the direction predicted by the changing CR (decreasing for CR < 1, stable for CR ≈ 1, increasing for CR > 1). A non-monotonic result indicates a feedstock composition error or material transfer bug.

**T7.3 — EG24 equilibrium establishment**
The simulation must demonstrate that a stable EG24 equilibrium is reached, applying the same criteria as T7.1. Report the number of transition cycles required. This is a primary result of the paper.

**T7.4 — No inadvertent criticality excursion**
At no point during the transition may BOEC k-eff exceed 1.010 (allowing for the ±200 pcm control band). An excursion above this threshold triggers a halt and review of the TRU enrichment adjustment algorithm.
