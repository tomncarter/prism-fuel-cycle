# sources/

Reference documents behind the model's design data. **The PDFs are not committed**
(`sources/*.pdf` is git-ignored): their redistribution terms are either unstated
(DOE contractor report) or restrictive (journal paper, book). Download your own copy
from the links below; the SHA-256 lets you confirm it is the same file.

| File (local) | Citation | Where to get it | Used for |
|---|---|---|---|
| `134264.pdf` | C. Grandy & R. Seidensticker (eds.), *Advanced Burner Reactor 1000MWth Reference Concept*, ANL-AFCI-202 (ANL-ABR-4), Argonne National Laboratory, Sept. 2007 | https://publications.anl.gov/anlpubs/2017/04/134264.pdf (blocks scripted downloads; use a browser) | Assembly pitch 16.142 cm, gap 0.432 cm, HT9 duct 0.394 cm; pin designs `anl_abr_driver` (Table II.1-2), `anl_sprism_driver` and `anl_hcr_driver` (Table III.6-10), `anl_blanket` (Table III.6-11); S-PRISM performance comparison (Table III.6-12) |
| (not stored) | B. S. Triplett, E. P. Loewen & B. J. Dooies, "PRISM: A Competitive Small Modular Sodium-Cooled Reactor," *Nuclear Technology* 178(2), 186–200 (2012) | Journal (Taylor & Francis / ANS) | Core configurations: assembly counts (Table III), fuel-cycle parameters (Table II), UNF-recycle core map (Fig. 8) |
| (not stored) | C. E. Till & Y. I. Chang, *Plentiful Energy: The Story of the Integral Fast Reactor* (2011), ISBN 978-1466384606 | https://www.thesciencecouncil.com/pdfs/PlentifulEnergy.pdf | Pin designs `tc_d090`–`tc_d120` (Table 14-3, p. 311; the 1.0 cm row's fuel fraction prints as 0.415, used as 0.313); Fig. 14-5 |

SHA-256 of the local copies:

```
85d7ab51c59c06d29641170859b3dd27f0bfb31f08dba2c6722e1da0c98f0ed2  134264.pdf
```

Where each number enters the code: `materials/tier1_homog.py` (`PIN_DESIGNS`),
`geometry/params.py` (`AssemblyGeom`), and the `source` / `description` fields of
`core_loadouts/prism_*.json`.
