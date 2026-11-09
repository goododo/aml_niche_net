# g07_sim: simulation and derivation scripts behind the A1-G07 rules

Decision 9.8 of `scripts/ws/PREREGISTRATION_A1_G07.md` (section 8.4). These are the scripts that set the G07 rules
(removal readings, the stacked T x X model and rho_X, the S5b thinning, floor and shape rule, the delay-list rule). They
are copies, kept for provenance; G07's results come from `43_a1_g07.R` and `44_a1_g07_cells.R`, which never run them.

- **Verbatim copies.** Each file is byte-identical to the one listed in 8.4, so its md5 matches the frozen list. The
  paths inside are the original absolute paths and were not made relative (dated note, section 11). To rerun one, read
  the original location in the table and the input paths at the top of the script.
- **Original location.** `/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/` plus the subfolder shown here. That folder
  is outside the repository.
- **Derived inputs** (QC-level tables, not committed). Copied read-only, with `MD5SUMS.txt`, to
  `/LARGE1/gr10634/gaozy/aml_niche_net/ws_a1_g07/sim_inputs/`:
  - `01_hallmark_coverage.csv` e9585821
  - `real_background_profile_nonF402.rds` 262491be
  - `real_T_percell_umi.csv.gz` 6d037ce4
  - `real_sample_summary.csv` 86656d46

| folder | scripts (md5, first 8) | what they set |
|---|---|---|
| `inventory/scratch_D1/` | 01_inventory_facts 335b46fe; 02_sim_geneset 94730713; 03_sim_crosscomp 30079a4d; 04_summarise a2159db9; 05_crosscomp_readout 76028d73 | camera on the S5 fit; the cross-compartment read (D1) |
| `inventory/scratch_D2/` | 01_real_structure fd59f762; 02_real_tables e3e2ba66; 10_engine 11573063; 11_validate 21da522c; 12_validate_ucell f7723e19; 13_sim bf631cdb; 14_summarise 28c46716; 16_rule_table 61e976a8; 17_final_rule e52077fd | the S5b engine, thinning and the shape rule (D2) |
| `review1/L1/` | 03_membership_and_design b789289d; 04_sim_k5_stacked 7c7b216c; 04b_summarise_k5 9d505d53; 05_sim_removal 09788f43 | K5 and the removal reading (review 1) |
| `review1/L2/` | t2_k6b_sim c83a43c5; t6_k6b_sim_p3c 24e7ab9f | K6b (review 1) |
| `review2/R2b/` | 02_savage_rule c9792664; 03_sim_removal_v2 02e44a9d; 04_k6b_v2 9a2a5259; 04b_gate_shape_bias 59f7f09b; 07_sim_rhoX f6b79e5a; 09_oetjen_axis 92f18f26 | the delay-list rule, removal v2, graded K6b, rho_X, Oetjen's gene axis (review 2) |
| `drafting/` | detection_floor 81cdc674; membership_counts.py 21602863; removal_counts.py 1df874e5; savage_celltypes 87bace41 | the detection floor; the counts quoted in 1.4 and 2.3 |
