# Stage 0a: do the P01 hashtag doublets behave as random doublets? (plan v2, 2026-11-11)

Status: FINAL PLAN after review. Waiting for the user's OK on the file split; nothing has been coded or run yet.
- Supersedes STAGE0A_plan_draft.md, which is kept for the audit trail.
- Review: workflow wf_fc1f5a63-f6a, with reviewers for statistics, pre-registration/masking, code feasibility and public sources, plus a consolidation step.
  - All four reviews returned "fixable"; none found anything fatal.
  - The main session re-ran the consolidation's decisive simulations. The same seeds gave identical numbers.
- The authoritative frozen text will be the header of scripts/dbl/03b_hto_desk.R at git tag prereg-dbl-0a-v1. This folder is not under git.

## 1. Purpose (unchanged)
The share redesign (REDESIGN_share_plan.md) needs a real-doublet check, K1. K1 uses the cross-donor hashtag doublets of the P01 pool. These formed after pooling, so they are random by construction.

Stage 0a decides, from untyped hashtag data only, whether K1 is usable. If it is not, the line closes with BN-1. Nothing in 0a types a doublet.

## 2. What changed from the draft (and from the plain-language rules the user approved)
| Item | Draft | v2 | Reason |
|---|---|---|---|
| G1 | chi2/df > 2 fails | unchanged | |
| G2 numerator | all in-QC cross_donor calls | in-QC calls with exactly two tags, of two donors | Stricter. It matches G1/G3, and triplets no longer inflate the median |
| G2 bar | median ratio < 1.4 fails | unchanged | The contamination estimate f_hat is reported, not gating (see 13.1) |
| G3 | q75 doublet-like cut, band [0.67, 1.5] | unchanged; pi is exact, not 20,000 draws | q75 is robust if real doublets carry less RNA than the summed singlets; q90 is not (section 5) |
| Lanes | three lanes pooled | unchanged; count-1 alone and count-2+3 alone are reported | Pooling is what the user approved on 11-11. This amends redesign section 7, which said "count-1" |
| Known-answer worlds | 20 replicates, bars 0.8-0.9, sub-pool pair unnamed | 200 replicates, every bar >= 0.95, every world at a true rate >= 0.99 | 20 replicates fail correct code 7-21 % of the time |
| Failure handling | "read only if all pass" | fixed run order and failure table (section 3) | Prevents rerun-shopping and lane arguments made after the fact |
| Close text | one BN-1 text | one clause per gate (section 7) | "The cross-donor check failed" would be false for a G3-only failure |
| Code | 03b ~180 lines + 01 +6 + config +15 | two steps: config +4 / 01 +8, then config +14 / 03b ~280-300 lines | Validity machinery, provenance guard and write safety |

## 3. Run order and what each failure leads to
| Step | Content | On failure |
|---|---|---|
| (0) | Unit tests of decide() on 7 cases: all pass; G1 failing alone; G2 failing alone; G3 > 1.5 alone; G3 < 0.67 alone; one excluded lane (df 24); one failed world | Fix before the freeze (pre-freeze bug) |
| (a) | (1) md5 of config_dbl.R and config_paths.R equal the literals in the 03b header. Checked before any data is read. | Stop. Fix and resubmit; this is not a rerun |
| | (2) Count-1 regression (section 4.6). (3) dbl_checks_01_p01.csv all PASS. (4) Per P01 donor over the 3 lanes, B_union equals the QC report's n_doublet (51/196/925/71/581) and demuxed barcodes equal n_raw (2,186/4,597/10,054/2,499/8,350) | Rerun rule (below) |
| (b) | Per lane: (1) stage 03's four known answers: own-top share >= 0.80; class-A single calls name the author's donor >= 0.95; class-A cross_donor share < 0.05; in-QC class C cross share > class A. (2) >= 50 class-A single calls per tag. (3) Every tag in >= 1 in-QC exactly-two-tag cross call | A count-2/3 lane is excluded from the worlds and the gates (df 12 per remaining lane); this is not a close. If both are excluded, the gates run on count-1 alone, redesign section 7's original definition. Count-1 failing counts as (a) |
| (c) | The five gating worlds (section 5) on the remaining lanes | Rerun rule |
| (d) | Real G1-G3 computed, printed and written. A non-converged GLM, or a fitted-vs-observed tag-margin difference > 1e-6 x the lane total, counts as (c) | Verdict (section 7) |

Rerun rule:
- One rerun is allowed, only for a demonstrated code bug: a minimal reproduction showing the code departs from the header.
- It uses the same seeds and gets a dated amendment. Earlier outputs are moved to tables/superseded_<jobid>/, not deleted.
- Without a demonstrated bug the outcome is BN-1: "the cross-donor check could not be validated".

Until (c) passes, 03b writes only dbl_checks_03b_p01.csv and dbl_hto_desk_known_p01.csv: booleans, thresholds and synthetic results, never a real tag-pair count.

## 4. Sets and gates (literal)
4.1 **Positivity.** Stage 03's rule, copied verbatim (03_call_hashtags.R lines 35-55).
- A tag is positive iff count > max(quantile(bg, 0.995, type 7), 5).
- bg is the tag's counts on author-assigned (demuxed) barcodes of other donors, in any class, in or out of line QC.
- The AML0160 tags therefore need >= 6 UMI.

4.2 **Sets.**
- Base: in_counts == TRUE and exactly two positive tags (n_pos == 2).
- Cross: the two tags name different donors (19 tag pairs per lane).
- Two-tag: {AML0160.1, AML0160.2} or {AML3762.1, AML3762.2}.
- Barcodes with >= 3 tags enter no gate.
- nCount = bc$nc (RNA). Assert no NA nc among in-QC barcodes.

4.3 **G1.**
- Per lane, glm.fit(Z, y, family = poisson()) with default control.
- Z is the 19 x 7 tag-incidence matrix (rank 7); y is the 19 cross counts, zeros included.
- X2 = sum((y - mu)^2 / mu). G1 = sum of X2 over lanes / (12 x lanes).
- Fail iff G1 > 2.

4.4 **G2.** Median nc of in-QC exactly-two-tag cross calls, pooled over lanes, divided by median nc of class-A single calls, pooled over lanes. Fail iff G2 < 1.4.

4.5 **G3.**
1. t_l = quantile(nc of the lane's class-A single calls, 0.75, type 7). A barcode is doublet-like iff nc >= t_l. The cut is lane-pooled by design: AML3762 singles pass at 0.37, AML0160 at 0.22.
2. pi_jk = the exact share of all pairs (a of tag j, b of tag k) of the lane's class-A single calls with nc_a + nc_b >= t_l. Computed with sorted findInterval; no Monte Carlo.
3. Refit G1's model on the doublet-like cross counts with offset log(pi_jk).
4. Per split donor d: P_l,d = exp(eta_d.1 + eta_d.2) x pi_d.1,d.2, and O_l,d = the doublet-like in-QC two-tag count.
5. G3 = sum O / sum P over the passing lanes and both split donors, unrounded. Fail iff G3 < 0.67 or G3 > 1.5; equality passes.

4.6 **Count-1 regression** (calls recomputed in memory against the stored 11-09 files; only booleans are printed).
- Thresholds: all.equal to dbl_hto_map_smoke.csv.
- Cell order, call, donors_pos, class, donor and in_counts: identical() to the stored calls file.
- n_pos and hto_total: identical() after as.integer().
- The class x in_counts x call counts: all.equal to dbl_hto_summary_smoke.csv.
- Calls always use the in-memory thresholds.

## 5. Known-answer worlds (200 replicates each; every gating bar >= 0.95)
**Pseudo-lanes.**
- One pseudo-lane per passing real lane, made of N_l synthetic doublets. Each doublet is two class-A single calls of that lane, drawn uniformly with replacement. Same-tag and same-donor pairs stay in the draw.
- HTO uses the "one-ambient" construction: H = H_a, then H[t_b] += max(H_b[t_b] - m(t_b), 0). Here m(t) is the mean count of tag t on the lane's class-A singles of other donors. Background is droplet-level, so the plain sum creates 15 % false third tags.
- nCount = nc_a + nc_b. Calls use the real lane's thresholds, held fixed.
- N_l = round(n_l / r_l). n_l is the lane's in-QC exactly-two-tag cross count. r_l is the share of 100,000 random pairs called exactly-two-tag cross; for count-1, r = 0.749, so N_l is about 3,900.
- The reference (t_l, pi, G2's denominator) is the lane's class-A singles, unmodified.

| World | Construction | Gating bar | Measured on count-1 (3 pseudo-lanes) |
|---|---|---|---|
| W1 random | pseudo-lanes as above | P(all three gates pass) >= 0.95 and abs(mean log G3) <= 0.05 | all pass 1.000; mean log G3 -0.009 |
| W2 excess x2 | W1 plus Poisson(n_d) extra (d.1, d.2) doublets per split donor; n_d = the replicate's count of (d.1, d.2) pairs among W1's draws | P(G3 > 1.5) >= 0.95 | 1.000 |
| W3 borrowed x1.75 | W1 plus Poisson(0.75 n_d) single calls of donor d, drawn uniformly over both aliquots; sibling tag set to floor(threshold) + 1; own nc | P(raw O/P >= 1.5) >= 0.95 and P(0.67 <= G3 <= 1.5) >= 0.95 | raw 0.990; G3 in band 1.000 (mean 1.286) |
| W4 sub-pool x2 | W1 plus Poisson(n_AB) extra AML0160-AML3133 doublets; n_AB = that pair's count among W1's draws (the smallest donor pair) | P(G1 > 2) >= 0.95 | 1.000 (every pair is 1.000 at x2) |
| W5 contamination | W1; each in-QC exactly-two-tag cross call's nc is replaced with probability 0.7 by a random class-A single nc. G2 only | P(G2 < 1.4) >= 0.95 | 1.000 (median 1.236) |

Reported only, never deciding:
- excess x1.5 and x1.75 (G3 detection 0.50 and 0.985);
- the other 29 donor-pair x {1.3, 1.5, 2} sub-pool rows (G1 power at x1.5: 0.625-0.995; at x1.3: 0.105-0.455);
- borrowed singlets drawn in proportion to nc;
- W1 and W2 with the plain HTO sum.

Seeds: DBL_SEED_0A = SEED + 961.
- set.seed(DBL_SEED_0A) before the sizing draw.
- + i before gating world Wi.
- + 10 + k before the k-th reported world, in table order.
- + 99 before the bootstrap.
- Replicates run sequentially within a world. The range 961-1060 is free; existing streams use 801-952 plus the +1000/+5000/+7000 ranges.

## 6. Reported, never deciding
- Per lane and per split donor: O, P and O/P, raw (G1 fit) and doublet-like.
- G3 restricted to count-1, and to count-2+3 pooled (a blind replicate of the exposed lane), each with its CI.
- G3's 95 % CI from 2,000 draws. Each lane's 21 doublet-like counts are redrawn as Poisson(observed) and refitted with the same offsets. Type-7 2.5/97.5 % percentiles; non-converged refits are counted and dropped. The reference is not resampled.
- G3 at q90, and with pi computed at kappa = 0.8 (doublet nCount = 0.8 x sum).
- G3 with pi from class-A plus in-QC B_mad single calls (AML3762 is 27.8 % B_mad).
- G2 on all in-QC cross_donor calls; the >= 3-tag count.
- kappa_hat and f_hat:
  - kappa_hat = q75 of in-QC exactly-two-tag cross nc / q75 of the exact mixture of nc_a + nc_b over the observed cross pairs.
  - f_hat = sum_l n_l (F_x,l - F_s,l) / sum_l n_l (0.25 - F_s,l). F_x,l is the share of the lane's cross calls at or below the lane's single-call q25. F_s,l is the mixture CDF at q25 / clip(kappa_hat, 0.6, 1).
- Class (A / B_mad / B_union / C) of two-tag calls; nCount medians by call among in-QC barcodes.
- r_low deciles for two-tag vs cross calls, each split by doublet-like. r_low is the minimum over a call's two tags of count / the median own-tag count of that tag's class-A singles. A real doublet has both tags near 1; a borrowed tag gives one near threshold.
- Pooled standardized Pearson residual per cross pair, with the AML3762 pairs listed next to G3.
- Leave-one-donor-pair-out G3: 10 refits, each of rank 7, as a min-max range. A pass whose range crosses 0.67 or 1.5 is reported as "pass, model-sensitive". The verdict stays "K1 usable".
- Restricted-fit X2/df; exp(eta_j) / sum exp(eta) beside each tag's class-A single share.
- Bootstrap failure count.

## 7. Pre-written verdict clauses
- All gates pass: "K1 usable"; go to plan 0b.
- G1 fails: BN-1, "cross-donor doublets do not follow random pairing after pooling (G1 = x)".
- G2 fails: BN-1, "cross-donor calls do not carry doublet-level RNA (median ratio x; kappa_hat y; f_hat z)". If f_hat <= 0.25, the reading is that real doublets carry less RNA than summed synthetic doublets, so the summed-synthetic calibration K1 checks would not transfer.
- G3 fails (BN-1b): "same-donor doublets of the P01 pool depart from random formation after the sibling aliquots met (doublet-like two-tag O/P = x [95 % CI]); aggregation or a tag artefact cannot be separated from pairing, so a share result could not be read as pairing".
- Validity failure without a demonstrated bug: BN-1, "the cross-donor check could not be validated".

The typed two-tag share stays a reported negative control in redesign Stage 4, whatever G3 gives. These consequences hold whatever .1/.2 later turn out to mark.

## 8. What a pass or a failure can and cannot mean (header statements)
1. The gates pool the three P01 lanes, as the user approved on 11-11; count-1 alone is reported.
2. What .1/.2 mark is not stated in any public record: GEO, BioSample, the paper and its supplement, the authors' code, or SCP1987.
   - The authors map both tags to one patient. Supplementary Table 1 lists one diagnosis bank sample per donor.
   - Most likely they are two separately hashed tubes of one diagnostic specimen that met at pooling. AML3762.2 stains about 1.9x brighter than AML3762.1, so the tubes were stained separately.
   - G3 therefore tests same-donor association that formed after the sibling tubes met. A pass does not exclude aggregation before pooling: within-tube aggregates carry one tag.
3. G1's 7 free tag effects absorb loading, hashtag sensitivity and generic per-donor stickiness, so no gate tests generic stickiness.
4. The class-A reference was already cleaned by the authors. They removed HTODemux multi-tag barcodes, including .1+.2, and their recoverDoublets neighbours.
5. The verdict holds only for this rule set: the 0.995-quantile / 5-UMI rule, per-lane thresholds and name-based tag mapping. 03b replaces redesign section 9's planned "03 --desk". 03b writes an md5 of paste(cell, call) per lane, so a later generalised 03 can prove it makes the same calls.

## 9. Masking and write safety
**Masking.**
- Right after each readRDS, 03b keeps only x$hto (reordered to bc$cell) and x$bc[, .(cell, donor, demuxed, class, in_counts, nc)], then removes x.
- It reads no labels file and no author metadata. The ann element of the count-2/3 lane objects is not read before freeze 1.
- No output carries a cell type or a composition by tag or call.

**Changes to 01.**
- KEYS = c(DBL_SMOKE_LANES, DBL_P01_LANES[c("p01_2", "p01_3")]).
- Stop on an NA key (p01_1 is never a key, so the count-1 smoke object is never re-extracted), on mixing smoke and P01 keys, and on any P01 request other than exactly p01_2,p01_3.
- The table suffix is _smoke or _p01, following the key set.
- The class and check tables are written only when every requested lane was extracted in this run. This also fixes the current code, which shrinks dbl_checks_01_smoke.csv when lanes are skipped.

**config_dbl.R.** DBL_P01_LANES, plus stopifnot(identical(DBL_P01_LANES[["p01_1"]], DBL_SMOKE_LANES[["pooled"]])). The P01 lanes are never appended to DBL_SMOKE_LANES, which eleven stages loop over.

**03b.**
- Copies none of 03's writes. Stage 04 reads the count-1 calls file, and the _smoke map/summary tables are 0a's own regression reference.
- Writes no per-barcode file and nothing named *_smoke.
- Writes only dbl_checks_03b_p01.csv, dbl_hto_desk_known_p01.csv, dbl_hto_desk_map_p01.csv (per-lane thresholds, own-top shares, the calls md5), dbl_hto_desk_pairs_p01.csv and dbl_hto_desk_gates_p01.csv.
- Stops if any of these exists and --force is not given. With --force, it first moves them to tables/superseded_<jobid>/.
- Asserts that every output path starts with DBL_SMOKE_DIR. results/tables/*.csv is not git-ignored, and the 23:30 cron pushes to a public remote.

## 10. Provenance
- The 0a constants live in config_dbl.R section 11a.
- The 03b header holds the rule text and the md5 literals of config_dbl.R and config_paths.R. 03b refuses to run on a mismatch, and records both md5s and the md5 of its dated copy in dbl_checks_03b_p01.csv.
- Commits are by explicit path only: git add -- <new file>; git commit -m ... -- <paths>. Never -A or -a. Nothing under scripts/ws or results/tables/ws is touched (another workstream commits to main today).
- Both commits land before 23:30. If 03b is not ready by then, it leaves the tree for the night.
- Tag prereg-dbl-0a-v1 on commit 2. SMOKE_REPORT section 6h records the tag, the commit, the blob hashes and the md5 of the dated copy that ran.
- config_dbl.R is not edited between the tag and the start of the 03b job. The md5 guard enforces this.

## 11. Code and runs (two steps, 03_call_hashtags.R unchanged)
**Step 1, commit 1.**
- config section 11a: DBL_P01_LANES and stopifnot (+4 lines). 01: +8 net lines.
- Run `sbatch scripts/dbl/submit_dbl.sbatch 01 --lanes=p01_2,p01_3`. Expected about 3 min and 18-19 GB (job 3687229: 2:44, 18.1 GB).
- Read dbl_checks_01_p01.csv. A failure here is fixed before 03b exists, so it is not a rerun.

**Step 2, commit 2 + tag.**
- config section 11a, +14 lines: DBL_0A_G1_MAX 2; DBL_0A_G2_MIN 1.4; DBL_0A_DL_Q 0.75; DBL_0A_G3_BAND c(0.67, 1.5); DBL_0A_MIN_SINGLES 50L; DBL_0A_N_SIZE 1e5; DBL_0A_N_REP 200L; DBL_0A_BAR 0.95; DBL_0A_BIAS_MAX 0.05; DBL_0A_MULT c(excess = 2, borrow = 1.75, subpool = 2); DBL_0A_SUBPOOL c("AML0160", "AML3133"); DBL_0A_CONTAM 0.7; DBL_0A_N_BOOT 2000L; DBL_SEED_0A.
- New scripts/dbl/03b_hto_desk.R, about 280-300 lines, a third of it frozen header text.
- Self-test before the tag: `03b --selftest`.
  - It runs (0) and W1-W5 on count-1 class-A single calls only: 3 pseudo-lanes, N = 3,900, the sizing the reviews used.
  - It prints the world table and writes nothing to DBL_TAB_DIR.
  - It reads no real pair count, so it is a pre-freeze bug check, not a run.
- Commit, tag, then run `sbatch scripts/dbl/submit_dbl.sbatch 03b` (under 5 min).

There is no "results commit": outputs stay on /LARGE1, and SMOKE_REPORT and memory are outside the repo.

After the run:
- SMOKE_REPORT section 6h and memory are updated.
- Pass: plan 0b, which re-scores the smoke synthetics and the planning power.
- Fail: the BN-1 close-out in SMOKE_REPORT section 11.

## 12. Exposures before the freeze (copied into the 03b header)
- P01 count-1 hashtag-level tables of 11-09 (SMOKE_REPORT section 4):
  - thresholds 5, 5, 32, 29.275, 124.275, 14, 46; own-top shares 0.956-1.000;
  - class x in_counts x call: all classes cross 5,874, two-tag 602, negative 1,275;
  - in-QC: cross 2,945 = A 115 + B_mad 16 + B_union 16 + C 2,798; two-tag 244 = A 8 + B_mad 4 + B_union 1 + C 231.
- Ratios derived from those tables:
  - 244/2,945 = 0.083 vs <= 0.047 under a 50/50 split with class-A loads (>= 1.75x);
  - 244 vs about 109 (2.2x, redesign designs);
  - 1.32x with author-assigned A + B_mad + B_union loads;
  - purity-corrected 1.27-1.71x;
  - 0.122 (358/2,929) among out-of-QC barcodes.
- scDblFinder called fractions of count-1's hashtag sets, main/devel: cross 0.43/0.39, two-tag 0.25/0.24, class-A single 0.0145/0.0047.
- Author-assigned donor counts for count-1..3 (recon tables, qc_percell, QC report):
  - n_raw / n_doublet per donor;
  - author-kept donor shares per lane: AML0612 0.36-0.37, AML3762 0.29-0.31, AML0310 0.15-0.18, AML3133 0.08-0.10, AML0160 0.07-0.09.
- Count-2/3 barcode classes (the code reviewer's dry run of 01's class logic):
  - count-2: A 7,722 / B_mad 1,154 / B_union 554 / C 17,763;
  - count-3: A 8,143 / B_mad 981 / B_union 446 / C 16,794;
  - 04_annotated cells per donor.
  - No hashtag values.
- Structure only: metadata columns, HTO feature names and sequences, matrix dimensions, called-barcode counts.
- Review-round computations on count-1 class-A single calls and synthetic doublets only:
  - loads 357/281/704/771/821/2,374/1,105 (0160.1, 0160.2, 3133, 3762.1, 3762.2, 0612, 0310);
  - nCount q25/q50/q75/q90 3,599/5,236/7,011/8,966;
  - per-tag background and sibling-tag background;
  - own-tag medians;
  - pi at q75 (cross 0.79-0.95, two-tag 0.965/0.899);
  - every operating characteristic in section 5;
  - the G2 x kappa table.
- The 10-28 hashtag tally covers GSM5613747-52 only, not P01.

No session or agent has computed any of the following on any P01 lane beyond the count-1 class x call table: hashtag calls, tag-pair matrices, two-tag or >= 3-tag counts, or the nCount or HTO profiles of cross or two-tag barcodes. No hashtag call exists for count-2 or count-3.

## 13. Residual risks (stated, not fixed)
1. **G2 and low-RNA doublets.** G2 closes the line if real doublets carry <= about 70 % of the summed RNA of their two cells, even with no contamination. Simulated G2 is 2.04 at kappa 1, 1.63 at 0.8 and 1.43 at 0.7. The rule is user-approved; the G2 clause reports kappa_hat and f_hat, so the reading is explicit.
2. **Borrowed tags that track RNA size.** q75 would false-fail about 30 % of the time if borrowed tags tracked nCount. The data do not support that mechanism: own-tag HTO vs nCount Spearman 0.05; sibling background is flat across nCount quartiles and equals other donors' background.
3. **Pair stickiness that G1 misses.** G1's power at x1.3 is only 0.105-0.455. Such stickiness moves G3, and AML3762 pairs pull it down (x1.5: 0.74-0.89). The residuals and the leave-one-pair-out range are reported for this.
4. **Singlets with a stray second tag.** The worlds hold doublets only. Real singlets with a stray second tag above threshold, roughly 6-15 % of cross calls, are not modelled. They enter both the cross and the two-tag sets, and the doublet-like cut removes most of them. The class and r_low profiles show their size.
5. **Aggregation before pooling.** A G3 pass does not exclude it (section 8.2).

## 14. Pre-freeze amendments (2026-11-12, before the tag; no real gate value had been computed)
Implementation and tests on 2026-11-11 and 11-12:
- step 1: commit 76bcd37; 01 on count-2/3 gave 34/34 PASS;
- --selftest: 33 PASS;
- --dryrun=null: K1 usable, G3 1.002;
- --dryrun=excess: BN-1b, G3 2.016;
- the review workflow wf_edee2922-f6a was interrupted by a session end, and its interim notes were mined.

Changes found necessary and applied before the freeze:
1. **The worlds always use 3 pseudo-lanes.** A pseudo-cohort now has 3 pseudo-lanes, cycling over the passing lanes. Before, it had one pseudo-lane per passing lane.
   - With count-2 and count-3 both excluded (count-1 alone), W3's bars sat at 0.925 and 0.93 from noise alone. The pre-registered fallback "count-1 alone" would therefore always have ended INVALID.
   - With 3 passing lanes the worlds are unchanged.
2. **A zero-row grouping crash is fixed.** A lane with zero in-QC two-tag calls crashed the reported block, after the gates were computed but before anything was written. Reproduced on synthetic input.
3. **The verdict is written before the reported items.** The checks table, the known table and the core gates with the verdict are now written to disk before the reported items are computed. A crash in reporting can no longer lose the verdict.
4. **The blind-replicate label** names the lanes actually used.

Exposure, stated explicitly: the known table records the pseudo-lane sizes N_l = round(n_l / r_l). They reveal each lane's total in-QC exactly-two-tag cross count n_l. That total enters no gate.

Provenance note: the 23:30 cron of 11-11 auto-committed and pushed the pre-freeze draft as 06f9f87, together with the config section 11a constants. That draft's 03b md5 is 42a2db03. The frozen version is the tagged commit prereg-dbl-0a-v1, and config_dbl.R is unchanged since 06f9f87 (md5 d933a1b5).

Independent checks, kept in prefreeze_checks/:
- recompute_from_dryrun.R: G1, raw O/P, G3, restricted X2/df and the pooled residuals, recomputed with stats::glm, match 03b to 1e-14;
- pi_bruteforce.R: the exact pi equals brute force at kappa 1 and 0.8, and pair coding equals a naive loop;
- gen_lanes.R + run03b_real.R: synthetic REAL-branch scenarios.
