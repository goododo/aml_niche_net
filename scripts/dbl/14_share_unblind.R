#!/usr/bin/env Rscript
# 14_share_unblind.R ----
# Doublet-pairing redesign (share design), Stage 4: unblinding. STUB written at freeze 1 (tag prereg-dbl-share-f1); the code is
# written and approved only after freeze 2.
# Rules (fixed now): scripts/dbl/prereg_share/PREREGISTRATION_SHARE_F1.md sections 4-6 and REDESIGN_share_plan_2026-11-11.md
#   section 7, Stage 4:
#   - refuses to run unless the freeze-2 record matches (git rev-parse HEAD:scripts/dbl, md5 of the parameter table and of this file);
#   - U0: singlet rates calibrated inside this worker (10 discard-mode runs per lane, a disjoint tenth of class A each); a lane
#     above 0.008 makes its patients descriptive and n >= 10 is rechecked;
#   - U1: K1 = b-hat-0 on DOUBLET-LIKE cross calls only (freeze-1 amendment 1); its 95 % CI entirely below 0.80 -> D (BN-3);
#   - U2: outcomes A / B-small / B-undetermined / C / D (redesign section 7), psi_T qualifier, the reported list (section 5);
#   - seeds DBL_SEED_SHARE_UNB.
stop("14_share_unblind.R is a freeze-1 stub: it runs only after freeze 2 and the approval of its code")
