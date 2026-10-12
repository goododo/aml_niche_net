# GATES.md -- one entry per gate: date, numbers, pass or fail

<!-- W1.1 start -->
## GATE_W1.1 positive control -- 2026-10-12

Script: scripts/ws/20_t3c_positive_control.R | config: scripts/ws/config_ws.R (frozen 2026-10-12)

| row | ligand | sender -> receiver | receptor | arm | n | z | p_perm | fake P95 | recovered | diag rho |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | IFNG | T_NK -> Mono_DC | IFNGR2 | main | 28 | -0.55 | 0.5534 | 3.55 | FALSE | 0.38 |
| 1 | IFNG | T_NK -> Mono_DC | IFNGR2 | sens | 28 | 0.10 | 0.4316 | 3.56 | FALSE | 0.38 |
| 2 | IFNG | T_NK -> LMPP_GMP | IFNGR1 | main | 17 | 1.14 | 0.1588 | 1.71 | FALSE | -0.08 |
| 2 | IFNG | T_NK -> LMPP_GMP | IFNGR1 | sens | 17 | 1.64 | 0.0899 | 1.61 | FALSE | -0.04 |
| 3 | TNF | Mono_DC -> LMPP_GMP | TNFRSF1A | main | 18 | 0.35 | 0.4216 | 2.45 | FALSE | 0.53 |
| 3 | TNF | Mono_DC -> LMPP_GMP | TNFRSF1A | sens | 18 | -1.09 | 0.6743 | 2.99 | FALSE | 0.52 |
| 4 | TNF | B_Plasma -> Mono_DC | TNFRSF1B | main | 28 | 0.93 | 0.3886 | 3.90 | FALSE | -0.35 |
| 4 | TNF | B_Plasma -> Mono_DC | TNFRSF1B | sens | 27 | 0.55 | 0.4565 | 4.98 | FALSE | -0.35 |
| 5 | TGFB1 | T_NK -> Mono_DC | TGFBR2 | main | 16 | -1.37 | 0.8871 | 2.52 | FALSE | 0.51 |
| 5 | TGFB1 | T_NK -> Mono_DC | TGFBR2 | sens | 15 | -1.43 | 0.9451 | 1.77 | FALSE | 0.54 |
| 6 | TGFB1 | T_NK -> LMPP_GMP | TGFBR2 | main | 10 | -0.10 | 0.4885 | 1.08 | FALSE | -0.19 |
| 6 | TGFB1 | T_NK -> LMPP_GMP | TGFBR2 | sens | 10 | 0.15 | 0.4386 | 0.90 | FALSE | -0.19 |

Pairs recovered -- IFNG: FALSE | TNF: FALSE | TGFB1: FALSE | arms agree: TRUE
**Verdict: FAIL** -- anything else -> check the pipeline, do not proceed
<!-- W1.1 end -->
