# =============================================================================
# config_nodestatus.R   --   is a missing node missing BIOLOGICALLY or TECHNICALLY?
#
# THE PROBLEM THIS SOLVES. Unbalanced FGW absorbs mass differences between graphs.
# If a sample has no T_NK node, FGW will happily treat that as a topology
# difference -- which is correct when the marrow really has no lymphocytes
# (blast-dominated AML) and catastrophic when the protocol removed them
# (CD3-depleted library). Both look identical in a count matrix. The distinction
# is not recoverable from the data; it lives in the protocol, i.e. in the curated
# `sorting` / `cell_prep` / `tissue` fields.
#
# THREE STATES, per (sample, node):
#   present            cells are there
#   absent_biological  the protocol WOULD have captured this node and it is not
#                      there -> a real feature of this marrow, feed it to FGW
#   NA_technical       the protocol CANNOT contain this node -> mask it. Never
#                      zero-fill: a zero is a measurement, and this is not one.
#
# A FOURTH SIGNAL, orthogonal to the three: `abundance_unreliable`. Chen2023 and
# GSE253355 sort five FACS fractions and re-pool them, so every node EXISTS but
# its mass is a pipetting decision, not a biological abundance. The node is not
# NA -- dropping it would throw away real cells -- but the node MARGINAL, which is
# exactly what FGW transports, is not comparable to a whole-marrow sample. Flagged
# separately so the marginal can be renormalised or down-weighted without losing
# the node.
#
# EVERY RULE CITES ITS EVIDENCE. A rule that masks a node is deleting data, so it
# has to be arguable from the protocol, not from a hunch about what "sorted"
# implies.
# =============================================================================

NODE_STATUS_STATES <- c("present", "absent_biological", "NA_technical")

# The 8 hierarchy bins: 7 CCC-graph nodes + Stromal (in_ccc_graph = FALSE, kept as
# an auxiliary bucket -- see ROLE_TABLE Chen2023/GSE253355).
NODE_ALL <- c("HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte",
              "T_NK", "B_Plasma", "Stromal")

## -- RULE TABLE ---------------------------------------------------------------
# Matched against the resolved curated fields in 00_curated_manifest.csv.
# `mask` lists the nodes this rule sets to NA_technical.
#
# field/value are matched exactly (case-insensitive); several rules can fire for
# one sample and the masks are unioned -- a PB sample from a CD34-sorted library
# loses both sets.
NODE_MASK_RULES <- data.table::data.table(
  rule_id = c("PB", "PROG_GATE", "LIN_DEPLETED", "CD34_SORT", "CD34_CD38_SORT"),
  field   = c("tissue", "sorting", "sorting", "sorting", "sorting"),
  value   = c("PB", "blast_enriched_progenitor_gate", "blast_enriched_lin_depleted",
              "CD34pos", "CD34posCD38neg"),
  regex   = c(FALSE, FALSE, FALSE, FALSE, FALSE),
  mask    = c(
    # PERIPHERAL BLOOD. Stroma is a tissue compartment that is not in blood at all.
    # Erythroid precursors and megakaryocytes do not circulate -- normoblasts and
    # megakaryocytes are retained in the marrow by definition of the marrow-blood
    # barrier. HSC_MPP and LMPP_GMP are deliberately NOT masked: circulating
    # blasts are the whole reason PB is sampled in AML, so their absence in a PB
    # sample is informative rather than structural.
    "Stromal,Erythroid,Megakaryocyte",
    # CD117+/CD34+ POSITIVE selection, no lineage dump (GSE185991, curated
    # sorting_detail: "two libraries, Live CD117+ and Live CD34+; NO lineage
    # dump"). A positive gate on progenitor markers cannot yield mature T, B,
    # monocytes, megakaryocytes or normoblasts, and stroma is CD45-/CD34+ only for
    # endothelium. What survives is the progenitor compartment.
    "T_NK,B_Plasma,Mono_DC,Erythroid,Megakaryocyte,Stromal",
    # Same positive gate PLUS an explicit lineage dump (CD3-/CD19-/CD235a-).
    # The curation carries this as a CRITICAL instruction on 34 GSE185991 rows:
    # "negative selection removes T_NK, B_Plasma and Erythroid BY DESIGN - mask
    # these nodes as NA downstream, never zero-fill". The remaining three follow
    # from the positive gate, as above.
    "T_NK,B_Plasma,Mono_DC,Erythroid,Megakaryocyte,Stromal",
    # GSE116256 CD34+ sorted libraries (BM5-34p and siblings).
    "T_NK,B_Plasma,Mono_DC,Erythroid,Megakaryocyte,Stromal",
    # CD34+CD38- is a strictly tighter gate than CD34+.
    "T_NK,B_Plasma,Mono_DC,Erythroid,Megakaryocyte,Stromal"),
  evidence = c(
    "marrow-blood barrier: normoblasts and megakaryocytes do not circulate; stroma is not a blood compartment",
    "GSE185991 curated sorting_detail; GEO GSM5628163 characteristics 'cell type: Live CD117+'",
    "GSE185991 curator_note CRITICAL instruction; curated sorting_detail 'CD3- CD19- CD235a-'",
    "GSE116256 curated sorting",
    "GSE116256 curated sorting"))

## -- STROMAL: PRESENCE REQUIRES POSITIVE EVIDENCE -----------------------------
# Stroma is stated as a mask rule of its own because the default is the opposite
# of every other node. Mesenchymal and endothelial cells are adherent, low
# frequency, and are lost by density-gradient separation; a bone marrow aspirate
# loaded on a droplet platform contains essentially none unless the protocol
# specifically retained or enriched them. So the honest default for Stromal is
# NA_technical, and PRESENCE is what needs evidence -- not absence.
#
# Writing it the other way round ("mask when cell_prep is mononuclear") looked
# equivalent and was not: it silently left Stromal `present` for the 8 samples
# whose cell_prep the curation leaves blank (GSE207356 x3, Petti2019 x5). Those
# are exactly the rows with the least evidence, and they were getting the most
# permissive treatment. Under-specified metadata must not buy a node.
#
# A sample keeps Stromal only if its protocol positively supports it:
STROMAL_EVIDENCE_SORTING   <- c("multi_fraction_FACS_recombined")   # a stromal fraction was sorted and re-pooled
STROMAL_EVIDENCE_CELL_PREP <- "^whole_marrow"                       # no density gradient was applied

## -- ABUNDANCE-UNRELIABLE (orthogonal to the three states) --------------------
# The node exists and its cells are real; only its MASS is an experimental design
# choice. Chen2023 and GSE253355 sort five FACS fractions per donor and re-pool
# them, so the ratio between nodes is pipetted. FGW transports node marginals, so
# comparing these against a whole-marrow sample compares protocols.
ABUNDANCE_UNRELIABLE_SORTING <- c("multi_fraction_FACS_recombined")

## -- PRESENCE THRESHOLD -------------------------------------------------------
# Below this a node is reported absent. Whether that absence is `absent_biological`
# or `NA_technical` is decided by the rules above, NOT by the count -- which is the
# entire point: a zero cannot tell you why it is zero.
# Tracks CCC_MIN_CELLS_PER_NODE in config_ccc.R.
NODE_PRESENT_MIN_CELLS <- 10L

## -- FGW marginal handling ----------------------------------------------------
# How 07_fgw consumes the three states. Recorded here so the semantics live beside
# the rules that produce them rather than inside the solver.
#   present            node enters with its observed marginal
#   absent_biological  node enters with marginal 0 -- a real, transportable
#                      difference that unbalanced FGW is designed to absorb
#   NA_technical       node is REMOVED from this sample's graph and its marginal
#                      mass is redistributed over the surviving nodes, so the
#                      comparison is made on the sub-graph both samples actually
#                      observed. It is never entered as 0.
FGW_NA_POLICY <- "drop_and_renormalise"
