# coralesce 1.0.2

* Added a standalone vignette (html file)

# coralesce 1.0.1

* Additive update. Two new exported functions plus documentation clarifications.

* New functions

`collapseToGenets(dataset, genetAssignment, representative, drop_unassigned)`
Reduces a genotype data frame to one representative colony per genet, using
genet assignments from `runGenets()`. Picks the most-complete colony per genet
(lowest `pctNull`) by default; keeps unassigned/inadequate colonies as
individuals unless `drop_unassigned = TRUE`. Matches rows on `Coral_ID`, so
identifiers may contain any characters. Pure row filter — genotype values are
untouched.

`computeKinship(data, subset, targetN)`
Data-frame entry point to the kinship pipeline — the in-memory counterpart to
`runKinship()` (no `Data`/`Results` folders). Runs the same internal steps
(translate to IUPAC codes, set aside all-NA colonies, omit invariant loci,
compute kinship) and returns the same `PopAvgMKGD` / `MK_init` / `MK_final`
list, keyed by `Coral_ID`. Verified to return results identical to
`runKinship()` on the same data. Errors informatively if `Coral_ID` is missing
or if the input is single-letter coded data rather than paired alleles.

* Together these enable clone-corrected diversity in one line

eligible-pool diversity (all colonies, ramets included): computeKinship(raw)
clone-corrected (one colony per genet): 
ga <- runGenets(PctMatchThreshold = 90, PctNotNullThreshold = 50)$genetAssignments[[1]]
computeKinship(collapseToGenets(raw, ga))

# coralesce 1.0.0

* initial release

