# Dimer Matching

Dimer matching extends an MS1 search by also looking for **crosslinked dimers** which is two muropeptide stems
joined by a transpeptidase-mediated bond; rather than just monomers. See the [Dimer Matching
section](usage.md#dimer-matching) of the Usage guide for how to use this from the WebUI or command line. This
page documents the chemistry rules and how to define your own for a species that isn't pre - built in.

## How donors and acceptors are chosen

- **Donors** are selected from the muropeptides actually *detected* in your MS1 data. Specifically, from
  the detected monomers that are structurally valid donors for the crosslink type being considered. Only a
  detected, intact monomer can plausibly have reacted as a donor. Among that structurally eligible pool,
  donors are taken from most to least abundant (by summed intensity) until their combined intensity covers
  a configurable **donor abundance threshold** (`donor_abundance_threshold`, default 90%) of the eligible
  pool's total intensity, not of all detected monomers, so an abundant but structurally ineligible peak
  never dilutes the budget or gets selected as a donor itself. A threshold near 0% recovers the single most abundant donor
  behavior as a special case (at least one donor is always included whenever any are eligible); a threshold
  of 100% includes every structurally eligible donor, regardless of abundance.
- **Acceptors** are selected from the *theoretical* mass database you supply (e.g. one of the [built-in
  mass databases](data_dictionary.md#target-structures), or your own custom one) Every structure in the
  database that qualifies as an acceptor is used, since an acceptor doesn't need to have been detected
  itself.

Which donors were actually selected — their structures, intensities, and their share of the eligible
pool's cumulative intensity  is recorded in a `donors_used` output (`*_donors_used.csv`). This file is
produced automatically by the CLI alongside results, by the WebUI's **Generate Theoretical Dimers** preview,
and by **Run Analysis** when dimer matching is enabled.

During crosslink formation, the donor loses one water molecule (forming the new bond) and, depending on the
crosslink type, one terminal D-Ala residue (a free leaving group, not consumed into the new bond). Both
masses are subtracted using the monoisotopic convention used throughout PGFinder (see the
[Residues](pglang.md#residues) table in the PGLang reference — `A` = Alanine = `89.047678`, matching the
constant used here up to rounding).

### Modified monomers

A monomer carrying a [modification](data_dictionary.md#modifications) (e.g. Anhydro-MurNAc, Amidation,
Deacetylation) is still recognised as a donor or acceptor. Its modification tag (e.g. `"gm-AEJAA (Anh)|1"`)
is set aside before the stem is checked against the donor/acceptor rules, and re-attached afterwards. The
modification's mass is preserved automatically, since donor/acceptor masses come straight from the
already-modified `Theo (Da)` value, a crosslink involving a modified stem is therefore both correctly
*recognized* and correctly *priced*. 

Modifications that rewrite the glycan prefix itself rather than appending a tag (`Extra Disaccharide (+gm)`,
`Lactyl Peptides (Lac)`, `Loss of Disaccharide (-gm)`) are not yet handled and are excluded from donor/acceptor
candidacy, same as before, these are rarer than the tag-style modifications above.

## Built-in species rules

| Species | Code | Type | Crosslinks | Donor | Acceptor | Lost residue |
|---|---|---|---|---|---|---|
| _Escherichia coli_ | `ecoli` | DAP-type (`J`) | 4-3, 3-3 | **4-3**: pentapeptide (length 5), `J` at position 3, ends `...AA`. **3-3**: tetrapeptide (length 4), `J` at position 3 | `J` at position 3, length ≥ 3 | D-Ala |
| _Fusobacterium_ | `fusobacterium` | DAP-type (`J`) | 4-3, 3-3 | same as _E. coli_ | same as _E. coli_ | D-Ala |
| _Clostridioides difficile_ | `cdiff` | DAP-type (`J`) | 4-3, 3-3 | same as _E. coli_ | same as _E. coli_ | D-Ala |
| _Bacillus subtilis_ | `bsubtilis` | DAP-type or Lys-type | 4-3, 3-3 | either `J`- or `K`-at-position-3, same length rules as above | either `J` or `K` at position 3, length ≥ 3 | D-Ala |
| _Staphylococcus aureus_ | `saureus` | Lys-type (`K`) + pentaglycine bridge | 4-3-bridge only | pentapeptide (length 5), `K` at position 3, ends `...AA` | `K` at position 3, length ≥ 3, pentaglycine bridge (exactly 5 glycines in strict mode, 1-5 in permissive mode) | D-Ala |
| _Enterococcus faecalis_ | `efaecalis` | Lys-type (`K`) + D-Asp bridge | 4-3-bridge, 3-3 | **4-3-bridge**: pentapeptide, `K` at position 3, ends `...AA`, no D-Asp bridge. **3-3**: tetrapeptide, `K` at position 3, no D-Asp bridge | **4-3-bridge**: `K` at position 3, has D-Asp bridge. **3-3**: length 4, `K` at position 3, no D-Asp bridge | D-Ala for 4-3-bridge; nothing for 3-3 |

`J`/`K` here refer to the position-3 stem residue symbols from the [PGLang Residues
table](pglang.md#residues) — `J` is meso-DAP, `K` is Lysine.

Only `ecoli` and `cdiff` currently have bundled mass databases under `lib/pgfinder/masses/` (see [Target
Structures](data_dictionary.md#target-structures)); the other species require you to supply your own
`--masses_file` / custom mass database.

## Donor trimming mechanism

Two transpeptidation mechanisms are modeled, both of which release one free D-Ala as a leaving group and
consume one water molecule forming the new bond and they differ only in which bond is cleaved:

- **4-3 (D,D-transpeptidase)**: the donor pentapeptide's D-Ala4-D-Ala5 bond is cleaved, releasing D-Ala5.
  The new bond forms between the donor's D-Ala4 carbonyl and the acceptor's position-3 side chain amine.
  The donor is trimmed from a pentapeptide (length 5) to a tetrapeptide (length 4) in the resulting dimer.
- **3-3 (L,D-transpeptidase)**: the donor tetrapeptide's mDAP3-D-Ala4 bond is cleaved, releasing D-Ala4.
  The new bond forms between the donor's mDAP3 carbonyl and the acceptor's position-3 side chain amine.
  The donor is trimmed from a tetrapeptide (length 4) to a tripeptide (length 3).

In both cases the **acceptor is unconstrained in length** as it only needs an accessible side-chain amine at
position 3, regardless of how many residues follow it.



## Bridge mechanisms

For some species, the crosslink doesn't bond the two stems directly as it passes through a peptide bridge
attached to the acceptor's side chain:

- **Pentaglycine bridge** (_S. aureus_): the donor's D-Ala4 bonds to the N-terminus of a chain of glycines,
  whose C-terminus is attached to the ε-amino group of the acceptor's Lys3 side chain. Bridge length is
  conventionally 5 glycines, written `[GGGGG]` in PGLang notation (see [PGLang
  Syntax](pglang.md#syntax)).
- **D-Asp/D-Asn bridge** (_E. faecalis_): a single D-aspartate (or D-asparagine, strain-dependent) residue
  bridges the same way, written `[D]`.

## Custom crosslink rules

If your species isn't one of the six built in, the WebUI lets you define a custom rule instead of using a
species code. This is **not currently available from the command line** — only via the WebUI.

### Donor / acceptor patterns

Both patterns are regular expressions matched against the bridge stripped stem sequence (e.g. `AEJAA`).
The residue letters from the [PGLang Residues table](pglang.md#residues), with any bridge bracket already
removed.

Anchoring matters:

- `^` anchors to the start of the stem, `$` anchors to the end. Use `$` whenever the rule depends on an
  *exact* stem length — for example, `^..JAA$` matches only a pentapeptide ending in mDAP-Ala-Ala (positions
  3, 4, 5), and nothing longer.
- Omit `$` when the rule only cares about a prefix — for example, `^..J` matches any stem with mDAP at
  position 3, regardless of how many residues follow. This is normally what you want for an **acceptor**
  pattern, since acceptors are unconstrained in length (see above).

### Donor loses terminal D-Ala

Toggle this on if your crosslink mechanism cleaves a terminal D-Ala from the donor (true for both the 4-3
and 3-3 mechanisms described above); leave it off if the donor's full matched structure is what ends up in
the dimer unchanged.

### Acceptor bridge

If your species crosslinks through a peptide bridge, choose the bridge type so the acceptor pattern alone
doesn't have to encode it:

- **None** (default) — bridge composition isn't checked; the acceptor pattern match is all that's required.
  Use this for direct, bridge-less crosslinks (the DAP-type mechanism above).
- **Glycine** — requires a polyglycine bridge whose length falls within the min/max range you set (e.g. 5
  to 5 for an exact pentaglycine bridge, or a wider range for permissive/novel-discovery searches).
- **D-Asp** — requires a single D-Asp/D-Asn bridge residue.

### Preview Matches

Before running a full analysis, click **Preview Matches** to see which structures in your loaded mass
database actually match your donor and acceptor patterns (and bridge constraint, if set). This is the best
way to catch a pattern mistake before it silently produces a wrong dimer mass — **Run Analysis** stays
disabled until a preview has run successfully with your current patterns.
