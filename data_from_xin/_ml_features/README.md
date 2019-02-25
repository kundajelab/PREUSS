# ML Feature Changes

## v2.3 - 2019_02_17

- Features changes
    - sim_nor_score - NEW, a "normalized similarity score" between the RNA isoform and the WT
    - site_1_1 - CHANGED, changed from "site_1_1_internal" which now also includes the "STEM" type

## v2.2 - 2019_02_05

- Add three new features
    - all_stem_length - include "site" if "site" is Stem
- The sequence of AJUBA_BC get "sliced"

## v2 - 2019_01_30

- Add three new features
    - free_energy - calc by "RNA Fold"
    - mut_ref_struct - "Struct type" @ corresponding "mutation position" in "wt"
    - site_1_1_internal - The "nt pair" of the interior loop if "editing site" is a 1:1 interior loop
- Change the wording from "interior" to "internal" for "ALL" feature names. 
