# ML Feature Changes

## v2.2 - 2019_02_05_

- Add three new features
    - all_stem_length - include "site" if "site" is Stem
- The sequence of AJUBA_BC get "sliced"

## v2 - 2019_01_30

- Add three new features
    - free_energy - calc by "RNA Fold"
    - mut_ref_struct - "Struct type" @ corresponding "mutation position" in "wt"
    - site_1_1_internal - The "nt pair" of the interior loop if "editing site" is a 1:1 interior loop
- Change the wording from "interior" to "internal" for "ALL" feature names. 
