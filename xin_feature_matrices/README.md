# ML Feature Changes

## v2.5 - 2020_04_18
- Features changes
    - editing_value - CHANGED. Update the value of "Neil1" dataset based on the 6 replicates instead of 2.  
    - Structure-related features - CHANGED. Update the value based on the updated "structure classification methods". 

## v2.4 - 2019_03_14
- Features changes
  - add new feature column ("probability_active_conf" - Ensemble Probability Confirmation) 
  to Neil1, AJUBA and TTYH2_Combined feature sets. 

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
