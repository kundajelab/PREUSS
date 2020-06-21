# RNA Lib Structure Summary - Neil1, Combined Dataset

> Date: 2018-10-15

## Dataset

The “combined” dataset is based on Neil1’s dataset - `2017-08-07`, which combines the “two” replicates. 

## Summary Structure

The summary result is in the `JSON` format, which has the following structure: 

```json
{
	// Basic info of RNA Library, like “title”, “code”, etc
	“indel_rna_id”: str // The info about “intel” RNA items.
	“wt_id”: number // The “RNA ID” of Wild-Type item
	“items”: [  // List of RNA items
		{
			// Basic Info of RNA Item, like “RNA ID”, “Sequence”.
			“reactivity”: list  // The list of reactivity of each “nt”.
			“reference_structure”: str  // The structure only based on “sequence”
			“reference_bpp”: 2-D array  // The “Base Pair Probability” of each “nt” pair - 2-dimensional array. 
			“inferred_structure”: str  // The structure based on “sequence” + “reactivity”.
			“inferred_bpp”: 2-D array  // The “Base Pair Probability” of each “nt” pair - 2-dimensional array. 
			“bootstrap_structures”: The inferred structure for each “bootstrap run” (It has total 1000 runs).  
			“A-to-I_editing_site”: number  // The “editing position” based on the “Sequence”.
			“A-to-I_editing_level”: float  // The “editing value”.
		} 
	]
}
```
