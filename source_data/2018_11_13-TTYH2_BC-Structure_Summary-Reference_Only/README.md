# RNA Lib Structure Summary - TTYH2_BC Reference Only

> Date: 2018-11-13

## Dataset

No experimental dataset. Only predicted the "reference" structure of each RNA Library item. 

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
			“reference_structure”: str  // The structure only based on “sequence”
			“reference_bpp”: 2-D array  // The “Base Pair Probability” of each “nt” pair - 2-dimensional array. 
			“A-to-I_editing_site”: number  // The “editing position” based on the “Sequence”.
			“A-to-I_editing_level”: float  // The “editing value”.
		} 
	]
}
```
