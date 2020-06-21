# RNA Library Analysis

This repo documents all the "scripts" that are used to help run the analysis of "RNA Library" projects. 


## Repo Structure

The repo includes two major parts: 

- Scripts - All the scripts to help run the analysis. 
- Notebook - Organize the "analysis code" in a favor of "notebook". 

## Dependencies

Besides the dependencies defined in `requirements.txt`, it also depends a "local" library - neoRNA helps interact 
with DNA(RNA mainly) elements - sequence, structure, etc as well as some basic operations. 

Here is the way to import the project `neoRNA` 

```bash
# Include the "path" of "neoRNA" into the `PYTHONPATH`. Ex: 
export PYTHONPATH=$PYTHONPATH:$HOME/Documents/_code_repo/geno/neoRNA
```

## Usage

### Machine Learning Related

#### Generate "Feature" File

```bash

# `cd` into "repo root"
./scripts/bash_scripts/gen_features.sh 
```

#### Generate `rmarkdown` Notebook for Random Forest Method

```bash

# `cd` into "repo root"
Rscript ./scripts/r_scripts/gen_notebook.R 

```


## References

### Jupyter (a.k.a IPython Notebook)

Jupyter provides a "interactive" env. for writing Pyhton code and viewing data analysis result in one place. 

 