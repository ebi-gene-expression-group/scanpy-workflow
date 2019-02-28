# Run the steps of the Scanpy workflow

This is a Nextflow Workflow leveraging the [scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts) package to run individual steps of the Scanpy workflow.

## Setup

### Conda/ Bioconda

Workflow dependencies are managed via Conda and Bioconda, so you'll need to set that up, see instructions [here](https://bioconda.github.io/#install-conda). 

### Nextflow

Obviously you'll need Nexflow itself. If you don't have it already you can install via Conda:

```
conda install nextflow
```

You may well want to do this within a Conda environment you create for the purpose.

## Run the workflow

### Inputs

Expected inputs are:

 * A .zip file containing a single directory with the files matrix.mtx, barcodes.tsv and genes.tsv - i.e. the sparse MTX format.
 * A GTF file with gene IDs specifying the features you wish to filter by. 
 
 ### Parameters
 
By default, the workflow will run using parameters from the default configuration file, and with the 'local' executor- i.e. one process at at time. 

You can copy the default configuration, edit the Scanpy and other parameters, and provide it to Nextflow to override any of the settings. See the [Nexflow documentation](https://www.nextflow.io/docs/latest/executor.html) for executor settings.
 
 ### Execution

The workflow can be run directly from the repository:

```
nextflow run -config <your nextflow.config> ebi-gene-expression-group/scanpy-workflow --matrix <mtx zip> --gtf <gtf> --resultsDir <final results dir>
```

This will download the workflow, create any necessary environments, and run the workflow with the specified innputs. Future executions will use a cached copy of the pipeline, should you wish to update the code in future, you can do so like:

```
nextflow pull ebi-gene-expression-group/scanpy-workflow
```

### Outputs

Outputs will be placed in the directory defined as WORKFLOW_RESULTS_DIR under 'env' in nextflow.config ('results' by default). Outputs include:

 * Cluster definitions
 * t-SNE embeddings
 * UMAP coordinates
 * Marker genes

