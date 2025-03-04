# The workflow

We have written a workflow based on snakemake to cover the end-to-end analysis. Users can provide the raw methylation files and the sample sheet as a starting point, and the pipeline will automatically:

1. Preprocess the raw files to extract the beta values with [`minfi`](https://doi.org/doi:10.18129/B9.bioc.minfi).
2. If reference samples available, call CNVs with [conumee](https://doi.org/doi:10.18129/B9.bioc.conumee) to blacklist potentially disruptive CpGs.
3. Call fluctuating CpG sites (If no fCpG site list is provided).
4. Create the XML file required by BEAST (alternative starting point if beta-values are already given - _trees mode_).
5. Run [`BEAST`](https://beast.community/).
6. Automatic model selection
7. Summarise results and draw the phylogenetic trees.

## End-to-End Methylation Analysis Workflow

### Overview

Welcome to our end-to-end methylation analysis workflow powered by Snakemake! This guide will walk you through processing raw methylation files to generate phylogenetic trees.

### Prerequisites

Ensure you have the following installed:

* Python
* Snakemake
* R with `minfi` and `conumee` packages
* BEAST

### Getting Started

#### Input Files

You will need:

* Raw methylation files
* Sample sheet

#### Steps

1. **Preprocessing**: Extract beta values from raw files using [`minfi`](https://doi.org/doi:10.18129/B9.bioc.minfi).
2. **CNV Calling**: Call CNVs with [`conumee`](https://doi.org/doi:10.18129/B9.bioc.conumee) if reference samples are available to blacklist potentially disruptive CpGs.
3. **CpG Site Calling**: Identify fluctuating CpG sites. If no fCpG site list is provided, the pipeline will generate one.
4. **Generate XML for BEAST**: Create the XML file required by BEAST for beta-value inputs (_trees mode_).
5. **Run BEAST**: Execute [`BEAST`](https://beast.community/) for phylogenetic analysis.
6. QC: mixing and convergence checks.
7. **Model Selection**: The pipeline will automatically select the appropriate model.
8. **Computation of integral solution**. Integrate over each number of stem cell to compute the integrate.&#x20;
9. **Result Summarization**: Summarize results and draw phylogenetic trees.

### Running the Pipeline

To run the pipeline, execute:

```bash
phyfum run --input exampleBeta.csv --patientinfo metadata.csv --output test  --nchains 2 --stemcells 3-6-1 --workdir test --mle-ps 
```



### Support

For any issues or questions, please reach out to our support team or check the documentation.

Enjoy your analysis!

***

This tutorial should help you navigate our methylation analysis workflow seamlessly. Happy analyzing!

<figure><img src=".gitbook/assets/f1.png" alt=""><figcaption><p>PHYFUM's workflow</p></figcaption></figure>



