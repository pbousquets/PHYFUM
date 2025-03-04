# Quick usage

#### Complete mode

If you're working with raw data (IDAT files), you can run phyfum in **complete** mode. It will preprocess the files with [minfi](https://bioconductor.org/packages/release/bioc/html/minfi.html). If needed and if both tumor and normal samples are available, it will also run a copy number analysis with [rascal](https://github.com/crukci-bioinformatics/rascal). This will allow to blacklist fCpGs that won't fluctuate as the model expects.&#x20;

An example run for this workflow would look like this:

{% code overflow="wrap" %}
```bash
phyfum run --input ${input_dir}/epic_array_dir \
    --output experiment1 \
    --workdir experiment1 \
    --patientinfo /path/to/epic_array_dir/sample_sheet.csv \
    --patient-col patient \
    --age-col age \
    --patient-col Patient \ 
    --sample-col Sample \
    --sample-type-col group \
    --stemcells 3-10-3 
```
{% endcode %}

#### Trees mode

In case you already have the beta values, you can run phyfum in **trees** mode. The pipeline will simply deploy the XMLcreator tool to format the input data as expected by [BEAST](https://beast.community/) and run the inference.

{% code overflow="wrap" %}
```bash
phyfum run --input /path/to/exampleBeta.csv \
    --output onlybetas \
    --workdir onlybetas \
    --patientinfo /path/to/meta.csv \
    --patient-col patient \
    --age-col age \
    --sample-type-col group \
    --stemcells 3-10-3
```
{% endcode %}

Phyfum auto-detects what kind of input is provided and selects automatically the optimal workflow.

{% hint style="info" %}
A custom set of fCpG sites can be provided. Otherwise, an additional step to select the fluctuating sites will be ran (at least 2 individuals are required).&#x20;

CpGs fluctuate differently across tissues. Therefore, we recommend not reusing fCpG sets accross different tissues
{% endhint %}
