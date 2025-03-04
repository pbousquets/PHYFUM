# PHYFUM

[![PyPI - Version](https://img.shields.io/pypi/v/phyfum.svg)](https://pypi.org/project/phyfum)
![PyPI - Downloads](https://img.shields.io/pypi/dm/phyfum)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/phyfum.svg)](https://pypi.org/project/phyfum)
[![GitBook](https://img.shields.io/badge/GitBook-3884FF?logo=gitbook&logoColor=fff)](https://phyfum.gitbook.io/tutorial/)
[![Docker Image Version](https://img.shields.io/docker/v/pbousquets/phyfum?logo=docker&link=)](https://hub.docker.com/r/pbousquets/phyfum)
[![Docker Pulls](https://img.shields.io/docker/pulls/pbousquets/phyfum?logo=docker)](https://hub.docker.com/r/pbousquets/phyfum)


#### Visit our [GitBook](https://phyfum-1.gitbook.io/tutorial/) for a detailed tutorial of Phyfum
--- 

Phyfum is a tool for inferring phylogenetic trees on methylation-based studies. We harness fluctuating CpG (fCpG) sites of methylation arrays to study the clonal evolution of samples. You can read more about fCpGs in the [original paper](https://www.nature.com/articles/s41587-021-01109-w). 

We have implemented a phylogenetic model within [BEAST v.1.8.4] based on the original model described in the paper above. Mind that this repo contains only the BEAST-based program to run the phylogenetic inference. For the full suite of tools and workflow, please refer to the [PHYFUMflow](https://github.com/pbousquets/PHYFUMflow) repo. It covers the IDAT preprocessing, fCpG calling, automatic XML generation, BEAST inference and . Additionally, if both tumor and reference samples are available, CNVs are called to curate non-fluctuating CpGs.


## Build
```
git clone https://github.com/pbousquets/PHYFUM
cd PHYFUM
ant linux
```

## License

`phyfum` is distributed under the terms of the [CC-BY-NC-SA](LICENSE.txt) license.

