---
layout:
  title:
    visible: true
  description:
    visible: true
  tableOfContents:
    visible: true
  outline:
    visible: true
  pagination:
    visible: false
---

# Phyfum

Phyfum is a tool for inferring phylogenetic trees on methylation-based studies. We harness fluctuating CpG (fCpG) sites of methylation arrays to study the clonal evolution of samples. You can read more about fCpGs in the [original paper](https://www.nature.com/articles/s41587-021-01109-w).&#x20;

We have implemented a phylogenetic model within BEAST v.1.8.4 based on the original described in the above paper. We have also designed a snakemake-based pipeline, covering the IDAT preprocessing, fCpG calling, automatic XML generation and BEAST inference. Additionally, if both tumor and reference samples are available, CNVs are called to curate non-fluctuating CpGs.

Before you start using the tool, we recommend you reading [the-workflow.md](the-workflow.md "mention") section to fully understand how Phyfum works.
