# Installation

Phyfum consists of two modules: the [java package ](https://github.com/pbousquets/beast-mcmc-flipflop)containing BEAST and our phylogenetic model and the [python package ](https://github.com/pbousquets/phyfum)that we provide as a wrapper to facilitate running the program. Additionally, we have created a [docker image ](https://hub.docker.com/repository/docker/pbousquets/phyfum/general)where everything installed and ready to use. For this reason and the fact that docker images can be ran through singularity on HPC clusters, we encourage users to rely on this containerised version of phyfum.

To install docker, simply run:

```bash
docker pull pbousquets/phyfum:latest
```

To manually install the pipeline, follow the next steps:

{% code lineNumbers="true" %}
```bash
# Download the java package
wget https://github.com/pbousquets/beast-mcmc-flipflop/releases/download/flipflopV1/linux.zip
unzip linux.zip 
echo "export PATH=$PATH:$(pwd)/BEASTv1.0.0/bin" >> ~/.bashrc #To permanently add binaries to PATH

# Install the python package
pip3 install phyfum 
```
{% endcode %}

The first time phyfum is used, the runtime will take longer than usual, since it'll install missing R libraries. Please, make sure [R](https://cran.r-project.org/) is available too prior to running phyfum.
