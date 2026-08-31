# Nocosomu Singularity Containers

## Prerequisites
Building a container from a Singularity definition file natively requires root privileges (`sudo`). If you intend to run the Nocosomu pipeline on an HPC cluster, you must build these images on a local machine where you possess administrative rights. 

## Build Execution
Use the `singularity build` command to compile a definition file (`.def`) into a compressed Singularity Image File (`.sif`). 

```bash
sudo singularity build <image_name>.sif <definition_file>.def
```

## Example command for the Borzoi container:

```vash
sudo singularity build borzoi.sif borzoi.def
```

## Image Inventory
### Base Environments
- R_base.def
- python_base.def
- ubuntuPython_core.def

### Sequence-Based Models
- alphagenome.def
- borzoi.def
- chrombpnet.def
- enformer.def
- flashzoi.def
- sequence_based_models.def
- surecnn.def

### Utilities
- curzua_workRser.def
- helmsman_og.def
- nibbler_torch.def
- ocrmypdf.def
- pangolin_git.def
- pascalx.def
- puffin_resources.def
- spliceai_pip.def
- spring.def
- sra-toolkit.def

## Deployment
Transfer the compiled .sif files to your target cluster environment. Place the .sif files in the container directory specified by your nextflow.config file to ensure the pipeline resolves the required dependencies at runtime.
