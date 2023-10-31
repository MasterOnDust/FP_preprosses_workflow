# Snakemake Workflow for preprossing of FLEXPART and FLEXDUST output

The workflow expect the output format to be accoring to model setup using the [FLEXPART-script](https://github.com/MasterOnDust/FLEXPART-script) python command line tool.  

### Installation

1. Setup conda enviroment
```shell
conda env create -p ./dust -f environment.yml 
```
Activate conda enviroment

```shell
conda activate ./dust
```

### Running the workflow

The workflow can be configured in the `config/config.yaml` file. Most important is to set the paths correctly. 

1. Test that all the files are there:

```
snakemake -n make_source_contrib_march
```

2. If snakemake is able to resolve all the paths then run the workflow by:
```
snakemake -j2 make_source_contrib_march
```

Please have look at what the rules does in the Snakefile. 

If you have questions about how about this setup and want to setup a similar workflow the please get in touch at ovehaugv@outlook.com 