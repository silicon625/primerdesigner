# Primer designer version 0.1.2

## Introduction

This python package is designed for searching species-specific primers sets (including internal oligo probe) which have both high amplification efficiency and specificity to target species in genomic scope. Intended application scenario is for designing primers which can specifically detect certain pathogens whereas other genetically closely related pathogens may also appear in the same sample. 

Traditionally, the template gene is specificied before primer design. But in some condition, such a gene which is both suitable for primer design and shows adequate amount of genetical difference may not be so easy to find. Compared with other primer design program, the main advantage of this pipeline is that it automated the progress of finding a suitable gene.

Examples are provided in the "example" folder. Including the directory tree of raw dataset, steps for constructing database, and the usage of pipeline itself.

For detailed information, you can also check the source code directly.

# Installation

Pythonic dependencies and the program itself can be installed by pip command

```
pip install primerdesigner-0.1.2-py3-none-any.whl
```

This program also requires BLAST and mafft to be installed, please check their official release

BLAST: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

mafft: https://mafft.cbrc.jp/alignment/software/


# Dataset Directory Tree Structure and Database Construction
A dataset is a group of files downloaded from ncbi dataset website: https://www.ncbi.nlm.nih.gov/data-hub/genome/. Users can also create their own dataset using genomic sequence fasta and gff annotation, but they must let the dataset have the same structure with the official one.
    A dataset must have the following structure:
    - root_dir
        - assembly_data_report.jsonl
        - dataset_catalog.json
        - \<first assembly folder named with assembly accession\>
            - files
            - ...
        - \<second assembly folder named with assembly accession\>
            - files
            - ...
        - \<other assembly folders named with assembly accession\>
            - files
            - ...
dataset_catlog.json contains the mappings of all files in the folder, this file is also used for detecting a dataset. assembly_data_report.jsonl contains the meta data for all assemblies. They are both not omittable. Be sure that customized datasets contain both of them and their format are correct.

A database can be build from raw dataset by  primerdesigner.database.Database.build(input_dir = "your_input_dir", output_dir = "your_output_dir"). The input and output directory needs to be specified.

Once built, a database can be loaded again quickly by using primerdesigner.database.Database(config_path = "path_to_your_dataset_catlog"). You need to specify path to dataset_catalog.json of your database. 

# Pipeline Parameter Explaination

The main entry point for the pipeline is primerdesigner.primer_designer.find().

it has 6 parameters:

 - db: database object loaded and constructed from raw dataset.

 - include: a list of identifers for genome assemblies which the final primer set output must be able to amplify. Genus name, species name and organism id are allowed. The pipeline will find all assemblies matches these identifiers.

 - exclude: a list of identifers for genome assemblies which the final primer set output must NOT be able to amplify. The looking up stragegy is the same as "include" parameter. Note that if a assembly is both marked as include and exclude, it will be considered as exclude.

 - workers: thread number for BLAST and mafft.

 - pick_probe: whether picking internal oligo probe is necessary. Note that default parameters for primer3 is different when this flag is turned on or off. Processing logic is also a bit different. You should not turn on this option when designing primers without probe.

 - reference_id: In the pipeline, a reference genome will be selected from annotated genomes in "include" group according to the sequencing quality and coverage. The following homologous gene search and primer design is mainly based on the sequence of this assembly. You can also manually specify.


