# Description

This is a nanopore analysis pipeline for variant calling.

It is provided as the source code for the methods for this article: [TBA].

For [TBA] usage and coverage analysis, please see: [TBA].

# Usage

The pipeline was made to run on SLURM, and server specific settings have been removed.
`nanopore_slurm_pipeline.py` handles alignment and variant calling. The .sh files handle the rest of the .slurm files in `slurmScripts/`. Bash scripts have a header describing command-line arguments, and python scripts use `argparse`, `nanopore_slurm_pipeline.py` use the config files.

The config files should be edited to fit the input data. See config/config.sh for an example.

