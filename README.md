# DelSIEVE - Deletions-enabled SIngle-cell EVolution Explorer

DelSIEVE is a novel probabilistic model that from raw read counts of scDNA-seq data jointly infers cell phylogeny and somatic variants, including SNVs and deletions. DelSIEVE is an extension to our predecessor model, SIEVE (check out the [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02813-9) and [software](https://github.com/szczurek-lab/SIEVE)). DelSIEVE is implemented as a superset of SIEVE and available as a package of [BEAST 2](https://www.beast2.org/).

For more details, check out our paper in bioRxiv:

[DelSIEVE: cell phylogeny model of single nucleotide variants and deletions from single-cell DNA sequencing data](https://www.biorxiv.org/content/10.1101/2023.09.09.556903v2)

## Installation

There are multiple ways to install DelSIEVE (listed below). We recommend building a Singularity container, but users can also manually install it by following [install by hand steps](https://www.beast2.org/managing-packages/). If you choose to install DelSIEVE by yourself and SIEVE already exists, please remove it before installing DelSIEVE. 

### 1. Use Singularity container

It is highly recommended to build a Singularity container as defined in `exampels/Singularity.def`:

```bash
singularity build --fakeroot /path/to/container exampels/Singularity.def
```

### 2. Use built package

The required and up-to-date compressed file of DelSIEVE for installation can be found under `dist/`.

### 3. Build from source

It is also possible to build DelSIEVE from the source with `ant` according to `build.xml`. Simply run `ant build` under the root directory of the cloned project. Once it is done, you can find the compressed file under `dist/`.

### Compatibility with BEAST 2

Currently, DelSIEVE is compatible with BEAST 2.6. Later DelSIEVE will be adapted for BEAST 2.7.

## Input to DelSIEVE

DelSIEVE requires the following input files:

1. Raw read counts as an output of [DataFilter](https://github.com/szczurek-lab/DataFilter). An example can be found in `examples/data/read_counts.tsv`.
2. Cell names. An example can be found in `examples/data/cell_names`.
3. A template file which has been configured for the analysis.
   BEAST 2 requires a configuration file for the phylogenetic analysis. This configuration file can be edited directly or using Beauti, a GUI tool provided by BEAST 2. Although at this moment DelSIEVE does not support Beauti, we provide some configured examples under `examples/templates`, which were used in the benchmarking and real data analysis in the corresponding paper. Their names follow a pattern of `A_B_C_D.xml`, where:
   - `A`: either `rmc`, for relaxed molecular clock, which requires the [ORC](https://github.com/jordandouglas/ORC) package, or `smc`, for strict molecular clock.
   - `B`: either `mu`, for a subsittution rate matrix only considering mutations, or `mudel`, for a subsitution rate matrix considering both mutations and deletions.
   - `C`: either `ld`, for locus dropout ("LDO" mode), or `sa`, for single alleleic dropout ("ADO" mode).
   - `D`: either `stage_1`, for the first stage (without acquisition bias correction), or `stage_2`, for the second stage (with acquisition bias correction).

## Run DelSIEVE

There are multiple ways to run DelSIEVE (listed below).

We recommend using a two-stage strategy when accurate branch lengths of the cell phylogeny are of the concern. To is end, we have provided a stack of snakemake rules (`examples/Snakefile`) to execute this strategy automatically. A graphical illustration of the procedure can be found in `examples/dag.pdf`. To take advantage of the snakemake pipeline, please create a conda environment using `examples/environment.yml`:
```bash
mamba env create -f examples/environment.yml
``` 

Before running DelSIEVE, please set up the `CONFIGURATION` section in `Snakefile`.

### 1. Use Snakemake with Singularity container

### 2. Use Snakemake with manually installed DelSIEVE 

### 3. Step-by-step guide

The two-stage strategy mentioned above can also be executed step-by-step as follows.

1. DataCollector
   
   DataCollector is used to integrate the selected data from DataFilter with the configured template for the first stage. It can be called using BEAST 2's `appluancher`:
   
   ```bash
   $ /path/to/applauncher DataCollectorLauncher -help
     Usage: datacollector [-help] [-prefix <output_file_prefix>] [-cell <cell_names_file>] [-sciphi] [-data <data_file>] [-datatype <i>] [-ignoreSex] [-template <template_configuration_file>] [-exclude <excluded_cell_names>] [-out <output_file>] [-sample <i1> <i2>] [-bgcs <i1> <i2> <i3>] [-bgfsc <i1> <i2> <i3> <i4> <i5>] [-miss <r>]
      -help option to print this message -> OPTIONAL
      -prefix specifies the prefix of output files (a folder must be ended with '/') -> OPTIONAL
      -cell specifies a blank spaces separated document containing cell names -> MANDATORY
      -sciphi specifies whether the file of cell names is compatible with SCIPhI (default incompatible) -> OPTIONAL
      -data provides a data document ending with .tsv -> MANDATORY
      -datatype specifies input datatype (0: Coverage-Support, 1: Full supports-Coverage; default: 1) -> OPTIONAl
      -ignoreSex ignores candidate mutated sites from sex chromosomes; usually turned on for males -> OPTIONAL
      -template provides a configuration document ending with .xml -> MANDATORY
      -exclude specifies names of the cells to be excluded -> OPTIONAL
      -out specifies the configuration file integrating with the input data -> OPTIONAL
      -sample samples a part of the input data for test purposes; the first number defines the number of sampled cells, and the second number defines the number of sampled loci; by default, all the data will be loaded -> OPTIONAL
      -bgcs specifies the order of background information for "Coverage-Support" datatype w.r.t. (coverage 0, variant 1, normal 2); default: 0 1 2; working with "-datatype 0" if specified -> OPTIONAL
      -bgfsc specifies the order of background information for "Full support-Coverage" datatype w.r.t. (0 - variant1, 1 - variant2, 2 - variant3, 3 - normal, 4 - coverage); default: 0 1 2 3 4; working with "-datatype 1" if specified -> OPTIONAL
      -miss specifies the threshold of missing data percentage; set to 1 to turn off the selection; the programme will choose as many sites as possible without exceeding the specified threshold -> OPTIONAL
   ```
   
   Note that options marked by "MANDATORY" are required for running the command, while those by "OPTIONAL" are not.
   
2. BEAST 2
   
   With the output of step 1, now it is the time to run the phylogenetic analysis with the following command:
   
   ```bash
   $ /path/to/beast -overwrite -threads {NUM} -prefix {PREFIX} {CONFIGURED FILE}
   ```
   
3. TreeAnnotator
   
   TreeAnnotator is used to summarise a maximum clade credibility tree (MCC) from the sampled trees in step 2.
   
   ```bash
   $ /path/to/applauncher SCSTreeAnnotatorLauncher -help
     Usage: scstreeannotator [-prefix <output_file_prefix>] [-heights <keep|median|mean|ca>] [-burnin <i>] [-b <i>] [-limit <r>] [-target <target_file_name>] [-help] [-forceDiscrete] [-lowMem] [-hpd2D <r>] [-nohpd2D] [-noSA] [-simpleTree] [-noTrunk] <input-file-name> [<output-file-name>]
      -prefix specifies the prefix of output files (a folder must be ended with '/') -> OPTIONAL
      -heights an option of 'keep' (default), 'median', 'mean' or 'ca' -> OPTIONAL
      -burnin the percentage of states to be considered as 'burn-in' -> MANDATORY, the same as -b"
      -b the percentage of states to be considered as 'burn-in' -> MANDATORY, the same as -burnin
      -limit the minimum posterior probability for a node to be annotated -> OPTIONAL
      -target specifies a user target tree to be annotated -> OPTIONAL
      -help option to print this message -> OPTIONAL
      -forceDiscrete forces integer traits to be treated as discrete traits -> OPTIONAL
      -lowMem use less memory, which is a bit slower -> OPTIONAL
      -hpd2D the HPD interval to be used for the bivariate traits -> OPTIONAL
      -nohpd2D suppress calculation of HPD intervals for the bivariate traits -> OPTIONAL
      -noSA interpret the tree set as begin from a not being from a sampled ancestor analysis, even if there are zero branch lengths in the tree set -> OPTIONAL
      -simpleTree simple output tree only containing tree heights (the output file will be labeled with 'simple'), default disabled -> OPTIONAL
      -noTrunk trees processed do not have a trunk connecting the root of samples and the root of entire tree, default disabled -> OPTIONAL
   ```
   
4. VariantCaller
   
   VariantCaller is used to call variants given the outputs of step 1-3.
   
   ```bash
   $ /path/to/applauncher VariantCallerLauncher -help
     Usage: variantcaller [-help] [-threads <i>] [-prefix <output_file_prefix>] [-runconfig <run_config_file>] [-burnin <i>] [-b <i>] [-cached <cached_estimates_file>] [-estimates <mean|median|mode|all>] [-kdedist <gaussian|epanechnikov|rectangular|triangular|biweight|cosine|optcosine>] [-mcmclog <mcmc_log_file>] [-config <config_file>] [-tree <tree_file>] [-details] [-cells <i>] [-medianrate] [-userate]
      -help option to print this message -> OPTIONAL
      -threads  specifies the number of threads (default 1); recommending to use more -> OPTIONAL
      -prefix specifies the prefix of output files (a folder must be ended with '/') -> OPTIONAL
      -runconfig specifies a modified configuration file and calls variants -> OPTIONAL, having priority over other flags except for -threads and -prefix
      -burnin specifies the percentage of samples to be considered as 'burn-in' -> MANDATORY (or -b)
      -b the same as 'burn-in' -> MANDATORY (or -burnin)
      -cached specifies the cached estimates file -> OPTIONAL, conflicting with -estimates, -kdedist, -mcmclog, -allelic
      -estimates specifies which kind of estimates of samples will be used to perform variant calling (one of 'mean', 'median' (default), 'mode', and 'all', where the last option compares likelihoods between all situations and chooses the largest one) -> OPTIONAL, conflicting with -cached
      -kdedist specifies which KDE distribution will be used if mode estimates is selected (one of 'gaussian' (default), 'epanechnikov', 'rectangular', 'triangular', 'biweight', 'cosine', 'optcosine') -> OPTIONAL, conflicting with -cached
      -mcmclog specifies the log file containing MCMC samples -> OPTIONAL, conflicting with -cached
      -config specifies the configuration file used to performing phylogenetic analisys (either xml or json) -> MANDATORY
      -tree specifies the tree summarized from ScsTreeAnnotator -> MANDATORY
      -details marks whether to save details, including inferred genotypes and ternary matrix, if specified -> OPTIONAL
      -cells specifies the number of mutated cells used to filter invariant sites, default 1 -> OPTIONAL
      -medianrate use median rate rather than mean rate from the input tree -> OPTIONAL
      -userate use both branch length and rate from the input tree in variant calling -> OPTIONAL
   ```
   
5. Update the configuration file for the second stage. 
   
   Step 1-4 are for the first stage. Given the outputs of the above steps and a template for the second stage, we can update such a template with a python script provided by us (`examples/scripts/set_up_stage2_from_stage1.py`):
   
   ```bash
   $ python3 set_up_stage2_from_stage1.py -h
   usage: set_up_stage2_from_stage1.py [-h] [--tree TREE] --estimates ESTIMATES --results RESULTS --template1 TEMPLATE1 --template2 TEMPLATE2 --out OUT
   
   Setup configuration files for stage 2 based on results of stage 1.
   
   optional arguments:
     -h, --help            show this help message and exit
     --tree TREE           starting tree for stage 2 with the MCC tree from stage 1
     --estimates ESTIMATES
                           estimated parameter values from stage 1
     --results RESULTS     a file containing the a list of files of variants called from stage 1
     --template1 TEMPLATE1
                           template configuration file for stage 1
     --template2 TEMPLATE2
                           template configuration file for stage 2
     --out OUT             modified template configuration file for stage 2
   ```
   
6. With the updated template for the second stage, it is ready to run by repeating step 1-4. 

   Note that in order to map mutations (see [the following section](#mutation-mapping)) later, option `-details` in VariantCaller (step 4) must be explicitly specified.

## Mutation mapping

After getting variant calling results, SIEVE is able to map the mutations to branches. But first, we need to build mutation maps using tools like [Annovar](https://annovar.openbioinformatics.org/en/latest/), with the output of which SIEVE is compatible. For example, in the SIEVE paper we called Annovar with the following [command](https://davetang.org/wiki2/index.php?title=ANNOVAR):

```bash
$ perl path/to/table_annovar.pl path/to/vcf/file path/to/database/directory -buildver hg19 -out {PREFIX} -remove -protocol refgene,ensGene,dbnsfp30a,ljb26_all -operation gx,gx,f,f -nastring . -vcfinput -polish -xref path/to/gene_fullxref.txt
```

and we used four databases *refGene*, *ensGene*, *dbnsfp30a*, and *ljb26_all* of the *hg19* reference genome. The list of available databases of Annovar can be found [here](http://annovar.openbioinformatics.org/en/latest/user-guide/download/).

Besides, SIEVE is able to process another format of mutation maps where each line contains a chromosome label, a starting position number, an ending position number and a gene name, separated by any white spaces except new line.

Moreover, a list of disease-related genes from sources like [COSMIC](https://cancer.sanger.ac.uk/cosmic) can also be used to filter out unrelated/uninterested mutations. 

GeneAnnotator in SIEVE is developed for mutation mapping:

```bash
$ /path/to/applauncher GeneAnnotatorLauncher -help
  Usage: geneannotator [-help] [-details <variant_caller_details>] [-subst <i>] [-map <mutation_map_file>] [-anv <i>] [-anvfiltercol <i>] [-anvfilterkws <annovar_filtering_keywords>] [-filter <filtering_genes_file>] [-sep <i>] [-col <i>] [-out <out_tree_file>]
    -help option to print this message -> OPTIONAL
    -details specifies a directory containing the output of VariantCaller called with option '-details'. -> MANDATORY
    -subst specifies the substitution model used to infer phylogeny and call variants (0: ScsFiniteMuExtendedModel, 1: ScsFiniteMuDelModel). -> MANDATORY
    -map specifies the mutation map. -> MANDATORY
    -anv specifies when the value of --map option is an output of Annovar (0: white-space-separated, 1: comma-separated). -> OPTIONAL
    -anvfiltercol specifies which column of the output of Annovar should be used to filter genes; only working when --anv is specified; setting to -1 for non-filtering (0-based; default: -1). -> OPTIONAL
    -anvfilterkws specifies the comma-separated keywords to keep genes, with the first keyword being either 'f' for 'fixed' or 'r' for 'regular expression'. -> OPTIONAL
    -filter specifies filtering genes file (headers required). -> OPTIONAL
    -sep specifies the separator of filtering genes file (0: tab-separated; 1: comma-separated; default: 0). -> OPTIONAL
    -col specifies the column index of gene names in the filtering genes file (0-based; default: 0). -> OPTIONAL
    -out specifies the output tree file with genes annotated. -> OPTIONAL
```

## Troubleshooting

1. `java.lang.OutOfMemoryError: Java heap space`
    This could happen when the input data contains many candidate SNV sites (e.g., > 5k). 
    On Windows, macOS, and Linux with GUI, check the solutions [here](https://www.beast2.org/increasing-memory-usage/). 
    If this error appears on Linux server without GUI, the easiest solution without root privilege is to simply run `export _JAVA_OPTIONS="-Xmx10g"`, `10g` can be increased as you need as long as it is within the server's capacity.

