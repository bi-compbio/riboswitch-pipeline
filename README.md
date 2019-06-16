# Riboswitches pipeline

Sergi Sayols Puig, Holger Klein

Computational Biology Group

(c) 2016-2019, Boehringer Ingelheim Pharma GmbH & Co KG

License: GPL v3 (see LICENSE)





Process fastq files and report differences in abundance.

## TOC

* [Description](#description)
* [Usage instructions](#usage-instructions)
  * [For the impatients](#for-the-impatients)
  * [For the rest](#for-the-rest)
  * [If something went wrong](#if-something-went-wrong)
  * [Running the pipeline in the cluster](#running-the-pipeline-in-the-cluster)
* [Configuring ribosw.groovy](#configuring-ribosw.groovy)
  * [Targets](#targets)
  * [Contrasts](#contrasts)
  * [Model matrix](#model-matrix)
  * [Spike in](#spike-in)
  * [Patterns](#patterns)
* [Outputs](#outputs)
* [Current limitations](#current-limitations)

## Description

The pipeline has the following structure:

```
.
├── config
├── deps
│   └── JSON-2.94
├── modules
├── pipelines
└── tools
```

- the **`config`** folder contains configuration files to customize the pipeline for the current project/run. Here one will setup the patterns to look for, the targets and contrasts for differential expression, the parameters for the cluster resource manager, etc.
- the **`deps`** folder contains external dependencies, like the Perl library to read JSON.
- the **`modules`** folder contains the bpipe wrappers to call the tools. The wrappers are basically rules that specify how to get the output files from the input files, usually just calling a bash script.
- the **`pipelines`** folder contains the definition of the pipeline (what steps go after what). It has it's own folder as one could create different pipelines.
- the **`tools`** folder contains the actual tools. Tools already available in the cluster, like Perl or R, are not here. This folder is meant to store tools specific for the pipeline, like the fastq_processor or the reports.

More in detail, currently the pipeline contains:

```
.
├── config                      --> Pipeline configuration files
│   ├── bpipe.config            --> SGE resources needed to run the tasks
│   ├── ribosw.json             --> definition of RE (patterns) and edgeR parameters (targets, contrasts, model matrix...)
│   ├── ribosw.pipeline.groovy  --> customization parameters for the pipeline (location of directories and files)
│   └── tool.versions.groovy    --> tool's versions to be loaded via module load
├── deps                        --> external dependencies
│   └── JSON-2.94               --> the fastq_processor_riboswitch needs json to read the RE (patterns)
│       ├── <...> non relevant files here
│     <...>
├── modules                     --> the pipeline is modular. Here are the modules (programs) the pipeline calls
│   ├── collectBpipeLogs.groovy --> get the logs (stderr) from ./.bpipe and move them to the ./logs directory
│   ├── fastqc.groovy           --> call fastq
│   ├── fastq_processor_riboswitch.groovy --> count RE (patterns) in a fastq file
│   └── report.groovy           --> do the QC and differential expression analysis, and present the results as HTML
├── pipelines                   --> different Riboswitches-related pipelines
│   └── ribosw.groovy           --> currently only 1 pipeline, which just counts patterns and does the DE analysis
└── tools                       --> here lie the real tools the pipeline calls (apart from those loaded via module load)
    ├── fastq_processor_riboswitch.pl --> count RE (patterns) in a fastq file
    ├── report.css              --> CSS stylesheet for the reports (both use the same)
    ├── report.helpers.R        --> functions library for the reports
    ├── reportQC.Rmd            --> do the QC analysis, and present the results as HTML
    └── report.Rmd              --> do the differential expression analysis, and present the results as HTML
```


## Usage instructions

### For the impatients

- Create an empty project folder, and get in.
- Create a `rawdata` folder and copy there all your fastq(.gz) files.
- Clone this repository, with `git clone git@git.eu.boehringer.com:bibc_compbio/Riboswitches_pipeline.git`.
- Configure the pipeline. The file `./Riboswitches_pipeline/config/ribosw.json` contains most of the relevant configuration:
  - Define the patterns to look for in the `Patterns` section. Spikein names _must_ start with "Spikein". Side note HK: the patterns can be derived from the amplicons defined by Benjamin Strobel in the project sheet. Usually I take the 5bp upstream and downstream of the two variable regions defined by consecutive N, with a spacer adjusted to their distance. For the spikein I use 20 consecutive bp specific for the spikein sequences (to have the same error probabilities for spikeins as for the 4x5 bp in the riboswitches).
  - Define your samples in the `Targets` section.
  - Define your contrasts of interest in the `Contrasts` section.
  - Define the model matrix to use in edgeR in the `ModelMatrix` section.
  - If you'd like to normalize using a Spike-in, define it in the `UseSpikein` section.
- Configure the project folder, and also customize the project subtleties in `ribosw.pipeline.groovy`. Usually the only relevant parameter is the project folder. Look for the string ADJUSTME
- Tell the pipeline where the modules it has to call are located: `./Riboswitches_pipeline/pipelines/ribosw.groovy`. Look for the string ADJUSTME.
- Edit the static text of the reports in `./Riboswitches_pipeline/tools/report.Rmd` and `./Riboswitches_pipeline/tools/reportQC.Rmd`.
- Run the pipeline from the root of the project folder with:
```
$ module load Bpipe
$ bpipe run -n 8 ./Riboswitches_pipeline/pipelines/ribosw.groovy ./rawdata/*.fastq.gz
```

The `-n` parameter defines how many files in parallel will be processed. If everything run smooth, you should have a `results` folder with the raw and normalized counts and the results of the differential expression analysis, and a `reports` folder with a couple of html reports from the QC and differential expression analysis.

### For the rest

Read the instructions again. And have a look at the [Configuring ribosw.groovy](#configuring-ribosw.groovy) section below.

### If something went wrong

Don't panic.

Usually it's just a configuration issue.

- Check the targets are all set. Missing targets or expected counts files will break the pipeline.
- Check the contrasts are properly named. `groupA` and `groupB` must follow the edgeR expected pattern: \<targets_tuple_field_name\>\<factor\> pasted together.
- Check the output counts in the `./results/ribosw` folder. If the RE pattern was wrongly written could be you didn't match anything, and drive edgeR crazy.
- Check your `ribosw.json` file is still a valid json file.

Other sources of help, in order to trace back the error:

- Bpipe creates a file called `commandlog.txt`, which contains all the programs command lines called for every step and sample. It's a mess, since Bpipe runs several files in parallel and the contents of this file are sequential, but it may help to trace thing going wrong.
- Bpipe creates a `./.bpipe/commandtmp` folder with several subfolders, which contain the shell command, stdout and stderr called for each step and sample. The number of the task (folder) can be retrieved from the actual step+sample.properties files located in the `./.bpipe/outputs` folder (just "cat" the file).

### Running the pipeline in the cluster

Two options:

  - The easy way: wrap the `bpipe run` call within a script, to be submitted to SGE. The `bpipe run -n 32` parameter controls how many tasks in parallel bpipe will run: use it to set the required cpu resources in SGE.
  - The right way: let bpipe talk to SGE. Copy the file `./Riboswitches_pipeline/config/bpipe.config` to the root folder of the project and run bpipe normally. Bpipe will prepare one shell script per file+step to run in SGE. Adjust the task cpu/mem requirements in the `bpipe.config` file. This is the desidered way of running the pipeline, but didn't have the time to test it in the BI cluster.

## Configuring ribosw.groovy

This is a valid json file containing sections that describe the input files, the patterns to look for, and the behaviour of edgeR to do the differential abundance analysis.

**Remember:** JSON doesn't allow comments.

The sections need to be carefully filled out prior to running the pipeline. Follows a *brief* description of what you'll find there.

### Targets

This is a vector of tuples that describe each sample that will be used in the edgeR step. Not all fields need to be filled out, but **they're expected to be there**. Super important properties that need to be filled out are sampleName, treatment and filePrefix.

```
  "Targets" :
  [
    {
      "projectID" : "531",
      "sampleID"  : "1",
      "sampleName": "S1",
      "group"     : "HEK293",
      "treatment" : "untreated",
      "subject"   : "1",
      "filePrefix": "531_1_S1_R1_001"
    },
    {
      "projectID" : "531",
      "sampleID"  : "2",
      "sampleName": "S2",
      "group"     : "HEK293",
      "treatment" : "tetracycline_30uM",
      "subject"   : "1",
      "filePrefix": "531_2_S2_R1_001"
    }
  ]
```

### Contrasts

This is a vector of tuples that describe each comparison that edgeR has to carry out.

The `Name` field is used in several places of the differential expression analysis, basically to name the results. It's free text.

`groupA` and `groupB` are the 2 factors that edgeR will test pairwise (A-B). It has to **strictly** follow the naming convention that edgeR understands: \<targets_tuple_field_name\>\<factor\>.

```
  "Contrasts":
  [
    {
      "Name"  : "tetracyclyne_30uM_vs_untreated",
      "groupA": "treatmenttetracycline_30uM",
      "groupB": "treatmentuntreated"
    },
    {
      "Name"  : "tetracyclyne_100uM_vs_untreated",
      "groupA": "treatmenttetracycline_100uM",
      "groupB": "treatmentuntreated"
    },
    {
      "Name"  : "tetracyclyne_100uM_vs_tetracycline_30uM",
      "groupA": "treatmenttetracycline_100uM",
      "groupB": "treatmenttetracycline_30uM"
    }
  ]
```

### Model matrix

This is a single character variable that contains a valid formula that edgeR understands. It's used to describe the linear model edgeR uses to test the contrasts.

```
  "ModelMatrix": "~ treatment + 0"
```

### Spike in

Defines the spike-in pattern to use to normalize the counts before edgeR. Currently, it supports **only** 1 spike-in. Has to be one of the pattern names defined in the `Patterns` section.

```
  "UseSpikein": "Spikein1"
```

### Patterns

This is a vector of tuples that describe each pattern to look for. Contains valid Perl RE.

- They're searched in parallel.
- The names must be unique, otherwise you'll increment the same counter for two different RE.
- The fastq_processor stops looking for more RE once one matches. Thus, the order in this file matters.
- You can define a Spikein sequence, which is output in the global \*\_stats.txt counts.
- Spikein names _must_ start with "Spikein".

```
  "Patterns":
  [
    {
      "Name": "Orig",
      "RE"  : "GATTC([ACGT]{5})AAAAC.{37}CACCT([ACGT]{2})TACAT"
    },
    {
      "Name": "TypI",
      "RE"  : "CGAAA([ACGT]{3})TCGGG.{62}GGGAT([ACGT]{3})TCCTG"
    },
    {
      "Name": "TypII",
      "RE"  : "GATTC([ACGT]{5})TCGGG.{62}GGGAT([ACGT]{2})TACAT"
    },
    {
      "Name": "Spikein",
      "RE"  : "AGTTGTTTTTGTTTTTAATT"
    }
  ]
```

A good advice is to try them with a toy sequence before running the pipeline.

## Outputs

The reports folder will contain a report with the QC, and another with the differential abundance analysis (and some other QC metrics that need from this analysis).

The results has the tables and R objects for the downstream analysis:

```
results
├── counts_raw.xls                         --> matrix variants x samples with the raw counts
├── counts_counts.xls                      --> matrix variants x samples with the edgeR normalized counts
├── counts+results.xls                     --> megamatrix with all counts for all motifs and results of the DE analysis
├── edger_de_results.RData                 --> list (of motifs) with the DGEList=edgeR object, LRT=test results for each contrast
├── edger_de_results.xls                   --> only the LRT results for each contrast/motif, for human beings.
└── ribosw                                 --> results saved by the fastq_processor:
    ├── 531_10_S10_R1_001.txt              --> Pattern + Var1 + Var2 match for each sequence of the fasta
    ├── 531_10_S10_R1_001_full.txt         --> Pattern + Var1 + Var2 + Fasta for each sequence of the fasta
    ├── 531_10_S10_R1_001_counts.txt       --> counts per variant (Pattern + Var1 + Var2)
    ├── 531_10_S10_R1_001_stats.txt        --> summary of matches per pattern
    └── 531_10_S10_R1_001_unmatched.txt    --> sequences from the fasta with no match
```

## Current limitations

- To keep the differential abundance analysis easy, currently it accepts only pairwise comparisons of the style A-B. More complex formulas for the linear model can be provided via the `ModelMatrix` pipeline parameter, although the default `ModelMatrix = "~ treatment + 0"` should fine for the standard case.
- The interactive plots can be damned slow if there are many contrasts or combinations of switches. If this is the case, edit the report and substitute the calls to `ggplotly()` to `print()`, and the call to `pairsD3()` to `pairs()`.
- The Project Info section of the report needs to be edited by hand.
