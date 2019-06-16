PIPELINE_FOLDER="/data/cbprojectarchive/ADJUSTME/Riboswitches_pipeline"

load PIPELINE_FOLDER + "/config/ribosw.pipeline.groovy"
load PIPELINE_FOLDER + "/config/tool.versions.groovy"

load PIPELINE_FOLDER + "/modules/fastqc.groovy"
load PIPELINE_FOLDER + "/modules/rnafold.groovy"
load PIPELINE_FOLDER + "/modules/fastq_processor_riboswitch.groovy"
load PIPELINE_FOLDER + "/modules/report.groovy"
load PIPELINE_FOLDER + "/modules/collectBpipeLogs.groovy"

run {
  "%.fastq.gz" * [ FastQC , fastq_processor_riboswitch ] + RNAfold + report + collectBpipeLogs
}
