// Here you have the bpipe specific configuration of the pipeline. It basically
// defines locations of files and directories relative to the project root.
//
// PIPELINE_FOLDER is the only global variable defined in the pipeline description file (./Riboswitches_pipeline/pipelines/ribosw.groovy)
//
// Don't forget to describe the samples/contrasts/patterns in ./ribosw.json

//
// project directories
//
PROJECT = "/data/cbprojectarchive/ADJUSTME"
LOGS    = PROJECT + "/logs"
QC      = PROJECT + "/qc"
REPORTS = PROJECT + "/reports"
RESULTS = PROJECT + "/results"
TMP     = PROJECT + "/tmp"

//
// tool locations (only tools not loaded via "module load"; look at the tool.versions.groovy for them)
//
TOOL_FASTQ_PROCESSOR_RIBOSWITCH = PIPELINE_FOLDER + "/tools/fastq_processor_riboswitch.pl"
TOOL_ALLVARIANTS      = PIPELINE_FOLDER + "/tools/allvariants.pl"
TOOL_REPORTQC         = PIPELINE_FOLDER + "/tools/reportQC.Rmd"
TOOL_REPORT           = PIPELINE_FOLDER + "/tools/report.Rmd"
TOOL_REPORT_HELPERS   = PIPELINE_FOLDER + "/tools/report.helpers.R"
TOOL_REPORT_CSS       = PIPELINE_FOLDER + "/tools/report.css"
TOOL_VERSIONS         = PIPELINE_FOLDER + "/config/tool.versions.groovy"

//
// module specific configuration
//
RIBOSW_CFG  = PIPELINE_FOLDER + "/config/ribosw.json"
FASTQC_OUT  = QC + "/fastqc"       //where the output files from FastQC will go
RNAFOLD_OUT = RESULTS              //where the output files from RNAfold will go
FASTQ_PROCESSOR_RIBOSW_OUT = RESULTS + "/ribosw" // where the outputs from the fastq processor will go
