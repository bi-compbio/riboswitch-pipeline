report = {
    doc title: "report",
    desc:  "makes the required pipeline environment variables available to the  markdown report",
    constraints: "reports must be called report.Rmd and reportQC.Rmd; needs pandoc > 1.12",
    bpipe_version: "tested with bpipe 0.9.8.7",
    author: "Sergi Sayols"
    
    output.dir = REPORTS
    
    produce("report.html", "reportQC.html") {
        exec """
            cp $TOOL_REPORT $REPORTS         &&
            cp $TOOL_REPORTQC $REPORTS       &&
            cp $TOOL_REPORT_HELPERS $REPORTS &&
            cp $TOOL_REPORT_CSS $REPORTS     &&

            PROJECT_NAME=\$(basename ${PROJECT})                                 &&
            sed -i "2,2s/REPORT_PROJECT/\${PROJECT_NAME}/" $REPORTS/report.Rmd   &&
            sed -i "2,2s/REPORT_PROJECT/\${PROJECT_NAME}/" $REPORTS/reportQC.Rmd &&

            echo "PROJECT=${PROJECT}" >  $REPORTS/report.conf &&
            echo "QC=${QC}"           >> $REPORTS/report.conf &&
            echo "LOGS=${LOGS}"       >> $REPORTS/report.conf &&
            echo "REPORTS=${REPORTS}" >> $REPORTS/report.conf &&
            echo "RESULTS=${RESULTS}" >> $REPORTS/report.conf &&
            echo "RIBOSW_CFG=${RIBOSW_CFG}"   >> $REPORTS/report.conf &&
            echo "FASTQC_OUT=${FASTQC_OUT}"   >> $REPORTS/report.conf &&
            echo "RNAFOLD_OUT=${RNAFOLD_OUT}" >> $REPORTS/report.conf &&
            echo "FASTQ_PROCESSOR_RIBOSW_OUT=${FASTQ_PROCESSOR_RIBOSW_OUT}" >> $REPORTS/report.conf &&

            module load Bioconductor/${Bioconductor_VERSION} &&
            module load pandoc/${pandoc_VERSION} &&
            Rscript -e 'setwd("$REPORTS"); rmarkdown::render("report.Rmd")'   &&
            Rscript -e 'setwd("$REPORTS"); rmarkdown::render("reportQC.Rmd")'
        ""","report"
    }
}

