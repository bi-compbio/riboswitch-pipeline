fastq_processor_riboswitch  = {
    doc title: "fastq_processor_riboswitch",
    desc:  "Find patterns in a fastq file and report stats",
    constraints: "",
    bpipe_version: "tested with bpipe 0.9.8.7",
    author: "Sergi Sayols"

    output.dir = FASTQ_PROCESSOR_RIBOSW_OUT

    transform(".fastq.gz") to(".txt") {
        exec """
            export PERL5LIB="${PIPELINE_FOLDER}/deps/JSON-2.94/share/perl5:\$PERL5LIB" &&
            perl $TOOL_FASTQ_PROCESSOR_RIBOSWITCH --infile $input --outdir $FASTQ_PROCESSOR_RIBOSW_OUT  --outprefix \$(basename $input.prefix.prefix) --config $RIBOSW_CFG &&
            OUTFILE=$output &&
            sort $output | uniq -c | sed 's/^\\s\\+//' | tr -s ' ' '\\t' > \${OUTFILE%.txt}_counts.txt &&
            gzip \${OUTFILE%.txt}_full.txt \${OUTFILE%.txt}_unmatched.txt   
        ""","fastq_processor_riboswitch"
    }

    forward input
}
