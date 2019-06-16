RNAfold  = {
    doc title: "RNAfold",
    desc:  "generate an exhaustive list of all possible variants, and predict the secondary structure",
    constraints: "ViennaRNA is installed as a part of LocARNA",
    bpipe_version: "tested with bpipe 0.9.8.7",
    author: "Sergi Sayols"

    output.dir = RNAFOLD_OUT

    produce("rnafold.txt") {
      exec """
        module load LocARNA/${LOCARNA_VERSION} &&
        export PERL5LIB="${PIPELINE_FOLDER}/deps/JSON-2.94/share/perl5:\$PERL5LIB" &&
        perl $TOOL_ALLVARIANTS --config $RIBOSW_CFG | RNAfold --noPS | paste -d"\\t" - - - | sed -E 's/^(.+)\\s/\\1\\t/' > $output
      ""","RNAfold"
    }
}
