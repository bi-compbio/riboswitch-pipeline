use warnings;
use strict;
use Getopt::Long;

binmode STDOUT, ":utf8";
use utf8;

use JSON;

my $configfile = "";
my $help       = 0;

##
## Get command line arguments
##
sub usage {
    my $message = $_[0];
    my $command = $0;   # get program name from command line
    $command =~ s#^.*/##;
    
    print STDERR (
        $message, 
        "Usage: $command [options]\n" . 
        "Options:\n" . 
        "  --config patterns.json\n" . 
        "  --help this help screen"
    );
    
    die("\n")
}

GetOptions('config=s' => \$configfile,
           'help'     => \$help)
           or usage("Invalid commmand line options.\n");

# validate arguments
usage("") if($help);
usage("Incompatible options.\n") if($help && $configfile);
usage("All arguments are mandatory.\n") if ($configfile eq "");

##
## read JSON config
##
my $json;
{
  local $/; #Enable 'slurp' mode
  open my $fh, "<", $configfile or die "File $configfile not found";
  $json = <$fh>;
  close $fh;
}
my $config = decode_json($json);

##
## Do the magic
##
sub f {
  my $amplicon=shift;
  my $wildcard=shift;
  my $name=shift;
  my $pattern=shift;

  if($amplicon =~ m/$wildcard/) {
    foreach ("a", "c", "t", "g") {
      f($amplicon =~ s/$wildcard/$_/r, $wildcard, $name, $pattern . $_);
    }
  } else {
    print(">", $name, ":", $pattern, "\n");
    print($amplicon, "\n");
  }
}

foreach (@{$config->{'Amplicons'}}) {
  f($_->{"Amplicon"}, $_->{"Wildcard"}, $_->{"Name"}, "");
}

