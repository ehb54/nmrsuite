#!/usr/bin/perl

$notes = "usage: $0 dir

finds run jobs under specified dir from _log,
organizes by executable
and pretty prints json inputs

need jq

";


$d = shift || die $notes;

$cmd = qq[find $d -name "_input_*"];
print "$cmd\n";

@f = `$cmd`;

grep chomp, @f;

# print join ( "\n", @f ) . "\n";
$dashes = '-'x80;
$eqs    = '='x80;

while ( $f = shift @f ) {
    ( my $id ) = $f =~ /_input_(.*)$/;
    my $cf = $f;
    $cf =~ s/_input_/_cmds_/;
    next if !-e $cf;
    my $cn = `awk '{ print \$1 }' $cf`;
    chomp $cn;
    push @cns, $cn if !$cnused{$cn}++;
    my $tc = "jq . $f";
    my $pjson = `jq . $f`;
    $rpt{$cn} .= qq[$id:
$pjson
$dashes
];
}

    
for my $k ( sort { $a cmp $b } @cns ) {
    print qq[$eqs
$k
$eqs
$rpt{$k}
];
}

    
