#! /usr/bin/perl -wT

# cgi setup
use strict;
use CGI;

# debug message
my $message = "";

# ordering
my %order_actions = (
    first => "[",
    before => "<",
    after => ">",
    "last" => "]",
);
my @order_actions = qw(first before after last);

# server id
my $server = $ENV{SERVER_NAME};

# hidden choices
my %hidden = map { ($_, 1) } qw(none timour breast micro);
if ($server eq "wigtop2.cshl.edu" || $server eq "wigtop1.cshl.edu") {
    %hidden = ("none", 1);
}
    
# view choices
my %choices = (
    chr => {
        choices => [ qw(all
        1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y all) ],
    },
    range => {
        choices => [ qw(common sample) ],
    },
    scale => {
        choices => [ qw(atan log linear) ],
    },
    color => {
        choices => [ qw(red blue) ],
    },    
    n_x => {
        default => "5",
        choices => [ qw(1 2 3 4 5 6 7 8 9 12) ],
    },
    data => {
        choices => [ qw(jason chd breast micro timour sonam) ], 
        group_names => { breast => "tumors",
                         chd => "families",
                         jason => "runs" },
    },
    group => {
        choices => [ qw(all USER) ],
    },
    selection => {
        choices => [ qw(all) ],
    },
    action => {
        choices => [ qw(none tag delete), @order_actions ],
    },
    move => {
        choices => [ qw(none) ],
    },
    order => {
        choices => [ qw(no yes) ],
    },
    sort => {
        choices => [ qw(user name nsam) ],
    },
    user => {
       choices => [ qw(none USER) ], 
    },
    interactive => {
       choices => [ qw(false true) ], 
    },
);
my @header_choices = qw(chr range scale color n_x order group data user);
my @query_choices = @header_choices;
my @choices = (@header_choices, qw(action move interactive));

our ($chr, $range, $scale, $color, $n_x, $order, $sort, $data, $group, $action, $move,
    $user, $interactive);
$sort = "user";

# special ids to distinguish tags from groups from samples
my $tag_id = 1;
my $normal_id = 2;
my $sample_id = 3;
my $group_id = 4;

# get choices made from user query
my $q = CGI->new;
sub get_choice {
    my $choice = $_[0];
    return if exists $choices{$choice}{choice};

    # set up default and valid choices 
    $choices{$choice}{default} = $choices{$choice}{choices}[0]
        if ! exists $choices{$choice}{default};
    $choices{$choice}{valid}{$_} ||= $normal_id
        for @{$choices{$choice}{choices}};

    # carefully get choice
    no locale;
    my $user_input = $q->param($choice) || "";
    $user_input =~ /^([\w.\#-]+)$/;
    $user_input = $1;
    $choices{$choice}{user_input} = $user_input || "";
    $user_input = $choices{$choice}{default}
    if ! $user_input || ! exists $choices{$choice}{valid}{$user_input};
    $choices{$choice}{choice} = $user_input;
    no strict "refs";
    $$choice = $user_input if $choice !~ /^\d/;
}

# get just dataset from query
get_choice("data");
if ($data eq "breast") {
    push @{$choices{sort}{choices}}, "size"; 
    push @{$choices{sort}{choices}}, "age"; 
}

# read in dataset information
open SAMPLES, "data/$data/samples.txt";
while (<SAMPLES>) {
    chomp;
    my ($group, $sample) = /(\w+) ([\w.\#-]+)/;
    if (! exists $choices{group}{valid}{$group}) {
        $choices{group}{valid}{$group} = $group_id;
        $choices{all}{valid}{$group} = $group_id;
        push @{$choices{move}{choices}}, $group;
        push @{$choices{all}{show}}, $group;
    }
    $choices{group}{valid}{$sample} = $sample_id;
    push @{$choices{move}{choices}}, $sample;
    push @{$choices{$group}{show}}, $sample;
    push @{$choices{$sample}{show}}, $sample;
}
close SAMPLES;

# read in existing user names
open USERS, "users/users.txt";
while (<USERS>) {
    chomp;
    $choices{user}{valid}{$_} = $tag_id;
}
close USERS;

# get new user name
get_choice("user");
if ($choices{user}{user_input} ne $choices{user}{choice} &&
    $choices{user}{user_input} =~ /^[a-zA-z]+$/) {
    $choices{user}{valid}{$choices{user}{user_input}} = $normal_id;
    $user = $choices{user}{choice} = $choices{user}{user_input};
    push @{$choices{user}{choices}}, $choices{user}{user_input};
    open USERS, ">>users/users.txt";
    print USERS "$user\n";
    close USERS;
} elsif ($choices{user}{choice} ne "none") {
    push @{$choices{user}{choices}}, $choices{user}{choice};
}
if ($user eq "zetterberg" || $user eq "wigler" || $user eq 'mskcc') {
    delete $hidden{breast};
    delete $hidden{micro};
}
if ($user eq "timour") {
    delete $hidden{timour};
}

# read in groups, treat almost as normal groups
my @groups = ();
if (open(GROUPS, "$data/groups.$user.txt")) {
    while (<GROUPS>) {
        chomp;
        $choices{group}{valid}{$_} = $tag_id;
        push @{$choices{group}{choices}}, $_;
        push @groups, $_;
    }
    close GROUPS;
}

# get choice from query
get_choice($_) for keys %choices;

# get and validate new group name if necessary
if ($choices{group}{user_input} ne $choices{group}{choice} &&
    $choices{group}{user_input} =~ /^[a-zA-z]\w*$/) {
    mkdir "$data/$choices{group}{user_input}"
        if ! -e "$data/$choices{group}{user_input}";
    $choices{group}{valid}{$choices{group}{user_input}} = $tag_id;
    $group = $choices{group}{choice} = $choices{group}{user_input};
    push @{$choices{group}{choices}}, $choices{group}{user_input};
    open GROUPS, ">>$data/groups.$user.txt";
    print GROUPS "$group\n";
    close GROUPS;
}

my $is_a_tag = $choices{group}{valid}{$group} == $tag_id;

# delete group if necessary
if ($is_a_tag && $action eq "delete") {
    unlink "$data/$group/order.$user.txt";
    @groups = grep { $_ ne $group } @groups;
    open GROUPS, ">$data/groups.$user.txt";
    print GROUPS "$_\n" for @groups;
    close GROUPS;
    @{$choices{group}{choices}} =
        grep { $_ ne $group } @{$choices{group}{choices}};
    $group = $choices{group}{choice} = "all";
    $action = "none";
    $is_a_tag = 0;
    
}

sub read_order {
    my $group = $_[0];
    my $order = "data/$data/$group/order/order.$user.txt";
    my @order = ();
    if (open ORDER, $order) {
        chomp(my $line = <ORDER>);
        @order = split /\s/, $line; 
        close ORDER;
    }
    return @order;
}

my %groups;
for my $group (@groups) {
    %{$groups{$group}} = map { ($_, 1) } read_order($group);
}

# set up show list for tagged group
if ($is_a_tag) {
    @{$choices{$group}{show}} = read_order($group);
} elsif ($group ne "all") {
    push @{$choices{group}{choices}}, $group;
}

# html header
print qq'content-type: text/html

<html>
<head>
<title>$data CN Viewer</title>
<style>
body { margin:0; padding:0; font-size:1.25em; width:100%; height:100%; }
a { text-decoration:none; color:#55F; font-weight:bold; }
div#header { text-align:center; background-color:#EEE; font-weight:bold; }
div#header p { margin:0px; padding:0.5em; border-bottom: 4px solid #B00; }
div#header span { color:#B00; }
div#body { margin:0.5em; }
div.thumb { float:left; }
div.thumb img { display:block; width:100%; }
div#copy { clear:both; text-align:right; padding:1em; padding-top:50em; }
span.chosen { font-weight:bold; color:#B00; }
h1, h2, h3, pre { clear:both; }
h2, h3 { padding-top: 1em; }
iframe { border:none; width:100%; height:100%;}
</style>
</head>
<body>
';
# filter:hue-rotate(60deg);

# report view choices
if (0) {
    print "<p>$#choices\n";
    for my $choice (keys %choices) {
        print "$choice=$choices{$choice}{choice}= " .
            "user=$choices{$choice}{user_input}= " .
            "default=$choices{$choice}{default}= ( ";
        my @options;
        push @options, "$_" .
            (exists $choices{$choice}{valid}{$_} ? "" : ":BAD")
            for (@{$choices{$choice}{choices}}, "bad");
        print join(" | ", @options) . " )<br />\n";
    }
    print "</p>\n";
}

# base query to use for links
my $query = "?" . join("&", "view=cn", map { 
    $choices{$_}{choice} eq $choices{$_}{default} ?
        () : "$_=$choices{$_}{choice}" } @query_choices);
my $default_form = join("\n", map { 
    $choices{$_}{choice} eq $choices{$_}{default} ||
        $_ eq "group" ?
        () : qq'<input type="hidden" name="$_" value="$choices{$_}{choice}">' }
                        @query_choices);
my $form = qq'';

# non breakable spaces
my $spc = "&nbsp;";
my $spc2 = "${spc}${spc}";
my $dash = "${spc2}-\n${spc}";

# choosing view options
for my $choice (@header_choices) {
    (my $squery = $query) =~ s/&$choice=[^&]+//;
    $squery =~ s/&group=[^&]+// if $choice eq "data" || $choice eq "user";
    # my $sep = $choice eq "chr" ? "${spc}| " : "${spc}|${spc}";
    my $sep = "$spc|$spc";
    $choices{$choice}{select} = join(
        $sep, map {
            $_ eq $choices{$choice}{choice} ? qq'<span>$_</span>' :
                ( $_ eq "USER" ?
                  qq'<input style="width:5em;" type="text" name="$choice">' :
                  qq'<a href="$squery&$choice=$_">$_</a>' ) }
        map({exists $hidden{$_} ? () : ($_)} @{$choices{$choice}{choices}}));
}
my $choices = ${spc2} . join(
    $dash, map { scalar(@{$choices{$_}{choices}}) > 1 ?
                     "<b>$_</b>:${spc2}($choices{$_}{select})" :
                     () } @header_choices) .
    qq'$dash<a href="./?data=$data">set${spc}defaults</a>';
$choices .= qq'$dash<a href="breast.xls">breast.xls</a>' if $data eq 'breast';

print qq'<form style="display:inline!important;" action="./" method="get">
<div id="header">
<p>
$choices
$default_form
<input type="submit" value="submit" style="visibility:hidden;">
</p>
</div>
<div id="body">
';

# title to show for each profile
sub profile_title {
    my $name = $_[0];
    (my $sname = $name) =~ s/.*?\.//;
    my $num = ($group eq 'all' || $is_a_tag) ? "$spc:$spc$_[1]" : "";
    my $result = $spc;
    for my $action (@order_actions) {
        if ($action eq 'after') {
            my $title = qq'<b>$sname</b>$num';
            $title = qq'<a href="$query&move=$name&action=tag">$title</a>'
                if $is_a_tag && $order eq "yes";
            $result .= qq'$title$spc';
        }
        next if $_[3] == 0 || $_[2] == 1 || $order eq "no" || $group =~ /\./;
        $result .= qq'<a href="$query&move=$name&action=$action">' .
            qq'$order_actions{$action}</a>$spc';
    }
    return $result;
}

sub byNumInName {
    (my $a = $_[0]) =~ s/\D+//g;
    (my $b = $_[1]) =~ s/\D+//g;
    return $a <=> $b;
}

# ordering of view
sub order {
    my $group = $_[0];
    my @show = sort { byNumInName($a, $b); } @{$choices{$group}{show}};
    my @order;
    if ($sort eq "name") {
        @order = @show;
    } elsif ($sort eq "nsam" && ($group eq "all" || $is_a_tag)) {
        @order = sort {
            my $na = scalar @{$choices{$a}{show}};
            my $nb = scalar @{$choices{$b}{show}};
            if ($na == $nb) {
                byNumInName($a, $b);
            } else {
                $nb <=> $na;
            }
        } @show;
    } elsif (($sort eq "size" || $sort eq "age") &&
             ($group eq "all" || $is_a_tag)) {
        my %meta;
        open META, "$data/meta.tsv";
        chomp(my $header = <META>);
        # $message .= "header: $header\n";
        my @header = split /\t/, $header;
        my $key = $sort eq "size" ? "Size (mm)" : "Age";
        while (<META>) {
            chomp(my $line = $_);
            my @line = split /\t/, $line;
            $message .= "$group $line[0]\n";
            next if ! exists $choices{$group}{valid}{$line[0]};
            $meta{$line[0]}{$header[$_]} = $line[$_] for 0 .. $#line;
            $message .= "line: $#line : $line[0] : $meta{$line[0]}{$key} : $line\n";
        }
        @order = sort { $meta{$b}{$key} <=> $meta{$a}{$key} } keys %meta;
    } elsif ($sort eq "user") {
        @order = read_order($group);
    }

    my %order = map { ($_, 1) } @order;
    for (@show) {
        push @order, $_ if ! exists $order{$_};
    }
    # $message .= scalar(@show) . " " . scalar(@order) . "\n"; 
    return [ @order ];
}
my @order = @{order($group)};
if ($action ne "none" && $move ne "none") {
    @order = ($move, grep(!/^$move$/, @order)) if $action eq "first";
    @order = (grep(!/^$move$/, @order), $move) if $action eq "last";
    if ($action eq "before" || $action eq "after") {
        my ($index) = grep { $order[$_] eq $move } (0 .. @order - 1); 
        if (defined $index) {
            if ($action eq "before") {
                ($order[$index - 1], $order[$index]) =
                    ($order[$index], $order[$index - 1]) if $index;
            } elsif ($action eq "after" && $index < @order - 1) {
                ($order[$index], $order[$index + 1]) =
                    ($order[$index + 1], $order[$index]);
            }
        }
    }
    if ($is_a_tag && $action eq "tag") {
        my %ingroup = map { ($_, 1) } @order;
        if (exists $ingroup{$move}) {
            @order = grep { $_ ne $move } @order;
        } else {
            push @order, $move;
        }

    }
    my $order = "data/$data/$group/order/order.$user.txt";
    open ORDER, ">$order";
    print ORDER join(" ", @order) . "\n";
    close ORDER;
}

# breast tumor metadata, pre-prepared for html
if ($group ne "all" && $data eq "breast") {
    (my $real_group = $group) =~ s/\..+$//;
    if (open META, "$data/$real_group/meta.txt") {
        my $meta = <META>;
        print $meta;
        close META;
    }
}

if ($order eq "yes") {
    print qq'<h3>Order profiles by clicking on the symbols next to profile title</h3>'
}

# profiles not in group, for tags only
my @notgroup = ();
if ($is_a_tag) {
    my %ingroup = map { ($_, 1) } @order;
    for my $sample (@{$choices{all}{show}}) {
        push @notgroup, $sample if ! exists $ingroup{$sample};
    }
    if ($order eq "yes") {
        print qq'<h3>Move profile in and out of group $group by clicking its title</h3>';
    }
    print qq'<h3>Remove the group $group by clicking <a href="$query&action=delete">here</a></h3>'
        if scalar(@order) == 0;
}

my @show = (@order, @notgroup);

my $n_samples = scalar @show;
my $width = 100 / $n_x;
if ($n_samples <= $n_x) {
    if ($n_samples == 1) {
        $width = 100;
    } elsif ($n_samples > 1 && $n_samples <= 4) {
        $width = 50;
    } elsif ($n_samples > 4 && $n_samples <= 9) {
        $width = 33;
    } elsif ($n_samples >= 10) {
        $width = 20;
    }
}

# show profiles
my $commands = "";
(my $ngquery = $query) =~ s/&group=[^&]+//;
$chr = "chr$chr" if $chr ne "all";
my $in_group = 1;
for my $show (@show) {
    next if $show eq 'all';
    if (scalar(@notgroup) > 0 && $show eq $notgroup[0]) {
        last if $order ne "yes";
        print qq'<h3>Profiles not in group $group:</h3>\n';
        $in_group = 0;
    }
    my $is_a_sample = ($show =~ /\./) ? 1 : 0;
    my @sorder = $is_a_sample ? ($show) : @{order $show};
    my $show_sample = $sorder[0];
    (my $actual_group = $show_sample) =~ s/\..+$//;
    my $id = "${show_sample}_${chr}_${range}_${scale}_${color}";
    my $base = "./data/$data/$actual_group/$show_sample/thumbs/$id";

    my $font = sqrt($width * 0.02);
    $font = 0.75 if $font < 0.75;

    my $cbox = "";
    if (0 && $is_a_tag) {
        $cbox = qq'<input type="checkbox" name="' .
            ($in_group ? "remove" : "add") . qq'" value="$show">';
    }

    my $img = qq'<img src="$base/$id.png" title="$show" alt="$show" />';
    if ($interactive eq "true") {
        $img = qq'<iframe src="/ggraph_web/?data=cn/data/$data/$actual_group/$group/$group.txt&abspos_low=0&abspos_high=3095498750&cn_low=0&cn_high=1">';
    }
    $img = qq'<a href="$ngquery&group=$show">$img</a>' if $show ne $group;
    print qq'<div class="thumb" style="width:$width%;"><center style="font-size:${font}em;">' . profile_title($show, scalar @sorder, scalar @order, $in_group) . $cbox . qq'</center>$img</div>\n';     
    
    if ($choices{group}{valid}{$group} == $sample_id) {
        print qq'<h3>View this sample <a href="$query&interactive=true">interactively</a></h3>';
        print qq'<h3>Binned data file:
            <a href="data/$data/$actual_group/$show_sample/$show_sample.txt">$show_sample.txt</a></h3>\n';

        
    }

    # write command to run for ggraph
    if ($group ne "all" && open COMMAND, "$base/command.txt") {
        chomp(my $command = <COMMAND>);
        close COMMAND;
        $commands .= "\n$command";
    }
}

if (0) {
    my $sample= "moooo";
    print qq'<h2>Binned data file:</h2>
<p><a href="data/$data/$group/$sample/$sample.txt">$sample.txt</a></p>\n';
}


if (0) {
if (scalar @groups &&
    $choices{group}{valid}{$group} == $group_id ||
    $choices{group}{valid}{$group} == $sample_id) {
    (my $actual_group = $group) =~ s/\..+$//;
    print qq'<h3>Change group membership for $actual_group.';
    my @in;
    my @out;
    for my $group (@groups) {
        if (exists $groups{$group}{$actual_group}) {
            push @in, $group;
        } else {
            push @out, $group;
        }
    }
    print qq' In:';
    if (scalar @in) {
        print qq' <a href="$ngquery&group=$_&move=$actual_group&action=tag">$_</a>' for @in;
    } else {
        print " none";
    }
    print qq' Out:';
    if (scalar @out) {
        print qq' <a href="$ngquery&group=$_&move=$actual_group&action=tag">$_</a>' for @out;
    } else {
        print " none";
    }
    print qq'</h3>\n';
}
}

# ggraph view commands
if ($commands ne "" && $ENV{HTTP_HOST} =~ /wigclust/) {
    print qq'<h2>G-Graph view commands (only accessible from wigclust):</h2><pre>$commands</pre>\n';
}

# finish page
print qq'
<pre>
$message
</pre>
<p>
<br /><br />
<br /><br />
<br /><br />
<br /><br />
<br /><br />
</p>
</div>
</form>

<div id="copy">
<p>
&copy; 2018-2020
<a href="http://drpa.us/">Peter Andrews</a>
@
<a href="http://cshl.edu/">CSHL</a>
in the
<a href="http://cshl.edu/research/faculty-staff/michael-wigler/">Wigler Lab</a>
</p>
</div>
</body>
</html>
';

exit 0;

__END__


