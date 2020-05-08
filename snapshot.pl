#!/usr/bin/perl

use vDOM::NeatHTML;
use Bio::DB::Fasta;
use Bio::SeqIO;
use strict;

$| = 1;
my $STEP   = 10;
my $ROOT   = '/data/human_ref/hg19_sanger/';
my $SAMTOOLS_ROOT = '';
my $WINDOW = 100;
my $AIR    = 5;   # Minimum distance between two seqs in the same line
my $MAX_QUAL = 40;
my $LARGE_STEP = 10;
my $QUAL_SHIFT = 33;
# Load reference sequence


my $tpath = $ROOT."hg19.fa";
my %db = ();
my $tobj = tie %db,'Bio::DB::Fasta', $tpath;
my @gchr = ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
            '11', '12', '13', '14', '15', ',16', '17', '18', '19', '20',
            '21', '22', 'X', 'Y', 'MT');

my @qcolors = load_colors($MAX_QUAL);


my $bam = shift;
my $plist = shift;
my $window = shift;
$WINDOW = $window if ($window);
my @t = split(/[\/\\]/, $0);
my $pname = pop @t;
die("Use: $pname bam_file coord_file [window_size]\n")
  unless ($bam && $plist && -e $bam && -e $plist);
  # plist
  my @plist = ();
  my @t = split /[\/\\]/, $plist;
  my $bout = pop @t;
  $bout =~ s/\..+//;
  my @t2 = split /[\/\\]/, $bam;
  my $bout2 = pop @t2;
  $bout2 =~ s/\..+//;
  $bout .= "\_$bout2";
  open (IN, $plist) or die ("Could not open $plist: $!\n");
    while (<IN>){
      chomp;
      my ($chr, $pos) = split /[\t\s,]+/;
      push @plist, [$chr, $pos] if ($chr && $pos);
    }
  close IN;
  foreach my $pair (@plist){
    my ($chr, $pos) = @$pair;
    print "Reading chromosome $chr, position $pos...\n";
    $chr =~ s/chr//gi;
    my $uparams = {};
    $uparams->{'bam'}    = $bam;
    $uparams->{'chrID'}  = $chr;
    $uparams->{'offset'} = $pos;
    my $fname = $bout."\_$uparams->{'chrID'}\_$uparams->{'offset'}\.html";
    open (OUT, ">$fname");
    my $current_bam = 0;
    my $wpage = vDOM::NeatHTML->new({'title' => $bam});
    $wpage->add_style({
      'name' => '.hcont',
      'padding' => '0',
      'margin' => '0',
      'display' => 'inline',
      'cursor' => 'pointer'
    });
    $wpage->add_style({
      'name' => 'span',
      'padding' => '1px'
    });
    $wpage->add_style({
      'name' => '.hit',
      'color' => '#FFFFFF'
    });
    $wpage->add_style({
      'name' => '.highlighted',
      'padding' => '0',
      'border-right' => '1px solid #000000',
      'border-left' => '1px solid #000000'
    });
    $wpage->add_style({
      'name' => '.genomic',
      'cursor' => 'pointer',
      'background-color' => '#44ff00'
    });
    $wpage->add_style({
      'name' => '.result',
      'float' => 'left',
      'clear' => 'both',
      'padding' => '1em',
      'position' => 'relative',
      'font-family' => 'Courier New',
      'font-size' => 'x-small',
      'font-weight' => 'bold'
    });
    $wpage->add_style({
      'name' => '.result_frame',
  		'float' => 'left',
  		'clear' => 'both',
  		'padding' => '1em',
  		'position' => 'relative',
  		'border' => 'thick groove #000000',
  		'font-family' => 'Courier New',
  		'font-size' => 'x-small',
  		'font-weight' => 'bold'
    });
    $wpage->add_style({
      'name' => '.reads',
  		'float' => 'left',
  		'clear' => 'both',
  		'padding' => '1em',
  		'position' => 'relative',
  		'border' => 'thick groove #000000'
    });
    $wpage->add_style({
    	'name' => '.pop',
      'background-color' => '#CCCCCC',
      'position' => 'absolute',
      'overflow' => 'visible',
      'z-index' => '999',
      'display' => 'inline',
      'margin' => '0.2em',
      'border' => 'outset'
    });
    $wpage->add_style({
    'name' => '.popbar',
      'font-family' => 'Arial',
      'font-size' => 'x-small',
      'background-color' => '#AAAAAA',
      'border-bottom' => 'inset',
      'position' => 'static',
      'height' => 'auto'
    });
    $wpage->add_style({
    	'name' => '.close',
      'font-family' => 'Arial',
      'cursor' => 'pointer',
      'font-size' => 'small',
      'position' => 'absolute',
      'background-color' => '#555555',
      'color' => '#ffffff',
      'top' => '0',
      'right' => '0',
      'vertical-align' => 'middle',
      'border-style' => 'outset',
      'border' => 'thin'
    });
    $wpage->add_style({
    'name' => '.close:hover',
    'background-color' => '#994444'
    });
    $wpage->add_style({
    	'name' => '.popup_table',
      'border-style' => 'double',
      'width' => 'auto'
    });
    $wpage->add_style({
    	'name' => '.popup_table a',
      'color' => '#0000aa'
    });
    $wpage->add_style({
    	'name' => '.qgap',
    	'color' => '#000000',
    	'background-color' => '#aaffaa'
    });
    $wpage->add_style({
    	'name' => '.insertion',
    	'color' => '#000000',
    	'background-color' => '#ffaaaa'
    });
    $wpage->add_style({
    	'name' => '.popup_table td',
    	'font-size' => 'x-small',
    	'border' => 'thin solid #555555',
    	'padding' => '0.5em',
    	'color' => '#000000'
    });
    $wpage->add_style({
      'name' => '.sam',
      'float' => 'left',
      'padding' => '1em',
      'clear' => 'both',
      'margin-top' => '1em'
    });
    $wpage->add_script('
    var currentPopup = \'\';

    function th(col){
      var id = "l0p"+col;
      toggleHigh(id);
      for (var i = 0; i < nrows; i++){
        var id = "r"+i+"c"+col;
        toggleHigh(id);
      }
    }
    function toggleHigh(id){
      var p = new RegExp(/\s*highlighted/);
      var el = document.getElementById(id);
      if (el == undefined){return 0};
      var cl = el.className;
      if (p.test(cl)==true){
        var cl2 = cl.replace(p, "");
        el.className = cl2;
      }
      else{
        el.className = cl+" highlighted";
      }
    }

    function toggleInfo(myID){
    	if (currentPopup == myID){
    		destroy_tpopup();
    	}
    	else{
    		destroy_tpopup();
      	tpopup(myID, "Hit info", lcol, info[myID]);
      }
    }

    function spopup_crt(myID, code){
    	var el = document.getElementById(myID);
    	el.innerHTML += code;
    	alert(code);
    }


    function tpopup(myID, caption, ldata, rdata){
    	/*
    	Shows a popup table with two columns (from ldata and rdata arrays) which stems from myID.
    	For proper behaviour, this must be used with popup.css to provide non-disruptive styles.
    	The "myID" element must be enclosed in a parent element, where the popup will be created.
    	*/
    	var el = document.getElementById(myID);
      var l1 = \'\';
      for (var i = 0; i < ldata.length; i++){
      	l1 += \'<tr><td>\'+ldata[i]+\'</td><td>\'+rdata[i]+\'</td></tr>\';
      }
      var header = "<div class=\'popbar\'>"+ caption +
         "<span class=\'close\' onclick=\'destroy_tpopup()\'>X</span></div>";
      var code = "<div class=\'pop\' id=\'current_popup\'>"+header+\'<table class="popup_table">\'+l1+\'</table>\';
      el.parentNode.innerHTML += code;
      currentPopup = myID;
    }
    function destroy_tpopup(){
    	var el = document.getElementById(\'current_popup\');
      if (el){
      	el.parentNode.removeChild(el);
      }
      currentPopup = \'\';
    }
    ');
    set_result($wpage, $uparams);
    unless (exists($uparams->{'pos'}) && $uparams->{'pos'}){
      $uparams->{'pos'} = $uparams->{'offset'};
    }
    my $temp = {};
    $temp = {'onload'=>"th($uparams->{'pos'})"} if (exists($uparams->{'pos'}) && $uparams->{'pos'});
    my $html = $wpage->get_html($temp);
    #$html =~ s/\n\s+/\n/g;
    my $cpre = '<pre class="sam">.*?<\/pre>';
    my ($pre) = $html =~ /($cpre)/sm;
    if ($pre){
      $pre  =~ s/\n\s+/\n/g;
      $html =~ s/$cpre/$pre/sm;
    }
    print OUT $html;
    close OUT;
  }


sub set_result{
  my $wpage   = shift;
  my $uparams = shift;
  return '' unless (exists($uparams->{'chrID'}));
  my $offset = $uparams->{'offset'};
  if ($offset > $WINDOW){
    $offset -= $WINDOW >> 1;
  }
  my $end = $offset + $WINDOW;
  if (exists($uparams->{'pos'})){
    if ($uparams->{'pos'} < $offset || $uparams->{'pos'} > $end){
      $uparams->{'pos'} = '';
    }
  }
  my $current_bam = 0;
  if (exists($uparams->{'bam'})){
    $current_bam = $uparams->{'bam'};
  }
  #Load quality colors
  for (my $i = 0; $i < @qcolors; $i++){
    $wpage->add_style({
      'name' => ".q$i",
      'background-color' => $qcolors[$i]
    });
  }
  ##Result
  $wpage->add({
    'tag' => 'div',
    'class' => 'result_frame'
  });
  my $sam = '';
  $wpage->get_into();
    # Get genomic seq
    my $chr = $uparams->{'chrID'};
    my $ens_chr = $chr;
    $ens_chr =~ s/chr//i;
    $ens_chr = 'MT' if ($ens_chr =~ /^m/i);
    my $gstr = "$ens_chr\:$offset\,$end";
    my $gseq = {};
    my $genomic = genomic($db{$gstr}, $offset, $gseq);
    my @hits = ();
    my $lines = [];
    if ($offset > 0){
      my $start = $offset;
      my $rcoord = $offset + $WINDOW;
      $lines = get_reads($current_bam, $chr, $start, $rcoord);
      $wpage->add_script("var info = new Array();\n");
      my $lcol = join(', ',
        '"Id"',
        '"Chr"',
        '"Pos"',
        '"MapQ"',
        '"Cigar"',
        '"Mate_chr"',
        '"Mate_pos"',
        '"Ins_size"',
        '"Sequence"'
      );
      $wpage->add_script("var lcol = [$lcol];\n");
      my $c = 0;
      foreach my $line (@$lines){
        my $temp = $line;
        my @t = split(/\t/, $temp);
        cchomp(\$t[0]);
        my $l = join('", "',
          '"'.
          $t[0],
          $t[2],
          $t[3],
          $t[4],
          $t[5],
          $t[6],
          $t[7],
          $t[8],
          $t[9]
          .'"'
        );
        my $id = "seq$c";
        $wpage->add_script("info[\"$id\"] = [$l];\n");
        $c++;
        $sam .= "$temp\n";
      }
      my ($screen, $insertions) = set_screen($lines, $offset, $end, $gseq);
      # Insertions
      my $last_row = scalar @$screen;
      foreach my $ins (keys %{$insertions->{'pos'}}){
        my $igtext = '';
        # Genomic
        my $gid = "l0p$ins";
        for (0 .. $insertions->{'pos'}{$ins}-1){
          my $iid = $gid."i$_";
          $igtext .= "<span id=\"$iid\">-</span>";
        }
        my $gtext = "<span.*?id=\"$gid\"[^<>]*?>[^<>]*?</span>";
        $genomic =~ s/($gtext)/$1$igtext/;
        # Reads
        my $rid = "r\\d+c$ins";
        my $rtext = "<span[^<>]*?id=\"$rid\"[^<>]*?>[^<>]*?</span>";
        for (my $i = 0; $i < scalar @$screen; $i++){
          my $irtext = '';
          for (0 .. $insertions->{'pos'}{$ins}-1){
            my $iid = "r$i"."c$ins"."i$_";
            $irtext .= "<span id=\"$iid\" class=\"insertion\">-</span>";
          }
          $screen->[$i]{'text'} =~ s/($rtext)/$1$irtext/;
        }
      }
      foreach my $ins (@{$insertions->{'ins'}}){
        my $iid = 'r'.$ins->{'row'}.'c'.$ins->{'pos'}.'i'.$ins->{'n'};
        my $itext = "<span[^<>]*?id=\"$iid\"[^<>]*?>-</span>";
        my $newtext = $ins->{'code'};
        $screen->[$ins->{'row'}]{'text'} =~ s/$itext/$newtext/;
      }
      #####################
      # Add genomic
      $wpage->add({
        'tag' => 'div',
        'class' => 'genomic',
        'text' => $genomic
      });
      # Add reads
      foreach my $line (@$screen){
        $wpage->add({
          'tag'   => 'div',
          'class' => 'hit',
          'text'  => $line->{'text'}
        });
      }
      $wpage->add_script("
        var nrows = $last_row;".'
        var current = 0;
      ');
    }
  $wpage->get_out();
  if ($sam){
    $wpage->add({
      'tag' => 'pre',
      'class' => 'sam',
      'text' => $sam
    });
  }
}

sub get_reads{
  my $file = shift;
  my $chr  = shift;
  my $from = shift;
  my $to   = shift;
  my $command = $SAMTOOLS_ROOT."samtools view $file $chr:$from\-$to";
  my $lines = `$command`;
	warn ("$command\n");
  my @result = split(/\n/, $lines);
  return \@result;
}



sub set_screen{
  my $hits = shift;
  my $from = shift;
  my $to   = shift;
  my $gseq = shift;
  my $screen = [];
  my $temp = blank_line(0, $from, $to);
	push @$screen, $temp;
  my $c = 0;
  my $insertions = {};  # To keep track of insertions
	my $cov = [];         # To keep track of coverage
  foreach my $sam (@$hits){
    my @fields = split(/\t/, $sam);
    if ($fields[3] <= $to){
      my $cline = check_line($screen, $fields[3], $from, $to);
      my $id = "seq$c";
      add_hit(\@fields, $screen, $cline, $from, $to, $gseq, $id, $insertions);
      $c++;
    }
  }
  return ($screen, $insertions);
}

sub add_hit{
  my $hit    = shift;
  my $screen = shift;
  my $line   = shift;
  my $from   = shift;
  my $to     = shift;
  my $gseq   = shift;
  my $id     = shift;
  my $insertions = shift;
  my $pos = 0;
	my $seq = get_bam_seq($hit, $line, $from, $to, $gseq, $insertions);
  if (scalar @$seq > 1){
    ### First hit
    my $name = $hit->[0];
    my $fbase = shift @$seq;
    my $fb    = $fbase->{'text'};
    $pos  = $fbase->{'pos'};
    $screen->[$line]{'text'} =~ s/(<span[^<>]*?id=[^<>]*?c$pos[^<>]*?>\&nbsp\;<\/span>)/<div class=\"hcont\"><div class=\"hcont\" id=\"$id\" onclick=\"toggleInfo\(\'$id\'\)\" title=\"$name\">$fb/;
    ### Last hit
    my $lbase = pop @$seq;
    my $lb    = $lbase->{'text'};
    $pos  = $lbase->{'pos'};
    $screen->[$line]{'text'} =~ s/(<span[^<>]*?id=[^<>]*?c$pos[^<>]*?>\&nbsp\;<\/span>)/$lb\<\/div><\/div>/;
    ### Rest
    while (scalar @$seq){
      my $carryover = '';
      my $base = shift @$seq;
      my $b = $base->{'text'};
      $pos  = $base->{'pos'};
      my $icode = '';
      while (scalar @$seq && $seq->[0]{'type'} eq 'insertion'){
        my $ins = shift @$seq;
      }
      $screen->[$line]{'text'} =~ s/(<span[^<>]*?id=[^<>]*?c$pos[^<>]*?>\&nbsp\;<\/span>)/$b$icode/;
    }
  }
  elsif (scalar @$seq == 1){  # Only one base
    my $name = $hit->[0];
    my $fbase = shift @$seq;
    my $fb    = $fbase->{'text'};
    $pos  = $fbase->{'pos'};
    $screen->[$line]{'text'} =~ s/(<span[^<>]*?id=[^<>]*?c$pos[^<>]*?>\&nbsp\;<\/span>)/<div class=\"hcont\"><div class=\"hcont\" id=\"$id\" onclick=\"toggleInfo\(\'$id\'\)\" title=\"$name\">$fb<\/div><\/div>/;
  }
  $screen->[$line]{'pos'} = $pos;
}

sub check_line{
  my $screen = shift;
  my $pos    = shift;
  my $from   = shift;
  my $to     = shift;
  for (my $i = 0; $i < @$screen; $i++){
    if ($pos >= $screen->[$i]{'pos'} + $AIR){
      return $i;
    }
  }
  my $temp = blank_line(scalar @$screen, $from, $to);
  push @$screen, $temp;
  return scalar @$screen - 1;
}

sub blank_line{
  my $row  = shift;
  my $from = shift;
  my $to   = shift;
  my $result = {};
  my $text = '';
  for (my $i = $from; $i <= $to; $i++){
    $text .= "<span id=\"r$row"."c"."$i\">&nbsp;<\/span>";
  }
  $result->{'text'} = $text;
  $result->{'pos'} = 0;
  return $result;
}

sub sp{
  my $n = shift;
  return '' if ($n <= 0);
  my $result = '';
  for (1..$n){
    $result .= '<span>&nbsp;</span>';
  }
  return $result;
}

sub get_bam_seq{
  my $hit  = shift;
  my $n    = shift;
  my $from = shift;
  my $to   = shift;
  my $gseq = shift;
  my $insertions = shift;
  my $pos = $hit->[3];
  my $s   = $hit->[9];
  my $q   = $hit->[10];
  my $result = [];
  my @seq = split(//, $s);
  my @q   = split(//, $q);
  # Change sequence according to CIGAR
  my $cigar = $hit->[5];
  my @cig = $cigar =~ /(\d+\w)/g;
  foreach my $code (@cig){
    my ($ncode, $lcode) = $code =~ /(\d+)(\w)/;
    $lcode = "M" unless ($lcode =~ /[MIDNSHP]/);
    my $temp = {};
    #if ($pos + $ncode >= $from && $pos <= $to){
      if ($lcode eq 'M'){
        for (1..$ncode){
          my $b = shift @seq;
          my $qual = shift @q;
          if ($pos >= $from && $pos <= $to){
            $qual = ord($qual) - $QUAL_SHIFT;
            $qual = $MAX_QUAL if ($qual > $MAX_QUAL);
            my $g = $gseq->{$pos};
            if (uc($g) eq uc($b)){
              if (check_flag($hit->[1], 4)){
                $b = ',';
              }
              else{
                $b = '.';
              }
            }
            $temp = {
              'pos' => $pos,
              'type' => 'normal',
              'text' => "<span class=\"q$qual\" id=\"r$n"."c$pos\">".$b."</span>"
            };
            push @$result, $temp;
          }
          $pos++;
        }
      }
      elsif ($lcode eq 'D' || $lcode eq 'N'){  #Gaps
        for (1..$ncode){
          if ($pos >= $from && $pos <= $to){
            my $b = '-';
            my $qual = 'gap';
            $temp = {
              'pos' => $pos,
              'type' => 'normal',
              'text' => "<span class=\"q$qual\" id=\"r$n"."c$pos\">".$b."</span>"
            };
            push @$result, $temp;
          }
          $pos++;
        }
      }
      elsif ($lcode eq 'S'){  #Soft clipping
        for (1..$ncode){
          shift @seq;
          shift @q;
          #$pos++;
        }
      }
      elsif ($lcode eq 'I'){  #Insertion
				$pos--;
        unless (exists($insertions->{'pos'}{$pos}) && $insertions->{'pos'}{$pos} >= $ncode){
          $insertions->{'pos'}{$pos} = $ncode;
        }
        my $icode = '';
        for (0..$ncode-1){
          my $b = shift @seq;
          my $qual = shift @q;
          $qual = ord($qual) - $QUAL_SHIFT;
          $qual = $MAX_QUAL if ($qual > $MAX_QUAL);
          $icode = "<span class=\"q$qual\" id=\"r$n"."c$pos"."i$_\">".$b."</span>";
          my $temp = {
            'row'  => $n,
            'pos'  => $pos,
            'n'    => $_,
            'code' => $icode
          };
          push @{$insertions->{'ins'}}, $temp;
        }
				$pos++;
      }
    #}
  }
  return $result;
}

sub genomic{
  my $seq   = shift;
  my $start = shift;
  my $gseq = shift;
  my @letters = split(//, $seq);
  my $result = '';
  my $counter = 0;
  foreach my $l (@letters){
    my $pos = $start + $counter;
    $gseq->{$pos} = uc($l);
    $result .= "<span id=\"l0p$pos\" title=\"pos:$pos\" onclick=\"th($pos)\">$l</span>";
    $counter++;
  }
  return $result;
}

sub check_flag{
  my $flag = shift;
  my $pos  = shift;
  my $mask = 1 << $pos;
  my $r = $flag & $mask;
  my $result = 0;
  if ($r){
    $result = 1;
  }
  return $result;
}

sub load_colors{
  my $steps = shift;
  my @result = ();
  for (my $i = 0; $i < 256; $i = $i + 255/$steps){
    my $reach = sprintf("%02X", 255 - $i);
    my $geach = '00';
    my $beach = sprintf("%02X", $i);
    my $each = "\#$reach$geach$beach";
    push @result, $each;
  }
  return @result;
}


sub mynum{
  my $n = shift;
  my ($result) = $n =~ /(\-*[\.\e\d]+)/;
  return $result;
}

sub debug{
  my $text = shift;
  open (OUT, ">debug.txt");
    print OUT "$text\n";
  close OUT;
}

sub getFiles{
	#Returns an array with the files in a directory
	#Args: DirPath, [filter]
	#if present, filter must be a regex

	my $dpath  = $_[0] ? $_[0] : ".";
	my $filter = $_[1] ? $_[1] : "" ;
	my @temp = ();
	opendir DIR, $dpath or die "unable to open $dpath\n";
	while	(my $name = readdir(DIR)){
		if (($name =~ /$filter/)){
			push @temp, $name;
		}
	}
	closedir DIR;
	my @result = sort{$a cmp $b} @temp;
	return @result;
}

sub set_arrows{
  my $wpage = shift;
  my $uparams = shift;
  $wpage->add({
    'tag' => 'div',
    'id'  => 'arrows'
  });
  $wpage->get_into();
  {
    #################
    $wpage->add({
      'tag' => 'form',
      'action' => 'tview.cgi',
      'method' => 'get'
    });
    $wpage->get_into();
    {
      arrow($wpage, $uparams, $uparams->{'offset'} - $LARGE_STEP, '&lt;&lt;');
    }
    $wpage->get_out();
    ##############
    #################
    $wpage->add({
      'tag' => 'form',
      'action' => 'tview.cgi',
      'method' => 'get'
    });
    $wpage->get_into();
    {
      arrow($wpage, $uparams, $uparams->{'offset'} - 1, '&lt;');
    }
    $wpage->get_out();
    ##############
    #################
    $wpage->add({
      'tag' => 'form',
      'action' => 'tview.cgi',
      'method' => 'get'
    });
    $wpage->get_into();
    {
      arrow($wpage, $uparams, $uparams->{'offset'} + 1, '&gt;');
    }
    $wpage->get_out();
    ##############
    #################
    $wpage->add({
      'tag' => 'form',
      'action' => 'tview.cgi',
      'method' => 'get'
    });
    $wpage->get_into();
    {
      arrow($wpage, $uparams, $uparams->{'offset'} + $LARGE_STEP, '&gt;&gt;');
    }
    $wpage->get_out();
    ##############
  }
  $wpage->get_out();
}

sub arrow{
  my $wpage    = shift;
  my $uparams  = shift;
  my $move     = shift;
  my $caption  = shift;
  my $cpos = $uparams->{'offset'};
  if (exists($uparams->{'pos'}) && $uparams->{'pos'}){
    $cpos = $uparams->{'pos'};
  }
  $wpage->add({
    'tag' => 'div',
    'class' => 'hcont'
  });
  $wpage->get_into();
  {
    $wpage->add({
      'tag' => 'input',
      'type' => 'hidden',
      'name' => 'bam',
      'value' => $uparams->{'bam'}
    });
    $wpage->add({
      'tag' => 'input',
      'type' => 'hidden',
      'name' => 'chrID',
      'value' => $uparams->{'chrID'}
    });
    $wpage->add({
      'tag' => 'input',
      'type' => 'hidden',
      'name' => 'offset',
      'value' => $move
    });
    $wpage->add({
      'tag' => 'input',
      'type' => 'hidden',
      'name' => 'pos',
      'value' => $cpos
    });
    $wpage->add({
      'tag' => 'input',
      'type' => 'submit',
      'value' => $caption
    });
  }
  $wpage->get_out();
}

sub cchomp{
  my $str = shift;
  $$str =~ s/\00//g;
}

sub header{

  # Returns the header of the BAM file into a hash:
  # magic: four BAM-specific bytes.
  # l_text: length of SAM header
  # text: SAM header
  # n_ref: number of reference seqs. (usually chromosomes)
  # reference_info: (array with n_ref elements)
  #    l_name: length of the ref. name plus 1
  #    name: reference name
  #    l_ref: length of the reference sequence

  my $z = shift;
  my $MAGIC = 21840194;
  my $result = {};
  $z->seek(0, 0);
  my $magic = $z->get_int32();
  unless ($magic == $MAGIC){
    warn("Incorrect magic number!\n");
    return 0;
  }
  my $l_text = $z->get_int32();
  $result->{'l_text'} = $l_text;
  $result->{'text'}  = $z->get_chars($l_text);
  $result->{'n_ref'} = $z->get_int32();
  for (0..$result->{'n_ref'}-1){
    my $temp = {};
    my $l_name = $z->get_int32();
    $temp->{'name'} = $z->get_chars($l_name);
    $temp->{'l_ref'} = $z->get_int32();
    push @{$result->{'reference_info'}}, $temp;
    $z->{'contig_translator'}{$temp->{'name'}} = $_;
  }
  $z->{'header'} = $result;
  $z->{'current_contig'} = 0;
  return 1;
}

sub get_int32{

  # Gets a little-endian int32 from an IO::Uncompress::Gunzip object
  # @_ = (IO::Uncompress::Gunzip_object, starting_pos)
  # Returns the number stored in the following 4 bytes of the
  # object, in little-endian format. Starts at <starting_pos> if passed
  # or at the current cursor point otherwise

  my $z   = shift;
  my $pos = shift;
  my $temp = '';
  if ($pos){
    $z->seek($pos, 0);
  }
  my $st = $z->{'z'}->read($temp, 4);
  my $result = unpack("V*", $temp);
  warn("Failed int32 read!\n") unless ($st == 4);
  return $result;
}

sub get_chars{

  my $z   = shift;
  my $n   = shift;
  my $pos = shift;
  my $temp = '';
  return $temp unless ($n);
  if ($pos){
    $z->{'z'}->seek($pos, 0);
  }
  my $st = $z->{'z'}->read($temp, $n);
  warn("Failed string read!\n") unless ($st == $n);
  $temp =~ s/\000//g;
  return $temp;
}

sub get_string{
  my $z   = shift;
  my $pos = shift;
  my $result = '';
  if ($pos){
    $z->{'z'}->seek($pos, 0);
  }
  my $char = 1;
  while ($char > 0){
    my $c = $z->getc();
    $result .= $c;
    $char = unpack("C", $c);
  }
  return $result;
}

################# Debug ###################
sub hash_to_string {
  my $hash_ref = shift;
  my $tab      = shift;
  $tab = 1 unless ($tab);
  my $result = '';
  my @keys   = ();
  if (ref($hash_ref) eq 'HASH' ||
      ref($hash_ref) eq 'HTTP::Response' ||
      ref($hash_ref) eq 'HTTP::Headers' ||
      ref($hash_ref) eq 'HTTP::Request'
      ){
    @keys = sort { lc($a) cmp lc($b) } keys %$hash_ref;
  }
  elsif (ref($hash_ref) eq 'ARRAY'){
    @keys = @$hash_ref;
  }
  foreach my $key ( @keys ) {
    my $key_val = '';
    if (ref($hash_ref) eq 'HASH' ||
        ref($hash_ref) eq 'HTTP::Response' ||
        ref($hash_ref) eq 'HTTP::Headers' ||
        ref($hash_ref) eq 'HTTP::Request'
    ){
      $key_val = $hash_ref->{$key};
    }
    if (ref($hash_ref) eq 'ARRAY'){
      $key_val = $key;
    }
    if ( ref( $key_val ) eq 'HASH' ) {
      my $text = $key_val;
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .=
          sp($tab)
        . "$key => $text\n"
        . hash_to_string( \%{ $key_val }, $tab + 1 );
    }
    elsif ( ref( $key_val ) eq 'HTTP::Headers' ) {
      my $text = HTTP::Headers->new;
      $text = $key_val;
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .=
          sp($tab)
        . "$key => $text\n"
        . hash_to_string( $text, $tab + 1 );
    }
    elsif ( ref( $key_val ) eq 'HTTP::Response' ) {
      my $text = HTTP::Response->new;
      $text = $key_val;
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .=
          sp($tab)
        . "$key => $text\n"
        . hash_to_string( $text, $tab + 1 );
    }
    elsif ( ref( $key_val ) eq 'HTTP::Request' ) {
      my $text = HTTP::Request->new;
      $text = $key_val;
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .=
          sp($tab)
        . "$key => $text\n"
        . hash_to_string( $text, $tab + 1 );
    }
    elsif ( ref( $key_val ) eq 'ARRAY' ) {
      $result .= sp($tab)."$key => $key_val\n";
      $tab = $tab + 1;
      for (my $i = 0; $i < @{$key_val}; $i++){
        my $text = $key_val->[$i];
        my $t = sp($tab) . sp( length("$i => "), " " );
        $text =~ s/\n/\n$t/g;
        $result .=
            sp($tab)
          . "$i => $text\n";
        if (ref($key_val->[$i]) eq 'HASH'){
          $result .= hash_to_string( \%{ $key_val->[$i] }, $tab + 1 );
        }
        elsif (ref($key_val->[$i]) eq 'ARRAY'){
          $result .= hash_to_string( \@{ $key_val->[$i] }, $tab + 1 );
        }
      }
      $tab = $tab - 1;
    }
    else {
      my $text = $key_val;
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .= sp($tab) . "$key => $text\n";
    }
  }
  return $result;
}

sub sp {
  my $n   = shift;
  my $pad = shift;
  $pad = "\t" unless ($pad);
  return '' unless ($n);
  my $result;
  for ( my $i = 0 ; $i < $n ; $i++ ) {
    $result .= $pad;
  }
  return $result;
}

__DATA__
form{
float: left;
}
#hdr{
float: left;
clear: both;
padding: 1em;
margin: 1em;
}
.hcont{
padding: 0;
margin: 0;
display: inline;
cursor: pointer;
}
#arrows{
float: left;
clear: both;
}
.fileName{
  font-family: Times New Roman;
  color: #000;
  background-color: #ddffdd;
  font-size: larger;
  width: 100%;
}

#sam_text{
	width: 90%;
	height: 8em;
	overflow: scroll;
	overflow-y: scroll;
	overflow-x: scroll;
}
span{
padding: 1px;
}
.highlighted{
padding: 0;
border-right: 1px solid #000000;
border-left: 1px solid #000000;
}
.genomic{
cursor: pointer;
background-color: #44ff00;
}
.result{
float: left;
clear: both;
padding: 1em;
position: relative;
font-family: Courier New;
font-size: x-small;
font-weight: bold;
}
.result_frame{
float: left;
clear: both;
padding: 1em;
position: relative;
border: thick groove #000000;
font-family: Courier New;
font-size: x-small;
font-weight: bold;
}
#light{
	position: relative;
	float: left;
	clear: both;
	width: 1em;
	height: 1em;
}
.greenLight{
	background-color: #00ff00;
}
.redLight{
	background-color: #ff0000;
}
.reads{
float: left;
clear: both;
padding: 1em;
position: relative;
border: thick groove #000000;
}
.hit{
color: #FFFFFF;
}
.sam{
float: left;
padding: 1em;
clear: both;
margin-top: 1em;
}
.insertion{
  padding: 0;
  margin-left: -3px;
  border-left: 3px solid #00ff00;
}

/* popups  */
  .pop{
    background-color: #CCCCCC;
    position: absolute;
    overflow: visible;
    z-index: 999;
    display: inline;
    margin: 0.2em;
    border: outset;
  }
  .popbar{
    font-family: Arial;
    font-size: x-small;
    /*padding: 0.5%;*/
    background-color: #AAAAAA;
    border-bottom: inset;
    position: static;
    height: auto;
  }
  .close{
    font-family: Arial;
    cursor: pointer;
    font-size: small;
    position: absolute;
    background-color: #555555;
    color: #ffffff;
    top: 0;
    right: 0;
    vertical-align: middle;
    border-style: outset;
    border: thin;
  }
  .close:hover{
  background-color: #994444;
  }

  .popup_table{
      border-style: double;
      width: auto;
  }
  .popup_table a{
    color: #0000aa;
  }

  .qgap{
  	color: #000000;
  }

  .popup_table td{
  font-size: x-small;
  border: thin solid #555555;
  padding: 0.5em;
  color: #000000;
  }
