#!/usr/bin/perl

use strict;
use warnings;
use threads;
use Thread::Queue;
our %hs_arr;
my $filename = $ARGV[0];
my $num_threads = 20;  # You can adjust this based on your machine's capabilities

# Queues for threads
my $work_queue = Thread::Queue->new();
my $result_queue = Thread::Queue->new();

# Thread pool creation
my @threads;
for (1..$num_threads) {
    push @threads, threads->create(\&process_lines);
}

# Main thread reading the file and enqueuing work
open my $fh, '<', $filename or die "Cannot open file $filename: $!";
while (my $line = <$fh>) {
    chomp $line;
    $work_queue->enqueue($line);
}
close $fh;
$work_queue->end();

for my $thr (@threads){
	$thr->join();
}

#print "all threads completed\n";
# Collecting results
#my @results;
while (defined(my $result = $result_queue->dequeue_nb())) {
    print "$result";
}

#exit;
# Wait for all threads to complete
#threads->exit();
# Output the results

sub process_lines {
    while (defined(my $line = $work_queue->dequeue())) {
        # Process the line (this is your existing line processing logic)
        my $result = process_line($line);
        $result_queue->enqueue($result);
    }
}

sub process_line {
    my ($line) = @_;
    #my $i++;
    #print "$i\t$line\n";
    my %hs_qrc;	
    # Your existing line processing logic goes here
    my (@methy, @refmethy, @ar);
    my ($qlen, $direct, $rs, $i, $methyref);
    @ar=split(/\t/, $line);
    $qlen = qlencigar($ar[5]);
    $direct = ($ar[1] & 16) ? 1 : 0;
    #print "$qlen\t$direct\n";
    ($rs, %hs_qrc) = cigarparse($ar[5], $qlen, $ar[3], $direct);
    %hs_arr = ();
    @methy=split(/,/,$ar[6]);
    for($i=0;$i<=$#methy;$i++){
	if($hs_qrc{$methy[$i]}){
		push @refmethy,$hs_qrc{$methy[$i]};
	}else{
		push @refmethy,-1;
	}
    }
    $rs = $rs -1;
    $methyref = join(',', @refmethy);
    #print "$rs\t$ar[6]\t$methyref\n";
    my $output="$line\t$methyref\t$qlen\t$rs\t$direct\n";
    #print "$line\n";	
    return $output;  # return the processed line
}

sub cigarparse{
	my ($cigar,$qlen,$refstart,$direct)=@_;
	my ($q_s, $q_e, $r_s, $r_e);
	my ($num, $type);
	my %hs_q2r;
	$r_s = $refstart;
	if($direct==0){
		$q_s=1;
		while($cigar=~/(\d+)(\D)/g){
			$num = $1;
			$type = $2;
			if($type eq '='){
				$q_e = $q_s+$num;
				$r_e = $r_s+$num;
				%hs_q2r = hasharr($q_s, $q_e, $r_s, $r_e);
				}
			if($type eq 'X'){
				$r_e=$r_s+$num;
				$q_e=$q_s+$num;
				}
			if($type eq 'D'){
				$r_e=$r_e+$num;
				}
			if($type eq 'I'){
				$q_e=$q_e+$num;
				}
			if($type eq 'S'){
				$q_s=$q_s+$num;
				next;
				}
			$q_s=$q_e;
			$r_s=$r_e;
		}
	}else{
                $q_s=$qlen;
                while($cigar=~/(\d+)(\D)/g){
                        $num = $1;
                        $type = $2;
                        if($type eq '='){
                                $q_e = $q_s-$num;
                                $r_e = $r_s+$num;
                                %hs_q2r = hasharr($q_s, $q_e, $r_s, $r_e);
                                }
                        if($type eq 'X'){
                                $r_e=$r_s+$num;
                                $q_e=$q_s-$num;
                                }
                        if($type eq 'D'){
                                $r_e=$r_e+$num;
                                }
                        if($type eq 'I'){
                                $q_e=$q_e-$num;
                                }
                        if($type eq 'S'){
                                $q_s=$q_s-$num;
                                next;
                                }
                        $q_s=$q_e;
                        $r_s=$r_e;
                }
	
	}
	return ($r_s, %hs_q2r);
}

sub qlencigar{
	my ($cigar) = @_;
	#print "$cigar\n";
	my $qlen = 0;
	my ($num,$type);
	while($cigar=~/(\d+)(\D)/g){
		$num = $1;
		$type = $2;
		#print "$num\t$type\n";
		if($type=~/[=XSI]/){
			$qlen+=$num;
		}
	}
	return $qlen;
}

sub hasharr{
	my ($qs, $qe, $rs, $re)=@_;
	my ($i,$j);
	my (@keys, @values);
	#our %hs_arr;
	if($qs > $qe){
		for($i=$qs;$i>$qe;$i--){
			push @keys,$i;
		}
	}else{
		for($i=$qs;$i<$qe;$i++){
			push @keys,$i;
		}
	}
	if($rs > $re){
		for($j=$rs;$j>$re;$j--){
			push @values,$j;
		}
	}else{
		for($j=$rs;$j<$re;$j++){
			push @values,$j;
		}
	}
	for(my $m = 0; $m < scalar(@keys); $m++){
        	$hs_arr{$keys[$m]} = $values[$m];
    	}
	return %hs_arr;
}

