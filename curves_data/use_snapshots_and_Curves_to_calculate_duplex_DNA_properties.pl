#!/usr/bin/perl
$ARGV[0] && $ARGV[1] && $ARGV[2] or die "Give the beginning time and ending time you want to calculate and the time interval you want to use\n";
my $start=$ARGV[0];
my $end=$ARGV[1];
my $dt = $ARGV[2];
my $num=0;
my $snaps=0;

####generate the snapshots from the MD trajecotry file"
system "echo 1 > file_select.txt";
$cycle = int(($end-$start)/(1000*$dt));
if ($cycle == 0){
   system "trjconv -f md1.xtc -s topol.md1.tpr -o frame..pdb -dt $dt -b $start -e $end -split 1 -nice -19 -pbc whole < file_select.txt";
   print "Renaming the frame pdb files.....\n";
      foreach (glob "frame*"){
              @name=split /\./, $_;
              $num=sprintf ("%09d", $name[1]);
              system "mv $_ snapshot.$num.pdb";
              }
}
elsif ($cycle > 0){
   for($i=0; $cycle>$i; $i++){
   my $m = $start + 1000*$i*$dt;
   my $n = $start + 1000*($i+1)*$dt;
   system "trjconv -f md1.xtc -s topol.md1.tpr -o frame..pdb -dt $dt -b $m -e $n -split 1 -nice -19 -pbc whole < file_select.txt";
print "Totally $cycle cycles, currently renaming the frame pdb files for cycle $i.....\n";
      foreach (glob "frame*"){
              @name=split /\./, $_;
              $num=sprintf ("%09d", ($name[1]+1000*$i));
              system "mv $_ snapshot.$num.pdb";
              }
   }
}
system "rm file_select.txt";

####calculate DNA properties using Curves+
@pdbfile=glob "snapshot*";
###if there is a total_ file, delete it
   system "rm -r total_*";
print "Calculate DNA properties using Curves+, wait......\n";
$snaps = @pdbfile;
my $count =1;
foreach (@pdbfile){
   
   print "Totally $snaps snapshots, currently is $count.\n";
   system "sed -i 's/snapshot.pdb/$_/g' duplex_DNA_or_RNA_parameter.sh"; 
   system "sed -i 's/outname/out.$_/g' duplex_DNA_or_RNA_parameter.sh";
   system "./duplex_DNA_or_RNA_parameter.sh";
   system "cat out.$_\_b.pdb >> total_backbone_splines_b.pdb"; 
#   system "cat out.$_.cda >> total_analysis_for_Canal.cda"; 
   system "grep TITLE $_ >> total_DNA_properties_for_every_snapshot.lis"; 
   system "cat out.$_.lis >> total_DNA_properties_for_every_snapshot.lis"; 
   system "cat out.$_\_X.pdb >> total_helical_axis_X.pdb"; 
   
   system "rm out.$_\_b.pdb out.$_.cda out.$_.lis out.$_\_X.pdb";
   system "rm $_";
   system "sed -i 's/$_/snapshot.pdb/g' duplex_DNA_or_RNA_parameter.sh";
   system "sed -i 's/out.$_/outname/g' duplex_DNA_or_RNA_parameter.sh";
   $count++;
   }
print "The calculated DNA properties are in the file: total_DNA_properties_for_every_snapshot.lis\n"; 
