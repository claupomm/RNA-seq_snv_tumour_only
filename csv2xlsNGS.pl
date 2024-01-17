#!/usr/bin/perl
use Spreadsheet::WriteExcel;
use Getopt::Long;

#GetOptions(
#    	   "file=s", \$file,
#	   );


opendir(in,'./') || die $!;
@files=readdir(in);
close in;

opendir(DIR,'./') || die $!;
while($file = readdir(DIR)) { if($file =~/csv/)
{
    my @pvalue;
    my @symbol;
    my @logfc;
    my @Entrez;

@filo = split(/\//,$file);
$newfile = @filo[$#filo];
$newfile = substr($newfile,0,length($newfile)-3); 

my $workbook = Spreadsheet::WriteExcel->new($newfile."xls"); 
my $worksheet = $workbook->add_worksheet();
my $num1_format  = $workbook->add_format(num_format => '0.00%');
my $num2_format  = $workbook->add_format(num_format => '0.00E+00');
my $num3_format  = $workbook->add_format(num_format =>'[BLUE]0.00;[RED]-0.00');
my $num4_format  = $workbook->add_format(num_format => '0.00');
my $bold       = $workbook->add_format(bold => 1);


open(FH,"<".$file) or die "Cannot open file: $!\n";
my ($x,$y) = (0,0);

while ($line=<FH>){ 
 chomp;
 $line =~s/\n//g;
 $line =~s/\n//g;
 
 @list = split(/\t/, $line);


 foreach my $c (@list){
     if($x==0)
     {
	 if($c=~/FDR/||$c=~/padj/||$c=~/pvalue/||$c=~/P.Value/||$c=~/adj.P.Val/){
	     @pvalue = ($y,@pvalue);
	 }
	 elsif($c=~/symbol/||$c=~/external_gene_id/||$c=~/external_gene_name/){
	     @symbol = ($y,@symbol);
	 }
	 elsif($c=~/log2FC/||$c=~/log2FoldChange/){
	     @logfc = ($y,@logfc);
	 }
	 elsif($c=~/Entrez/||$c=~/^gi/||$c=~/entrez/||$c=~/^chr/||$c=~/^start/||$c=~/^end/||$c=~/^length/){
	     @Entrez = ($y,@Entrez);
	 }
	 
	 $worksheet->write($x, $y++, $c,$bold);	 
     }
     else{
	 if((grep {$_==$y;} @pvalue)>0)
	 {
	     $worksheet->write($x, $y++, $c,$num2_format);
	 }
	 elsif((grep {$_==$y;} @symbol)>0)
	 {
	     $worksheet->write($x, $y++, $c,$bold);   
	 }
	 elsif((grep {$_==$y;} @logfc)>0)
	 {
	     $worksheet->write($x, $y++, $c,$num3_format);   
	 }
	  elsif((grep {$_==$y;} @Entrez)>0)
	 {
	     $worksheet->write($x, $y++, $c);   
	 }
	 else {
	     $worksheet->write($x, $y++, $c,$num4_format);
	 }
     }
 }
 $x++;$y=0;
}
close(FH);
$workbook->close();
}
}
closedir(DIR);
