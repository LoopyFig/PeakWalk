#Author: Sami Teeny
#Date: Sep 23, 2021
#Modified original msconvert to work on linux server

use Switch;
use Cwd qw(abs_path);

# Default options
#1) raw files location F:\predictHD
$raw_files_loc=".";

#2) CDF file location
$cdf_files_loc=".";

#3) Conversion mode
$conversion_mode="profile"; #options: "centroid" or "profile"

#4)Number of batches
$numbatches=1;

#5)directory structure: hierarchical (batch1, batch2) or regular (all files in one folder)
$dirstructure="regular";

#6) select files
$select_files_name=""; #"F:\\predictHD\\hilicposfilenames.txt";

############################No changes needed below this line###########################

$j = 0;
$optionLen = $#ARGV;
while ($j <= $optionLen) {
	$option = $ARGV[$j];
	print $option;
	switch ($option) {
		case "-i" {$j++; $raw_files_loc=$ARGV[$j];}
		case "-o" {$j++; $cdf_files_loc=$ARGV[$j];}
		case "-c" {$j++; $conversion_mode=$ARGV[$j];}
		case "-b" {$j++; $numbatches=$ARGV[$j];}
		case "-d" {$j++; $dirstructure=$ARGV[$j];}
		case "-s" {$j++; $select_files_name=$ARGV[$j];}
		case "-h" {print "Converts Raw Xcallibur Output to XLML\nOptions:\n-i input directory\n-o output directory\n-c conversion mode {profile, centroid}\n-b number of batches\n-d directory structure {regular, hierarchical}\n-s select files\n";exit();}
	}
	$j++;
}

#change directory to raw files location
$raw_files_loc = abs_path($raw_files_loc);
$cdf_files_loc = abs_path($cdf_files_loc);
system("cd ".$raw_files_loc);
if (not (-d $cdf_files_loc)) {
	system("mkdir $cdf_files_loc");
}

%filenames = ();
if(length($select_files_name)>1){
	open(inf1, "<".$select_files_name);
	while($line=<inf1>){
		chomp $line;
		#print $line,"\n";
		$filenames{$line}=1;
	}
}

#open new file in the raw files folder to write file names
open(outf, ">".$raw_files_loc."/Rawfiles.txt");

for($i=1;$i<=$numbatches;$i++){
	if($dirstructure eq "hierarchical"){
		$raw_files_locsub=$raw_files_loc."/batch".$i;
		$raw_files_docsub="/raw/batch".$i;
	}else{
		$raw_files_locsub=$raw_files_loc;
		$raw_files_docsub="/raw";
	}
	#get file names of the files in the folder
	opendir(DIR, $raw_files_locsub) or die $!;
	while (my $file = readdir(DIR)) {
		if ($file =~ m/\.raw/){
			if(length($select_files_name)>0){
				if(exists($filenames{$file})){
					print $file,"\n";
					#print outf "$raw_files_locsub/$file\n";
					print outf "$raw_files_docsub/$file\n";
				}
			}else{
				print $file,"\n";
				#print outf "$raw_files_locsub/$file\n";
				print outf "$raw_files_docsub/$file\n";
			}
		}
  }
  closedir(DIR);
}
close(outf);

if($conversion_mode eq "centroid"){
	#Execute MSConvert in centroid mode with PeakPicking option set to true
	#$cmd="/usr/local/bin/pwiz/msconvert -f Rawfiles.txt -o $cdf_files_loc --64 --zlib --mzXML --filter \"peakPicking true 1-\"";
	$cmd="docker run -t --rm -e WINEDEBUG=-all -v ".$raw_files_loc.":/raw -v ".$cdf_files_loc.":/cdf chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert -f /raw/Rawfiles.txt -o /cdf --64 --zlib --mzXML --filter \"peakPicking true 1-\""
}else{
	#Execute MSConvert
	#$cmd="/usr/local/bin/pwiz/msconvert -f Rawfiles.txt -o $cdf_files_loc --64 --zlib --mzXML";
	$cmd="docker run -it --rm -e WINEDEBUG=-all -v ".$raw_files_loc.":/raw -v ".$cdf_files_loc.":/cdf chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert -f /raw/Rawfiles.txt -o /cdf --64 --zlib --mzXML"
}
print $cmd;
system($cmd);
