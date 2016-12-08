<?php
# read the html file
$qc_file = fopen("qualimapReport.html",'r') or die("Unable to open file!");
$contents = fread($qc_file, filesize("qualimapReport.html"));
#close the open handle
fclose($qc_file);
# remove the \n or \r 
$contents = preg_replace('/\r|\n/','',$contents);
preg_match_all('/<h3>Globals<\/h3><table class[^>]*?>(.*?)<\/table>/si', $contents, $G1);
preg_match_all('/<h3>Globals.+?<\/h3><table class[^>]*?>(.*?)<\/table>/si', $contents, $G2);
//echo $G1[1][0]."\n\n";  //echo $G2[1][0]."\n\n";   //echo count($G2)."\n\n";   //echo $contents;

//define a file name
$name='qualimapReport.html';
// check if a file exists;
if(!file_exists($name)){ 
    echo $name." not found!"; 
    exit; 
} 
//file function read file content
$lines_array = file($name);

//transform array to string
$lines_string=implode('',$lines_array);
//echo $lines_string;

$lines_string = preg_replace('/\r|\n/','',$lines_string);
preg_match_all('/<h3>Globals<\/h3><table class[^>]*?>(.*?)\/table>/si', $lines_string, $test_gbl);
preg_match_all('/<h3>Globals.+?<\/h3><table class[^>]*?>(.*?)<\/table>/si', $lines_string, $region);

preg_match_all('/<td class=column1>Reference size<\/td>[^>]*?>(.*?)<\/td>/si', $test_gbl[1][0], $res_global_rs);
preg_match_all('/Number of reads<\/td>[^>]*?>(.*?)<\/td>/si', $test_gbl[1][0], $res_global_reads);
preg_match_all('/Mapped reads<\/td>[^>]*?>(.*?)<\/td>/si', $test_gbl[1][0], $res_global_mapped);
preg_match_all('/\/percentage of reference<\/td>[^>]*?>(.*?)<\/td>/si', $region[1][0], $res_region);
preg_match_all('/Mapped reads<\/td>[^>]*?>(.*?)<\/td>/si', $region[1][0], $res_region_mapped);

$region_a = split(" \/ ", $res_region[1][0]);

print_r($res_region);
//echo count($region_a)."\n";
echo "Reference Size\t".$res_global_rs[1][0]."\n";
echo "Number of Reads\t".$res_global_reads[1][0]."\n";
echo "Mapped Reads Global\t".$res_global_mapped[1][0]."\n";
echo "Regions Size\t".$region_a[0]."\n";
echo "Percentage of Reference\t".$region_a[1]."\n";
echo "Mapped Reads Region\t".$res_region_mapped[1][0]."\n";

?>
