<?php
require_once 'class.Numerical.php';
require_once 'MyStatistics.php';
class ReferenceGenes
{
   var $expression;
   var $referenceAvg;
   var $genesNumber=0;
   var $treatmentNumber=0;   
   //best keeper
	var $bkGeomena=array();
	var $bkarmena=array();
	var $bkMin=array();
	var $bkMax=array();
	var $bkStdDev=array();
	var $bkCV=array();
	var $bkMinxFold=array();
	var $bkMaxFold=array();
	var $bkStdFold=array();
	var $bkBP=array();
	var $bkIndividualGeneValue=array();
	
	var	$bkcompareValue=array();
	var $bkcomparePValue=array();
	
	var $bkcompareBKValue=array();
	var $bkcompareBKPValue=array();
	
	var $bkFinalIndex;
	
	
	
	//normfinder
	var $normfinderResult=array();
	
	
	var $genormResult=array();
 
   function __construct($data)
	{
	     $arr = explode("\r",$data);
		 $i=0;
         foreach ($arr as $row)
		 {
		   $row=trim(str_replace(array("\r","\n"), '', $row));
		   if($row !="")
		   {		      
		      $this->expression[$i++]=preg_split("/\s+/",$row);
		   }		   
		 }
		 $this->genesNumber=count($this->expression[0]);
		 $this->treatmentNumber=$i-1;
 
		 $this->calculate();
		 $this->BestKeeper();
		 $this->normFinder();
		 $this->genorm();
	}

//ct approach	


  function calculate()
  {  
	  $this->referenceAvg=array();
	  for($gene =0; $gene< $this->genesNumber;$gene++)
	  {
			$stdevAvg=array();
			for ($i = 0 ; $i< $this->genesNumber;$i++)
			{	    
				If ($i <> $gene )
				{            
				   $substractValue=array();
					for ($row=1;$row<$this->treatmentNumber+1; $row++)
					{
					   array_push($substractValue, $this->expression[$row][$gene]- $this->expression[$row][$i] );              
					}
					array_push($stdevAvg,  Numerical::standardDeviation($substractValue)); 
			   }
		   }   
		  $this->referenceAvg[$this->expression[0][$gene]]=Numerical::mean($stdevAvg);  
	  }	  
  }

 function BestKeeper()
 {
    for($gene =0; $gene< $this->genesNumber;$gene++)
	{
	  $members=array();
	  $this->bkIndividualGeneValue[$gene]=array();
	  for ($row=1;$row<$this->treatmentNumber+1; $row++)
	  {
	     #$sum=$sum+$this->expression[$row][$gene];
		 array_push($members,floatval($this->expression[$row][$gene]) );  		 
		 array_push($this->bkIndividualGeneValue[$gene],floatval($this->expression[$row][$gene]) ); 
 
	  } 
	 
	  $geomean=LeonStat::GEOMEAN($members);
	  $mean=LeonStat::AVERAGE($members);
	  $min=LeonStat::MIN($members);
	  $max=LeonStat::MAX($members);
	  $stdev=LeonStat::AVEDEV($members);
	  $cv = $stdev / $mean *100;
	  $minfold= -pow(2,abs($min- $geomean));
	  $maxfold= pow(2,$max- $geomean);
	  $stdFold=pow(2,$stdev);
    	  
		array_push($this->bkGeomena,  $geomean); 
		array_push($this->bkarmena, $mean);  
		array_push($this->bkMin,  $min); 
		array_push($this->bkMax,  $max); 
		array_push($this->bkStdDev, $stdev); 
		array_push($this->bkCV,  $cv); 
		
		$this->bkFinalIndex[$this->expression[0][$gene]]=$stdev;  
		
		array_push($this->bkMinxFold,  $minfold); 
		array_push($this->bkMaxFold,  $maxfold); 	
		array_push($this->bkStdFold,  $stdFold); 
	}
	
	
	  for ($row=1;$row<$this->treatmentNumber+1; $row++)
	  {
	    $members=array();
	    for ($gene =0; $gene< $this->genesNumber;$gene++)
		{
		  array_push($members,floatval($this->expression[$row][$gene]) ); 
		}		
		array_push($this->bkBP,LeonStat::GEOMEAN($members));
	  }
	
	

	for($gene1 =0; $gene1< $this->genesNumber;$gene1++)
	{
		for($gene2 =$gene1+1; $gene2< $this->genesNumber;$gene2++)
		{
		  $this->bkcompareValue[$gene1][$gene2]=LeonStat::CORREL($this->bkIndividualGeneValue[$gene1],$this->bkIndividualGeneValue[$gene2]);
		  $this->bkcomparePValue[$gene1][$gene2]=LeonStat::TDIST(abs($this->bkcompareValue[$gene1][$gene2]/sqrt((1-$this->bkcompareValue[$gene1][$gene2] *$this->bkcompareValue[$gene1][$gene2])/($this->treatmentNumber -2))), ($this->treatmentNumber -2), 2);		
		}
	 $this->bkcompareBKValue[$gene1]=LeonStat::CORREL($this->bkIndividualGeneValue[$gene1],$this->bkBP);
	 $this->bkcompareBKPValue[$gene1]=LeonStat::TDIST(abs($this->bkcompareBKValue[$gene1]/sqrt((1- $this->bkcompareBKValue[$gene1] * $this->bkcompareBKValue[$gene1])/($this->treatmentNumber -2))), ($this->treatmentNumber -2),2);
	}
 
 
 
 
}

function reportBestKeeper()
{
	echo "<a name=\"detact\" id=\"bestkeeperAnchor\">&nbsp;</a><table id='bestkeeper'><tr align='center'><td colspan='".($this->genesNumber+1)."'>CP data of housekeeping Genes by BEST KEEPER</td></tr><tr><td>&nbsp;&nbsp;</td>"; 
	
	foreach ($this->expression[0] as $k)
	{ 
	  echo "<td>".$k."</td>";
	}   
	echo "</tr>";

	//n
	echo "<tr><td>n</td>"; 
	for($gene =0; $gene< $this->genesNumber;$gene++)
	{
	echo "<td>".$this->treatmentNumber."</td>";
	}   
	echo "</tr>";
 
	
	//geo Mean [CP]
	echo "<tr><td>geo Mean [CP]</td>"; 
	for($gene =0; $gene< $this->genesNumber;$gene++)
	{
	echo "<td>".number_format($this->bkGeomena[$gene],2, '.', '')."</td>";
	}   
	echo "</tr>";

	//ar Mean [CP]
	echo "<tr><td>AR Mean [CP]</td>"; 
	for($gene =0; $gene< $this->genesNumber;$gene++)
	{
	echo "<td>".number_format($this->bkarmena[$gene],2, '.', '')."</td>";
	}   
	echo "</tr>";	
	
	
    //min [CP]
	echo "<tr><td>min [CP]</td>"; 
	for($gene =0; $gene< $this->genesNumber;$gene++)
	{
	echo "<td>".number_format($this->bkMin[$gene],2, '.', '')."</td>";
	}   
	echo "</tr>";	
	
    //min [CP]
	echo "<tr><td>max [CP]</td>"; 
	for($gene =0; $gene< $this->genesNumber;$gene++)
	{
	echo "<td>".number_format($this->bkMax[$gene],2, '.', '')."</td>";
	}   
	echo "</tr>";	

	    //std dev [¡À CP]
	echo "<tr style='color: #0000FF'><td>std dev [+/- CP]</td>"; 
	for($gene =0; $gene< $this->genesNumber;$gene++)
	{
	echo "<td>".number_format($this->bkStdDev[$gene],2, '.', '')."</td>";
	}   
	echo "</tr>";	
	
	
	    //CV [% CP]
	echo "<tr><td>CV [% CP]</td>"; 
	for($gene =0; $gene< $this->genesNumber;$gene++)
	{   
	echo "<td>".number_format($this->bkCV[$gene], 2, '.', '')."</td>";
	}   
	echo "</tr>"; 

		    //min [x-fold]
	echo "<tr><td>min [x-fold]</td>"; 
	for($gene =0; $gene< $this->genesNumber;$gene++)
	{
	echo "<td>".number_format($this->bkMinxFold[$gene],2, '.', '')."</td>";
	}   
	echo "</tr>"; 
	
		    //max [x-fold]
	echo "<tr><td>max [x-fold]</td>"; 
	for($gene =0; $gene< $this->genesNumber;$gene++)
	{
	echo "<td>".number_format($this->bkMaxFold[$gene],2, '.', '')."</td>";
	}   
	echo "</tr>";
 
    //max [x-fold]
	echo "<tr><td>std dev [+/- x-fold]</td>"; 
	for($gene =0; $gene< $this->genesNumber;$gene++)
	{
	echo "<td>".number_format($this->bkStdFold[$gene],2, '.', '')."</td>";
	}   
	echo "</tr>";	
    echo "</table>";	
	
	
	//BEST KEEPER
echo "<br/><br/><table id='bestkeeper'><tr align='center'><td colspan='".($this->genesNumber+1)."'> Pearson correlation coefficient ( r ) by BEST KEEPER</td></tr><tr><td>&nbsp;&nbsp;</td>"; 
	foreach ($this->expression[0] as $k)
	{ 
	  echo "<td>".$k."</td>";
	}   
	echo "</tr>";	
	

   for($gene1 =1; $gene1< $this->genesNumber;$gene1++)	
   {
      if (($gene1%2)==0)
	  {
	    $style="style='color: #0000FF'";
	  }
	  else
	  {
	     $style="";
	  }

	    echo "<tr ". $style."><td>".$this->expression[0][$gene1]."</td>"; 
       for($gene2 =0; $gene2< $this->genesNumber;$gene2++)	
       {
	     if($this->bkcompareValue[$gene2][$gene1] =='')
		 {
		 echo "<td> - </td>"; 
		 }
		 else
		 {
          echo "<td>".number_format($this->bkcompareValue[$gene2][$gene1] ,3, '.', '')."</td>"; 
		  }
       }	   
	   echo "</tr>";
		
	    echo "<tr ". $style."><td>p-value</td>"; 
       for($gene2 =0; $gene2< $this->genesNumber;$gene2++)	
       {
	   	 if($this->bkcomparePValue[$gene2][$gene1]  =='')
		 {
		   echo "<td> - </td>"; 
		 }
		 else
		 { 
		   if($this->bkcomparePValue[$gene2][$gene1]  <0.001)
		   {
             echo "<td>0.001</td>"; 
		   }
		   else
		   {
		     echo "<td>".number_format($this->bkcomparePValue[$gene2][$gene1] ,3, '.', '')."</td>"; 
		   }
		  }
       }	   
	   echo "</tr>";	
		
   }
echo "</table>";

   echo "<br/><br/><table id='bestkeeper'><tr align='center'><td colspan='".($this->genesNumber+1)."'>Pearson correlation coefficient ( r )</td></tr><tr><td>BestKeeper  vs.</td>"; 
	foreach ($this->expression[0] as $k)
	{ 
	  echo "<td>".$k."</td>";
	}   
	echo "</tr>";	
	
      echo "<tr><td>coeff. of corr. [r]</td>"; 
       for($gene =0; $gene< $this->genesNumber;$gene++)	
       {
	   	 if($this->bkcompareBKValue[$gene]  =='')
		 {
		   echo "<td> &nbsp;</td>"; 
		 }
		 else
		 { 
		   if($this->bkcompareBKValue[$gene]  <0.001)
		   {
             echo "<td>0.001</td>"; 
		   }
		   else
		   {
		     echo "<td>".number_format($this->bkcompareBKValue[$gene]  ,3, '.', '')."</td>"; 
		   }
		  }
       }	   
	   echo "</tr>";	
	   
      echo "<tr><td>p-value</td>"; 
       for($gene =0; $gene< $this->genesNumber;$gene++)	
       {
	   	 if($this->bkcompareBKPValue[$gene]  =='')
		 {
		   echo "<td> &nbsp;</td>"; 
		 }
		 else
		 { 
		   if($this->bkcompareBKPValue[$gene]  <0.001)
		   {
             echo "<td>0.001</td>"; 
		   }
		   else
		   {
		     echo "<td>".number_format($this->bkcompareBKPValue[$gene]  ,3, '.', '')."</td>"; 
		   }
		  }
       }	   
	   echo "</tr>";     
    echo "</table>";	
}



function normFinder()
{
   // $InputRange,$MaxColumn,$MaxRow,$GeneCount,$SampleCount,$z,$varTest,$x,$y,$SheetName,
   // $NewSheetName,$SheetNumber,$NewSheet,$ws,$BestGene,$BestComb1,$BestComb2,
   // $k,$p,$GroupID_1,$GroupID_2,$i;
   
$groupNumber=0;
$DataSet;       #includes sample names, gene names and group identifiers
$Data=array();        #includes only data for (calculation
$ctrlgroup;     #Group identifier included
$ctrllog;       #Log transform data
$ctrlprint;     #Simple output only

$GeneNames;     #Gene names included
$SampleNames;   #Sample names included
 
$qtotal=array();# qualtity measure for (gene - global $variable
$qtotalord=array();# order for (quality of average - global $variable
$qtotalavg=array();# qualtity measure for (average - global $variable
$sig=array();# variance estimates for (gene and group - global $variable
$varsig=array();# variance on sig - global $variable
$gr=array(); # group label - global $variable
$grn; # number of groups - global $variable
$gamma; # estimate of variance component - global $variable
$dif=array();# residuals for (gavg - global $variable
$grnames=array(); #Groupidentifiernames


$i;
$j;
$m;
$sam;
$a1;
$a2;
$n; # number of samples
$k; # number of genes

$dat=array();# log transformed expression levels
$res=array();# residuals
$savg=array();# average for (each sample
$s1;
$s2; 
$s3;
$mi;
$di;
#$gr; # group label -  $variable
#$grn; # number of groups -  $variable
$grval=array(); # identifier for (each group
$gavg=array();# averages for (each gene and group
$gsavg=array();# averages for (each group
$ss=array();# sum of squared terms for (gene and group
#$sig;# variance estimates for (gene and group -  $variable
$ngr=array(); # number of samples in each group
$gavgrow=array();# row averages for (gavg
$gavgcol=array();# column averages for (gavg
$gavgtotal=array(); # total average for (gavg
#$gamma; # estimate of variance component -  $variable
#$dif;# residuals for (gavg -  $variable
$postmean=array();# posterior means
$postvar=array();# posterior variance
$qm=array();# quality measure for (gene and group
#$qtotal;# qualtity measure for (gene -  $variable
#$qtotalavg;# qualtity measure for (average -  $variable
#$qtotalord;# order for (quality of average -  $variable
$qmavg=array();# intermediate
$vardif=array();# variance on dif
#$varsig;# variance on sig -  $variable
$sigcol=array();# column sums for (sig/ngr
$sigrow=array();# row sums for (sig/ngr
$sigtotal=array(); # total sum for (sig/ngr
$tau=array();# for (calculation of varsig
$ome=array();# for (calculation of varsig
$pair=array();# used is an estimate become negative
$ome=array();

$SampleNames = 0;
$GeneNames = 1;
$ctrlgroup = 0;
$ctrllog = 1;
$ctrlprint=1;
$GeneCount = $this->genesNumber;
$SampleCount =  $this->treatmentNumber;
$genesNumber=$this->genesNumber;
$sampleNumber=  $this->treatmentNumber;

for($ii=1; $ii<=$SampleCount; $ii++)
{
  for($jj=0; $jj<$GeneCount;$jj++)
  {
     $Data[$jj][$ii-1]=$this->expression[$ii][$jj];
  }
}

if($ctrlgroup == 1) 
{ 
 
# there is more than one group

$n = $sampleNumber;
$k = $genesNumber;
# transforming to log values

if($ctrllog == 1)
{
for ($i = 0;$i<$k; $i++)
{
 for ($j =0; $j<$n;$j++)
 {
   $dat[$i][$j]  =log($Data[$i][$j]);
 }
}
}
else 
{
	for ($i = 0;$i<$k; $j++)
	{
		 for ($j =0; $j< $n; $j++)
		 {
		  $dat[$i][$j]  = $Data[$i][$j] ; 
		 }
	}
}

# finding the number of groups

$grn = 1;
$gr[0] = $Data[$k + 1-1][1-1];
$grnames[1-1] = $gr[1-1];
for ($j = 2-1;$j<$n; $j++)
{
 $gr[$j] = intval($Data[$k + 1-1][$j]);
 $a1 = 1;
 for ($i = 1-1; $i<$j - 1;$i++)
 {
  if(($gr[$i] == $gr[$j])) 
  {
   $a1 = 0;
  }
  }
 if(($a1 == 1)) 
 {
    $grn = $grn + 1;
    $grnames[$grn-1] = $gr[$j];
 }
}
 

# finding the group sizes

$grval[1-1] = $gr[1-1];
for ($m = 1-1; $m<$grn;$m++)
{
 $ngr[$m] = 0;
}

$a2 = 1;
for ($j = 2-1 ;$j< $n; $j++)
{
 $a1 = 1;
 for ($i = 1-1; $i<$j - 1; $i++)
 {
  if(($gr[$i]= $gr[$j])) 
  { 
   $a1 = 0;
  }
 if(($a1 == 1))
 { 
	 $a2 = $a2 + 1;
	 $grval[$a2] = $gr[$j]; 
 }
}
}

for ($j = 1-1; $j<$n; $j++)
{
 for ($m = 1-1;$i<$grn; $i++)
 {
  if(($gr[$j] == $grval[$m])) 
  { 
   $ngr[$m] = $ngr[$m] + 1;
  }
 }
}

# finding averages

for ($i =1-1;$i<$k;$i++)
{
 for ($m = 1-1;$m<$grn;$m++)
 {
  $s2 = 0;
  for ($j = 1-1;$j<$n;$j++)
  {
   if(($gr[$j] == $grval[$m])) 
   { 
    $s2 = $s2 + $dat[$i][$j];
   }
  $gavg[$i][$m] = $s2 / $ngr[$m];
 }
}
}

for ($m = 1-1;$m<$grn;$m++)
{
 $s1 = 0;
 for ($i = 1-1 ;$i<$k; $i++)
 {
  $s1 = $s1 + $gavg[$i][$m];
 }
 $gsavg[$m] = $s1 / $k;
}

for ($m = 1-1 ;$m<$grn; $m++)
{
 $gavgcol[$m] = $gsavg[$m];
}

for ($i = 1-1; $i<$k; $i++)
{
 $s1 = 0;
 for ($m = 1-1 ;$m<$grn; $m++)
 {
  $s1 = $s1 + $gavg[$i][$m];
 }
 $gavgrow[$i] = $s1 / $grn;
}

$s1 = 0;
for ($m = 1-1;$m<$grn;$m++)
{
 $s1 = $s1 + $gavgcol[$m];
}
$gavgtotal = $s1 / $grn;


for ($j = 1-1 ;$j<$n; $j++)
{
 $s1 = 0;
 for ($i = 1-1;$i<$k;$i++)
 {
  $s1 = $s1 + $dat[$i][$j];
 }
 $savg[$j] = $s1 / $k;
}

# calculation of residuals

for ($i = 1-1; $i<$k; $i++)
{
 for ($j = 1-1;$j<$n;$j++)
 {
  for ($m = 1-1; $m<$grn; $m++)
  {
   if(($gr[$j] == $grval[$m])) 
   { 
   $a1 = $m;
  }
  $res[$i][$j] = $dat[$i][$j] - $gavg[$i][$a1]- $savg[$j] + $gsavg[$a1];
 }
}
}

for ($i = 1-1; $i<$k;$i++)
{
 for ($m = 1-1; $m<$grn; $m++)
 {
  $dif[$i][$m] = $gavg[$i][$m] - $gavgcol[$m] - $gavgrow[$i] + $gavgtotal;
 }
}

# calculation of variances estimates

for ($i = 1-1; $i<$k; $i++)
{
 for ($m = 1-1; $m<$grn; $m++)
 {
  $s1 = 0;
  for ($j = 1-1;$j<$n;$j++)
  {
   if(($gr[$j] == $grval[$m])) 
   { 
     $s1 = $s1 + $res[$i][$j] *$res[$i][$j];
	}
  }
  $ss[$i][$m] = $s1 / (($ngr[$m] - 1) * (1 - 2 / $k));
 }
}

for ($m = 1-1; $m<$grn; $m++)
{
 $s1 = 0;
 for ($i =1-1;$i<$k; $i++)
 {
  $s1 = $s1 + $ss[$i][$m];
 }
 
 for ($i = 1-1;$i<$k; $i++)
 {
  $sig[$i][$m] = $ss[$i][$m]- $s1 / ($k * ($k - 1));
  #if((sig(i, m) < 0)) {
  #rettet d. 5-1-05 cgu/jlj
  if(($sig[$i][$m] <= 0)) 
  { #rettet d. 5-1-05 cgu/jlj
   for ($j = 1-1; $j<$k; $j++)
   {
    $s2 = 0;
    $s3 = 0;
    for ($sam = 1-1; $sam<$n; $sam++)
	{
     #if((gr(sam) = m)) { di = dat(i, sam) - dat(j, sam): s2 = s2 + di: s3 = s3 + di * di
     #rettet d. 5-1-05 cgu/jlj
     if(($gr[$sam] == $grval[$m])) 
	 { 
	    $di = $dat[$i][$sam] - $dat[$j][$sam];
		$s2 = $s2 + $di;
		$s3 = $s3 + $di * $di;
     }
	}
    $pair[$j] = ($s3 - $s2 * $s2 / $ngr[$m]) / ($ngr[$m] - 1);   
   }
   $s2 = 0;
   for ($j = 1-1; $j< $k; $j++)
   {
		if(($pair[$j] > $s2)) 
		{ 
		  $s2 = $pair[$j];
		}
   }
	
   for ($j = 1-1 ; $j< $k; $j++)
   {
    if(($pair[$j] > 0)) 
	{ 
	  if(($pair[$j] < $s2)) 
	  { 
	   $s2 = $pair[$j];
      }
	}
  }
   $sig[$i][$m] = $s2 / 4;
  }
 }
 
}

# calculation of variance component
$s1 = 0;
for ($i = 1-1; $i<$k; $i++)
{
 for ($m = 1-1;$m<$grn; $m++)
 {
  $s1 = $s1 + $dif[$i][$m] * $dif[$i][$m] / (($k - 1) * ($grn - 1)) - $sig[$i][$m] / ($ngr[$m] * $k * $grn);
 }
}

$gamma = 0;
if(($s1 > 0)) 
{
  $gamma = $s1;
}
# calculation of quality measure

for ($i = 1-1;$i<$k;$i++)
{
 for ($m = 1-1; $m<$grn;$m++)
 {
  $postmean[$i][$m] = $gamma * $dif[$i][$m] / ($gamma + $sig[$i][$m] / $ngr[$m]);
  $postvar[$i][$m] = $sig[$i][$m] / $ngr[$m] + $gamma * $sig[$i][$m] / ($ngr[$m] * ($gamma + $sig[$i][$m] / $ngr[$m]));
  $qm[$i][$m] = abs($postmean[$i][$m]) + sqrt($postvar[$i][$m]);
 }
}

for ($i = 1-1;$i<$k;$i++)
{
 $s1 = 0;
 for ($m = 1-1 ; $m<$grn;$m++)
 {
  $s1 = $s1 + $qm[$i][$m];
 }
 $qtotal[$i] = $s1 / $grn;
}

# Finding best genes and pair of best genes

$mi = $qtotal(1-1);
$a1 = 1;
for ($i = 1-1;$i<$k;$i++)
{
 if($qtotal[$i] < $mi) 
 { 
   $a1 = $i;
   $mi = $qtotal[$i];
 }
}
$qtotalavg[1-1] = $mi   ;                    # smallest index #
$qtotalord[1-1] = $a1  ;                 # gene with smallest index #


$mi = 100 * $qtotal[1-1] + 100;
for ($i = 1-1 ;$i<($k - 1);$i++)
{
 for ($j = ($i + 1);$j<$k;$j++)
 {
  for ($m = 1-1;$m<$grn;$m++)
  {
   $s1 = ($postmean[$i][$m] + $postmean[$j][$m]) / 2;
   $s2 = ($postvar[$i][$m] + $postvar[$j][$m]) / 4;
   $qmavg[$m] = abs($s1) * sqrt($k / ($k - 2)) + sqrt($s2);
  }
  $s1 = 0;
  for ($m = 1-1;$m<$grn;$m++)
  {
    $s1 = $s1 + $qmavg[$m];
  }
  
  $s1 = $s1 / $grn;
  if($s1 < $mi) 
  { 
    $a1 = $i;
	$a2 = $j;
	$mi = $s1;
  }
 }
}
$qtotalavg[2-1] = $mi;                      # smallest index based on two genes #
$qtotalord[2-1] = $a1  ;                    # gene 1 of two genes #
$qtotalord[3-1] = $a2   ;                   # gene 2 of two genes #

# calculation of variance on dif

$s2 = 0;
for ($i = 1-1;$i<$k; $i++)
{
   $s1 = 0;
	for ($m = 1-1; $m<$grn;$m++)
	{
	 $s1 = $s1 + $sig[$i][$m] / $ngr[$m];
	}
 $sigrow[$i] = $s1;
 $s2 = $s2 + $s1;
}
$sigtotal = $s2;

for ($m = 1-1; $m<$grn; $m++)
{
 $s1 = 0;
  for ($i = 1-1;$i<$k;$i++)
  {
    $s1 = $s1 + $sig[$i][$m];
  }
 $sigcol[$m]= $s1 / $ngr[$m];
}

for ($i = 1-1; $i<$k;$i++)
{
 for ($m =1-1; $m<$grn; $m++)
 {
  $vardif[$i][$m] = $k * $grn * ($k - 2) * ($grn - 2) * $sig($i, $m) / $ngr[$m] + $k * ($k - 2) * $sigrow[$i] + $grn * ($grn - 2) * $sigcol[$m] + $sigtotal;
  $vardif[$i][$m] = $vardif[$i][$m] / ($k * $k * $grn * $grn);
 }
}



# calculation of variance on the variance estimates

for ($m = 1-1;$m<$grn;$m++)
{
 $s1 = 0;
 for ($i = 1-1; $i<$k;$i++)
 {
  #tau(i) = (ngr(m) - 1) * ((1 - 2 / k) * sig(i, m) + sigcol(m) / (k * k)) / ngr(m)
  #rettet 5-1-05 cgu/jlj
  $tau[$i] = ($ngr[$m] - 1) * ((1 - 2 / $k) * $sig[$i][$m] + $sigcol[$m] * $ngr[$m] / ($k * $k)) / $ngr[$m];
  $tau[$i] = $tau[$i] * $tau[i];
  $s1 = $s1 + $tau[$i];
  for ($j = 1-1;$j<$k;$j++)
  {
   #ome(i, j) = -(ngr(m) - 1) * ((sig(i, m) + sig(j, m)) - sigcol(m) / k) / k / ngr(m)
   #rettet 5-1-05 cgu/jlj
   $ome[$i][$j] = -($ngr[$m] - 1) * (($sig[$i][$m] + $sig[$j][$m]) - $sigcol[$m] * $ngr[$m] / $k) / $k / $ngr[$m];
   $ome[$i][$j] = $ome[$i][$j] * $ome[$i][$j];
  }
 }
 
 
 for ($i = 1-1; $i<$k;$i++)
 {
  $s2 = 0;
  for ($j = 1-1 ;$j<$k; $j++)
  {
   $s2 = $s2 + $ome[$i][$j];
  }
  $s2 = $s2 - $ome[$i][$i];
  $s3 = 0;
  for ($j = 1-1; $j<$k;$j++)
  {
   for ($a1 = 1-1; $a1<$k;$a1++)
   {
    $s3 = $s3 + $ome[$j][$a1];
   }
   $s3 = $s3 - $ome[$j][$j] - $ome[$j][$i];
  }
  $varsig[$i][$m] = (1 - 2 / ($k * ($k - 1))) * $tau[$i] - (2 - 1 / ($k * ($k - 1))) * $s2 / ($k * ($k - 1)) + ($s1 + $s3) / ($k * $k) / (($k - 1) * ($k - 1));
  $varsig[$i][$m] = $varsig[$i][$m]* 2 * $ngr[$m] * $ngr[$m] / (($ngr[$m] - 1) * ($ngr[$m] - 1) * ($ngr[$m] - 1) * (1 - 2 / $k) * (1 - 2 / $k));
  $varsig[$i][$m] = sqrt($varsig[$i][$m] / $sig[$i][$m])/ 2;
 }
}

} 
else 
{ # there is one group only
$n = $sampleNumber;
$k = $genesNumber;
 
# transforming to log values
$ctrllog=0;
if($ctrllog == 1) 
{
	for ($i = 1-1; $i<$k;$i++)
	{
	 for ($j = 1-1; $j<$n;$j++)
	 {
	  $dat[$i][$j] = log($Data[$i][$j]);	  
	 }
 
	}
}  
else 
{
	for ($i =1-1;$i<$k; $i++)
	{
	 for ($j = 1-1;$j<$n;$j++)
	 {
	  $dat[$i][$j] = $Data[$i][$j];
	 }
	}
}

 

# finding averages

for ($i = 1-1 ;$i<$k;$i++)
{
 $s2 = 0;
 for ($j = 1-1 ;$j<$n;$j++)
 {
  $s2 = $s2 + $dat[$i][$j];
 }
 $gavg[$i] = $s2 / $n;
}



for ($j = 1-1; $j<$n; $j++)
{
 $s1 = 0;
 for ($i = 1-1 ;$i<$k; $i++)
 {
  $s1 = $s1 + $dat[$i][$j];
 }
 $savg[$j] = $s1 / $k;
}




$s1 = 0;
for ($j =1-1 ;$j<$n; $j++)
{
 $s1 = $s1 + $savg[$j];
}
$gavgtotal = $s1 / $n;


# calculation of residuals

for($i = 1-1 ;$i<$k; $i++)
{
 for  ($j = 1-1; $j<$n; $j++)
 {
  $res[$i][$j] = $dat[$i][$j] - $gavg[$i] - $savg[$j] + $gavgtotal;
 }
}




# calculation of variances estimates

for($i = 1-1 ;$i<$k; $i++)
{
  $s1 = 0;
  for ($j = 1-1; $j<$n; $j++)
  {
   $s1 = $s1 + $res[$i][$j] * $res[$i][$j];
  }
  $ss[$i] = $s1 / (($n - 1) * (1 - 2 / $k));
}
 
 
 
$s1 = 0;
for($i = 1-1 ;$i<$k; $i++)
{
  $s1 = $s1 + $ss[$i];
}


	

for($i = 1-1 ;$i<$k; $i++)
{
  $sig[$i] = $ss[$i] - $s1 / ($k * ($k - 1));
  #if((sig(i) < 0)) {
  #rettet d. 5-1-05 cgu/jlj
  if(($sig[$i] <= 0)) 
  {
   for ($j = 1-1; $j<$k;$j++)
   {
		$s2 = 0;
		$s3 = 0;
		for ($sam = 1-1; $sam<$n; $sam++)
		{
		 $di = $dat[$i][$sam] - $dat[$j][$sam];
		 $s2 = $s2 + $di;
		 $s3 = $s3 + $di * $di;
		}
		$pair[$j] = ($s3 - $s2 * $s2 / $n) / ($n - 1);
	
   }

	
   $s2 = 0;
   for ($j = 1-1; $j<$k;$j++)
   {
    if(($pair[$j] > $s2))
	{ 
	  $s2 = $pair[$j];
    }
   }
   
    for ($j = 1-1; $j<$k;$j++)
   {
    if(($pair[$j] > 0)) 
	 { 
	  if(($pair[$j] < $s2)) 
 	  { 
	    $s2 = $pair[$j];	  
      }
	}
	}
	$sig[$i] = $s2 / 4;
  }
   
 
}
$mi = $sig[1-1];
$a1 = 1-1;




for ($i = 2-1;$i<$k;$i++)
{
 if($sig[$i] < $mi)
 {
  $a1 = $i;
  $mi = $sig[$i];
 }
}
$qtotalord[1-1] = $a1;

# calculation of variance on the variance estimates

$s2 = 0;
for ($i = 1-1 ;$i<$k; $i++)
{
 $s2 = $s2 + $sig[$i];
}

 $s1 = 0;
 for ($i = 1-1; $i<$k; $i++)
 {
  $tau[$i] = ($n - 1) * ((1 - 2 / $k) * $sig[$i] + $s2 / ($k * $k)) / $n;
  $tau[$i] = $tau[$i] * $tau[$i];
  $s1 = $s1 + $tau[$i];
  for ($j = 1-1; $j<$k; $j++)
  {
   $ome[$i][$j] = -($n - 1) * (($sig[$i] + $sig[$j]) - $s2 / $k) / $k / $n;
   $ome[$i][$j] = $ome[$i][$j] * $ome[$i][$j];
  }
 }
 
 
 for($i = 1-1; $i<$k; $i++)
 {
  $s2 = 0;
  for ($j = 1-1; $j<$k; $j++)
  {
   $s2 = $s2 + $ome[$i][$j];
  }
  
  $s2 = $s2 - $ome[$i][$i];
  $s3 = 0;
  for ($j = 1-1;$j<$k; $j++)
  {
   for ($a1 = 1-1; $a1<$k; $a1++)
   {
    $s3 = $s3 + $ome[$j][$a1];
   }
   $s3 = $s3 - $ome[$j][$j] - $ome[$j][$i];
  }
  $varsig[$i] = (1 - 2 / ($k * ($k - 1))) * $tau[$i] - (2 - 1 / ($k * ($k - 1))) * $s2 / ($k * ($k - 1)) + ($s1 + $s3) / ($k * $k) / (($k - 1) * ($k - 1));
  $varsig[$i] = $varsig[$i] * 2 * $n * $n / (($n - 1) * ($n - 1) * ($n - 1) * (1 - 2 / $k) * (1 - 2 / $k));
  $varsig[$i] = sqrt($varsig[$i] / $sig[$i]) / 2;
 }
} #one or more groups

for($i=0; $i<$k; $i++)
{
 $this->normfinderResult[$this->expression[0][$i]]=sqrt($sig[$i]);
}
}


//genom algorithm
function genorm()
{
   $genormNormalizedDataByGene=array();
   $genes=array();
   $genormNormalizedDataByRow=array();
   $genormNormalizedData=array();
   
   $lastGeneStabilityValue=array();
   
   $minCTofGenes=array();
 
   
  for($j=0; $j<$this->genesNumber; $j++)
  {
    $temarray=array();
	$sum=0;
    for($i=1; $i<=$this->treatmentNumber; $i++)
	{
	  $temarray[$sum++]=floatval($this->expression[$i][$j]);	
	} 	
	$minCTofGenes[$j]=LeonStat::MIN($temarray);
  }    

   
 //stardalizing
  for($i=1; $i<=$this->treatmentNumber; $i++)
  {
    for($j=0; $j<$this->genesNumber; $j++)
	{
	   $genormNormalizedData[$i-1][$j]=pow(2, -(floatval($this->expression[$i][$j])-$minCTofGenes[$j]));
	} 	 
  }   
  
   
  //by row 
  for($i=0; $i<$this->treatmentNumber; $i++)
  {
    $temarray=array();
	$sum=0;
    for($j=0; $j< $this->genesNumber; $j++)
	{
	  $temarray[$sum++]=$genormNormalizedData[$i][$j];	
	} 
	$genormNormalizedDataByRow[$i-1]=$temarray;
  }

  //by gene or by column
  for($j=0; $j<$this->genesNumber; $j++)
  {
    $temarray=array();
	$sum=0;
    for($i=0; $i< $this->treatmentNumber; $i++)
	{
	  $temarray[$sum++]=$genormNormalizedData[$i][$j];	 
	} 
	$genormNormalizedDataByGene[$j]=$temarray;	  
	$genes[$j]=$this->expression[0][$j];

  } 
 

  //  spearMan ==> matrix_cv ==>avarge cv (averageM)
  $matrixs=array();  
  $newGenesNumber=$this->genesNumber;
  
  
 for($deleteI=1; $deleteI<= $this->genesNumber-1;$deleteI++ ) 
 {
   for($i=0; $i<$newGenesNumber; $i++)
   {
       for($j=$i+1; $j<$newGenesNumber; $j++)
       {
	     $logration=array();
	     if($i!=$j)
		 {
		    $twoLog=array();
		     for($ii=0; $ii< $this->treatmentNumber; $ii++)
			 {
			    array_push($twoLog, log($genormNormalizedDataByGene[$i][$ii] /$genormNormalizedDataByGene[$j][$ii]) /log(2) );
			 }		 
			$matrixs[$i][$j]= LeonStat::STDEV($twoLog);
			$matrixs[$j][$i]=$matrixs[$i][$j];			
		 }	
	   }  
		
   }
 
   $eachGeneAverageM=array();
   for($i=0; $i<$newGenesNumber; $i++)
   {
      $eachMatrixValue=array();
     for($j=0; $j<$newGenesNumber; $j++)
	 {
	   if(is_numeric( $matrixs[$j][$i]))
	   {
	   array_push($eachMatrixValue, $matrixs[$j][$i]);	 
	   }
	 } 
	 
	 array_push($eachGeneAverageM,LeonStat::AVERAGE($eachMatrixValue));
	 
   }
   

 
	  if($deleteI== $this->genesNumber-1)
	  {
		  $maximalOrder=$this->getMatchOrderNo($eachGeneAverageM);
		   $lastGeneStabilityValue[$genes[0]." | ".$genes[1]]=LeonStat::AVERAGE($eachGeneAverageM);	  
  
	  }else
	  {
		   $maximalOrder=$this->getMatchOrderNo($eachGeneAverageM);
		   $lastGeneStabilityValue[$genes[$maximalOrder]]=LeonStat::AVERAGE($eachGeneAverageM);	 
		  
	  }
 
   
   //delete the most unstable gene data
  // unset($genes[$maximalOrder]);
  // unset($genormNormalizedDataByGene[$maximalOrder]);
    
//delete gene and its data 
$data22=array();
$data33=array();
$sum=0;
for($i=0;$i<$newGenesNumber;$i++)
{
 if($i<>$maximalOrder)
 { 
  $data33[$sum]=$genes[$i];
  for($j=0;$j<$this->treatmentNumber;$j++)
  {
   $data22[$sum][$j]=$genormNormalizedDataByGene[$i][$j];
  } 
   $sum++;
 }
}
$genormNormalizedDataByGene=$data22;
$genes=$data33;
unset($data22);
unset($data33);
$newGenesNumber--;    
 }
   
  
  
  $this->genormResult=$lastGeneStabilityValue;
  

}

function getMatchOrderNo($data) #get the order of maximal digit in $data
{
  $m=$data[0];
  $order=0;
  for($i=0;$i < count($data);$i++)
  {
    if($m<$data[$i])
	{
	  $m=$data[$i];
	  $order=$i;
	}
  }
  return $order;

}
 
 
}
?>