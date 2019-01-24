<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
"http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<title>RefFinder</title>
<link rel="shortcut icon" href="images/favicon.ico">
<meta http-equiv="content-type" content="text/html;charset=utf-8">
<script type="text/javascript" src="js/jquery.min.js"></script>
<link href="css/topmenu.css" type="text/css" rel="stylesheet" />
<link href="css/cotton.css" type="text/css" rel="stylesheet" />
</head>
<body>

<?php 
error_reporting(E_ALL || ~E_NOTICE);
function drawPicture($data,$name, $id)
{
	$mydata='';
	$mycolors='';
	foreach ($data as $k=>$v)
	{
	   $mydata=$mydata."['$k', ".number_format($v,3, '.', '')."],";
	   $mycolors=$mycolors."'#81C714',";
	}
	$mydata=substr($mydata,0, strlen($mydata)-1);
	$mycolors=substr($mycolors,0, strlen($mycolors)-1);
	return "
	<div id=\"$id\"></div><script type=\"text/javascript\">
		var myData = new Array($mydata);
		var colors = [$mycolors];
		var myChart = new JSChart('$id', 'bar');
		myChart.setDataArray(myData);
		myChart.colorizeBars(colors);
		myChart.setTitle('$name');
		myChart.setTitleColor('#000000');
		myChart.setAxisNameX('<== Most stable genes   Least stable genes ==>');
		myChart.setAxisNameY('');
		myChart.setAxisColor('#000000');
		myChart.setAxisNameFontSize(12);
		myChart.setAxisNameColor('#000000');
		myChart.setAxisValuesColor('#000000');
		myChart.setBarValuesColor('#000000');
		myChart.setAxisPaddingTop(60);
		myChart.setAxisPaddingRight(140);
		myChart.setAxisPaddingLeft(150);
		myChart.setAxisPaddingBottom(40);
		myChart.setTextPaddingLeft(105);
		myChart.setTitleFontSize(12);
		myChart.setBarBorderWidth(1);
		myChart.setBarBorderColor('#C4C4C4');
		myChart.setBarSpacingRatio(50);
		myChart.setGrid(false);
		myChart.setSize(830, 321);
		myChart.setBackgroundImage('chart_bg.jpg');
		myChart.draw();
	</script>";

}

?>
<?php
require_once 'class/reference.php';
require_once 'refereceGene/graphicReference.php';
require_once 'class/class.Numerical.php';
require_once 'class/MyStatistics.php';
?>
<?php
	$type=trim($_GET["type"]);
	if($type =="reference")
	{	 
		$data=trim($_POST["data"]);		 			
	}
?>
<link rel="stylesheet" type="text/css" href="css/jqpopup.css"/>
<script type="text/javascript" src="js/jquery.bgiframe.min.js"></script>
 
<script type="text/javascript" src="js/jquery.jqpopup.min.js"></script>
<link rel="stylesheet" href="js/themes/base/jquery.ui.all.css">
<script src="js/ui/jquery.ui.core.js"></script>
<script src="js/ui/jquery.ui.widget.js"></script>
<script src="js/ui/jquery.ui.mouse.js"></script>
<script src="js/ui/jquery.ui.sortable.js"></script>
<script src="js/ui/jquery.ui.tabs.js"></script>
<script src="js/jschart.js"></script>
<div id="helpinfo" style="display:none" title="Information for you"> Just copy your orginial data (CT value) from excel file to the input box<br/>
  Keep data in right data format!<br/>
  Right data format can be referred by trying example.<br/>
  No space is allowed in gene names! </div>
<table  width="950px" align="center">
  <tr>
    <td></br></br><center><b> RefFinder </b></center>
	</br></br>Evaluating Reference Genes Expression<em>>></em> </td>
  </tr>
</table>
<div id="fragment-3">
  <form name="blast" method="post" action="type=reference">
    <table  align='center' border="0" width="950px">
	  <tr>
        <td align='left'><b>RefFinder</b>
          is a user-friendly web-based comprehensive tool developed  for evaluating and screening reference genes from extensive experimental datasets. It integrates the currently available major computational programs (<b>geNorm</b>, <b>Normfinder</b>, <b>BestKeeper</b>, and <b>the comparative ¦¤Ct method</b>) to compare and rank the tested candidate reference genes. Based on the rankings from each program, It assigns an appropriate weight to an individual gene and calculated the geometric mean of their weights for the overall final ranking. Please cite <b><a  href="http://www.ncbi.nlm.nih.gov/pubmed/22290409">F Xie, P Xiao, D Chen, L Xu, B Zhang. 2012. miRDeepFinder: a miRNA analysis tool for deep sequencing of plant small RNAs. Plant molecular biology 80 (1), 75-84.</a></b></br>
        </td>
      </tr>
      <tr>
        <td align='left'><b>Input your data: </b><br/>
          <textarea  id='data' name="data" rows="15" cols="130"><?php echo $data;?></textarea>
        </td>
      </tr>
      <tr>
        <td align='left'><input type='button' value='  Analyze  ' id='referenceAnalyze' >
          &nbsp;&nbsp;&nbsp;&nbsp; <a id="sample" href="#">Try example</a> &nbsp;&nbsp;&nbsp;&nbsp; <a id="removeSeq" href="#">remove data</a> 
          <script language='javascript'>   
		 
$("#sample").click(function(){		
$("#data").val("hBAct	hGAPDH	hSDHA	hTBCA	hTUBA1A	hRNU44	hU6	hRNU48	hRNU47	h18s\n19.3112	22.28325	24.8479	22.9217	24.7194	17.5574	14.46205	19.4794	16.4062	18.99305\n19.16265	22.63935	24.93535	22.8954	24.7734	17.58445	14.4329	19.5376	16.4733	19.33055\n19.14815	22.3895	24.56275	22.51135	24.4619	17.93015	14.4635	19.73925	16.5504	19.5449\n21.81065	24.6102	26.5362	23.36915	27.01725	18.1465	14.4691	20.0296	17.003	19.4468\n21.1704	24.0964	26.02375	23.5005	26.0287	17.7986	15.0001	19.6619	16.43175	19.5778\n23.4701	25.95015	27.0499	24.54845	28.30655	18.60915	16.04265	20.5171	17.3307	20.03305\n19.27045	23.49115	25.0835	22.84805	24.67245	17.7206	14.336	19.8189	16.5204	19.30995\n19.0253	22.8714	24.69045	22.7619	24.47635	17.8875	14.47215	19.87185	16.61655	20.05875\n19.16015	22.9632	24.68925	22.5935	24.49845	18.026	14.72145	19.98605	16.76375	20.56225\n20.23935	24.2292	25.4872	23.1425	25.45795	17.62315	14.73475	19.68395	16.3622	20.12155\n20.6476	23.9726	25.84975	23.4667	25.92005	17.91115	15.0755	19.7871	16.47465	20.0937\n22.8857	26.0722	27.2926	24.5212	27.9778	17.6749	15.2755	19.76915	16.386	20.35435\n19.96615	22.7419	25.27745	22.9304	25.04025	18.04825	14.99655	20.29905	16.9748	20.3836\n20.0786	22.61245	25.4461	22.79935	24.9942	17.74855	14.5316	20.155	16.67935	20.22445\n20.7771	23.82425	25.7362	22.70535	25.11675	16.88815	13.50115	19.1055	15.6059	18.39635\n21.58675	23.7839	26.3449	23.28645	26.0738	18.09565	15.0952	20.4421	17.02225	20.12955\n22.15435	24.16015	26.665	23.533	26.52845	17.21855	14.51215	19.70135	16.02825	18.68725\n24.07285	26.44245	27.4036	24.6452	29.01625	18.28	15.592	20.3794	16.7971	20.24645");
});


$("#removeSeq").click(function(){
$("#data").val("");	
$("#data").focus();					   
});		 
		 
$("#referenceAnalyze").click(function(){					
	if($.trim($("#data").val())=="")
	{
	   alert("data should be not empty");
	   $("#data").focus();
	   return false;
	}
	document.blast.submit();					
});	

   $("#openhelp").click(function () { 
      $("#helpinfo").jqpopup_open(this.id);
   });
   
 
			
	</script>
	
	
	<script type="text/javascript" src="js/highchart/highcharts.js"></script>
	<script type="text/javascript" src="js/highchart/modules/exporting.js"></script>
	
        </td>
      </tr>
      <tr>
        <td style="font:"Times New Roman", Times, serif"><?php	
if($data !="")
{ 
	$Reference=new ReferenceGenes($data);
	asort($Reference->referenceAvg);
	asort($Reference->bkFinalIndex);	
	asort($Reference->normfinderResult);	
	asort($Reference->genormResult);
 
	echo "<table width=90%  id='deltaCT'><tr align='center'><td colspan='".($Reference->genesNumber+1)."'><b>Ranking Order (Better--Good--Average)</b></td> <tr><td>Method</td>";
	$sum=1;
	for($gene =0; $gene< $Reference->genesNumber;$gene++)
	{
	  echo "<td>".$sum++."</td>";
	}   
	echo "</tr>";
	
	echo '<tr><td> <a href="#" id="detactAnchor">Delta CT</a></td>';
   foreach ($Reference->referenceAvg as $k=>$v)
   { 
     echo "<td>".$k."</td>";
   }	
	echo "</tr>";
	
	echo '<tr><td> <a href="#" id="bestkeeperAnchor">  BestKeeper </a></td>';
   foreach ($Reference->bkFinalIndex as $k=>$v)
   { 
     echo "<td>".$k."</td>";
   }	
	echo "</tr>";	
		
		
	echo '<tr><td> <a href="#" id="normdetect">Normfinder</a></td>';
   foreach ($Reference->normfinderResult as $k=>$v)
   { 
     echo "<td>".$k."</td>";
   }	
	echo "</tr>";	
	

	echo '<tr><td> <a href="#" id="Genormdetect">Genorm</a></td>';
	$sum=0;
   foreach ($Reference->genormResult as $k=>$v)
   { 
    $sum++;
    if($sum==2)
	{
	  echo   "<td> </td> <td>".$k."</td>";
	}
	else
	{
     echo "<td>".$k."</td>";
	}
   }	
	echo "</tr>";	

	 $finalGensRank22=array();
	//final ranking
	#delta ct $Reference->referenceAvg 
	#bestkeeper $Reference->bkFinalIndex
	#normfinder $Reference->normfinderResult
	#genorm $Reference->genormResult 
	$sum=0;
	foreach ($Reference->referenceAvg as $k=>$v) 
	{
	  $sum++;
	  $finalGensRank22[$k]=array();
	  array_push($finalGensRank22[$k], $sum);
	}
    $sum=0;
	foreach ($Reference->bkFinalIndex as $k=>$v) 
	{
	  $sum++;
	  array_push($finalGensRank22[$k], $sum);
	}

	$sum=0;
	foreach ($Reference->normfinderResult as $k=>$v) 
	{
	  $sum++;
	  array_push($finalGensRank22[$k], $sum);
	}	
	
	$sum=0;	
	foreach ($Reference->genormResult as $k=>$v) 
	{
	  $sum++;	  
	  if($sum==1)
	  {
	    #hRNU44 | hRNU47
	    $firstTwoGene= preg_split ("/\s+\|\s+/", trim($k));
		array_push($finalGensRank22[$firstTwoGene[0]], $sum);
		array_push($finalGensRank22[$firstTwoGene[1]], $sum);
		$sum++;
	  }
	  else
	  {	  
	     array_push($finalGensRank22[$k], $sum);
	  }
	}
	
	
	//get geomean
    foreach ($Reference->referenceAvg as $k=>$v) 
	{
     $finalGensRank22[$k]=LeonStat::GEOMEAN($finalGensRank22[$k]);
	}	
	
	//last ranking
	asort($finalGensRank22);	
	echo '<tr style="background-color:#D5F5FF;font-weight:bold"><td> <a href="#" id="comprehensiveRank">Recommended comprehensive ranking</a></td>';
	$displayFinalRnak=array();
	$displayFinalRnak=$finalGensRank22;
    foreach ($finalGensRank22 as $k=>$v) 
	{
       echo   "<td>".$k."</td>";
	}
	echo "</tr>";	
	echo "</table>";
	
		
   echo '<div id="tabs">
	  <ul> 	  
	    <li><a href="#tabs-1">Comprehensive Ranking</a></li>  
		<li><a href="#tabs-2">Delta CT</a></li>
		<li><a href="#tabs-3">BestKeeper</a></li>
		<li><a href="#tabs-4">Normfinder</a></li>
		<li><a href="#tabs-5">Genorm</a></li>
		
	</ul>';
	
 	echo '<div id="tabs-1">';
	echo "<table id='comRNAking' width='30%'>    
	<tr><td align='right'>Genes</td><td align='left' >Geomean of ranking values</td></tr>";   
	foreach ($displayFinalRnak as $k=>$v) 
	{
       echo "<tr><td align='right'>$k </td><td align='left' >".number_format($v,2, '.', '')."</td></tr>"; 
	}
	echo "</tr>";	
	echo "</table>";
	 echo drawPicture($displayFinalRnak, "Comprehensive gene stability","containerdeltComprehensive");
	echo '</div>';	
		
	// $myfile = 'class/images/gene.csv';
	// $file_pointer = fopen($myfile,"w");
	
	echo '<div id="tabs-2">';
	echo "<table id='deltaCT' width='30%'> 
	<tr><td align='right'>Genes</td><td align='left' >Average of STDEV</td></tr>";   
	foreach ($Reference->referenceAvg as $k=>$v)
	{ 
	 echo "<tr><td align='right'>$k  </td><td align='left' >".number_format($v,2, '.', '')."</td></tr>"; 
	// fwrite($file_pointer,$k.",".number_format($v,2, '.', '')."\n");
	}
	echo "</table><br/>";	
	
	 echo drawPicture($Reference->referenceAvg, "Gene stability by Delta CT method","containerdeltCT");
	// fclose($file_pointer);
	// echo "<table ><tr><td align='right'><img src='".generateGraphic('class/images/gene.csv')."'  border='0'></td></tr></table>";  
	echo '</div>';



	
	//report best keeper
	echo '<div id="tabs-3">';
	$Reference->reportBestKeeper();  
 echo drawPicture($Reference->bkFinalIndex, "Gene stability by BestKeeper","containerdeltBestKeeper");	
	echo '</div>'; 	

	//normfinder
	echo '<div id="tabs-4">';
	echo "<table id='normfinderTable' >
	 
	<tr><td>Gene name</td><td>Stability value</td><tr/>";
    foreach ($Reference->normfinderResult as $k=>$v)
    { 
     echo "<tr><td>".$k."</td>"."<td>".number_format($v,3, '.', '')."</td></tr>";
    }	 
	echo "</table>";
	 echo drawPicture($Reference->normfinderResult, "Gene stability by normFinder","containerNormfinder");
 
	echo '</div>'; 	
  
  //genorm
	echo '<div id="tabs-5">';
	echo "<table id='GenormdetectTable' >
	 
	<tr><td>Gene name</td><td>Stability value</td><tr/>";
	foreach ($Reference->genormResult as $k=>$v)
	{ 
	 echo "<tr><td>".$k."</td>"."<td>".number_format($v,3, '.', '')."</td></tr>";
	}	

	echo "</table>";
	 echo drawPicture($Reference->genormResult, "Gene stability by Genorm","containerdeltgenorm");
	echo '</div>'; 	
	


}	
?>
        </td>
      </tr>
      <tr>
        <td ><b> References </b> </td>
      </tr>
      <tr>
        <td valign="top"><ol>
          <li><b>BestKeeper</b>: Pfaffl MW, Tichopad A, Prgomet C, Neuvians TP. 2004. Determination of stable housekeeping genes, differentially regulated target genes and sample integrity: BestKeeper--Excel-based tool using pair-wise correlations. Biotechnology letters 26:509-515.</li>
          <li><b>NormFinder</b>: Andersen CL, Jensen JL, Orntoft TF. 2004. Normalization of real-time quantitative reverse transcription-PCR data: a model-based variance estimation approach to identify genes suited for normalization, applied to bladder and colon cancer data sets. Cancer research 64:5245-5250.</li>
          <li><b>Genorm</b>: Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F. 2002. Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome biology 3:RESEARCH0034.</li>
          <li><b>The comparative delta-Ct method</b>: Silver N, Best S, Jiang J, Thein SL. 2006. Selection of housekeeping genes for gene expression studies in human reticulocytes using real-time PCR. BMC molecular biology 7:33.</li>
          <ol>
        </td>
      </tr>
    </table>
  </form>
</div>

<script language="javascript">	
    $(function() {
		$( "#tabs" ).tabs().find( ".ui-tabs-nav" ).sortable({ axis: "x" });
	});
	
	$("#bestkeeper tr, #deltaCT tr, #normfinderTable tr, #GenormdetectTable tr, #comRNAking tr").hover(  		
		function(){
		$(this).addClass("detailedESThover");
		},
		function(){  
		$(this).removeClass("detailedESThover");  
		}   
	); 	



	$("#detactAnchor").click(function(){
	var $tabs = $('#tabs').tabs();  
	 $tabs.tabs('select', 1);  
    return false;	
	});
   
	$("#bestkeeperAnchor").click(function(){
	var $tabs = $('#tabs').tabs();  
	 $tabs.tabs('select', 2);  
    return false;	
	});
	
	$("#normdetect").click(function(){
	var $tabs = $('#tabs').tabs();  
	 $tabs.tabs('select',3);  
    return false;	
	});
	
	$("#Genormdetect").click(function(){
	var $tabs = $('#tabs').tabs();  
	 $tabs.tabs('select', 4);  
    return false;	
	});

	$("#comprehensiveRank").click(function(){
	var $tabs = $('#tabs').tabs();  
	 $tabs.tabs('select', 0);  
    return false;	
	});	
	
		
</script>
<script>

	</script>
 
<!-- footer start -->
	<link rel="stylesheet" href="js/themes/base/jquery.ui.all.css">

	<script src="js/ui/jquery.ui.core.js"></script>
	<script src="js//ui/jquery.ui.widget.js"></script>
	<script src="js/ui/jquery.ui.button.js"></script>
<br>
	<script>
	$(function() {
		$("input[type='button'], input[type='submit'] , input[type='reset']").button();
		//$( "input:submit")
	});
	</script>
<br>
<table align="center" width="950px">
<tr>
<td><div id="foot" >
  <p id="message">Copyright @  <a href="http://www.ecu.edu/cs-cas/biology/zhang_baohong.cfm"> Dr.Zhang's Lab.</a> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Developed by <a href="mailto:fulxie@gmail.com">Fuliang Xie</a>
</div></td>
</tr>
</table>
</body>
</html>