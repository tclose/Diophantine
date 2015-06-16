<!-- <html> -->
<!-- <body bgcolor="FFFFCCC"> -->
<?php
include ("check_input.php");
include ("library.php");
include ("lllhermite_.php");
error_reporting(E_ALL ^ E_DEPRECATED);
global $transposed;

// $matrix = "0 1 0 0 0 -1 0 -1 -1 1 1 0 0 1 0 2 -2 0 -3 -4 1 -4 -2 3 2 0 0  4 1 -3 0 0 0 4 0 3 3 -3 -3 0 0  -3 0 -1 1 1 0 2 0 2 2 -2 -2 0 0  -1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0";
$matrix = "0 0 1 -1 -1 0 0 0 2 0 0 2 -4 -4 0 0 0 6 0 1 -3 4 3 0 0 0 -7 1 0 -1 2 2 0 0 0 -3 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 0 0";
$a = split ( '[ ]+', $matrix );
$rows = 7;
$cols = count ( $a ) / 7;
$ii = "0";
for($i = "1"; le ( $i, $rows ); $i = bcadd ( $i, "1" )) {
	for($j = "1"; le ( $j, $cols ); $j = bcadd ( $j, "1" )) {
		$k = bcadd ( $ii, $j );
		$k = bcsub ( $k, "1" );
		$mat[$i][$j] = $a[$k];
	}
	$ii = bcadd ( $cols, $ii );
}

$m = bcsub ( $cols, "1" );
$t = test_zeromat ( $mat, $rows, $m );
// echo "Augmented matrix [A|B]=";
// printmat1 ( $mat, $rows, $cols );
// echo "<br>\n";
$transposed = transpose ( $mat, $rows, $cols );
$m1 = "1";
$n1 = "1";
$n = $rows;

$a = split ( '[ ]+', "-2 -1 9 1 2");
$b = split ( '[ ]+', "4 2 -5 7 -6");
$c = split ( '[ ]+', "8 -1 11 -1 5");
$d = split ( '[ ]+', "-5 1 3 -2 1");

// // Scalar functions
// for($i = "0"; lt ( $i, "5" ); $i = bcadd ( $i, "1" )) {
// //  	print $a[$i];
// //  	print $b[$i];
// //  	print $c[$i];
// //  	print $d[$i];
//   print "-------------- i = $i --------------\n";
//  	print "introot($a[$i], $b[$i], $c[$i], $d[$i]): " . introot(abs($a[$i]), abs($b[$i]), $c[$i], $d[$i]) . "\n";
//  	$out = egcd($a[$i], $b[$i]); 	
// 	print "egcd($a[$i], $b[$i]): " . $out . ", " . $multiplier1 . ", " . $multiplier2 . "\n";
// 	print "lnearint($a[$i], $b[$i]): " . lnearint($a[$i], $b[$i]) . "\n";
// 	ratior($a[$i], $b[$i], $c[$i], $d[$i]);
// 	print "ratior($a[$i], $b[$i], $c[$i], $d[$i]): " . $rationum . ", " . $ratioden . "\n";
// 	multr($a[$i], $b[$i], $c[$i], $d[$i]);	
// 	print "multr($a[$i], $b[$i], $c[$i], $d[$i]): " . $multnum . ", " . $multden . "\n";
// 	subr($a[$i], $b[$i], $c[$i], $d[$i]);
// 	print "subr($a[$i], $b[$i], $c[$i], $d[$i]): " . $subnum . ", " . $subden . "\n";
// 	addr($a[$i], $b[$i], $c[$i], $d[$i]);
// 	print "addr($a[$i], $b[$i], $c[$i], $d[$i]): " . $addnum . ", " . $addden . "\n";
// 	print "comparer($a[$i], $b[$i], $c[$i], $d[$i]): " . comparer($a[$i], abs($b[$i]), $c[$i], abs($d[$i])) . "\n";
// }


// 

// // Scalar functions
// for($i = "0"; lt ( $i, "5" ); $i = bcadd ( $i, "1" )) {
// 	//  	print $a[$i];
// 	//  	print $b[$i];
// 	//  	print $c[$i];
// 	//  	print $d[$i];
// 	print "-------------- i = $i --------------\n";
// 	print "introot($a[$i], $b[$i], $c[$i], $d[$i]): " . introot(abs($a[$i]), abs($b[$i]), $c[$i], $d[$i]) . "\n";
// 	$out = egcd($a[$i], $b[$i]);
// 	print "egcd($a[$i], $b[$i]): " . $out . ", " . $multiplier1 . ", " . $multiplier2 . "\n";
// 	print "lnearint($a[$i], $b[$i]): " . lnearint($a[$i], $b[$i]) . "\n";
// 	ratior($a[$i], $b[$i], $c[$i], $d[$i]);
// 	print "ratior($a[$i], $b[$i], $c[$i], $d[$i]): " . $rationum . ", " . $ratioden . "\n";
// 	multr($a[$i], $b[$i], $c[$i], $d[$i]);
// 	print "multr($a[$i], $b[$i], $c[$i], $d[$i]): " . $multnum . ", " . $multden . "\n";
// 	subr($a[$i], $b[$i], $c[$i], $d[$i]);
// 	print "subr($a[$i], $b[$i], $c[$i], $d[$i]): " . $subnum . ", " . $subden . "\n";
// 	addr($a[$i], $b[$i], $c[$i], $d[$i]);
// 	print "addr($a[$i], $b[$i], $c[$i], $d[$i]): " . $addnum . ", " . $addden . "\n";
// 	print "comparer($a[$i], $b[$i], $c[$i], $d[$i]): " . comparer($a[$i], abs($b[$i]), $c[$i], abs($d[$i])) . "\n";
// }


function printnp($matrix,$m,$n){
	$neg = 0;
	$maxlen = 0;
	for($i="1";$i<=$m;$i=bcadd($i,"1")){
		for($j="1";$j<=$n;$j=bcadd($j,"1")){
			if($matrix[$i][$j] < 0){
				$neg = 1;
			}
			if(len($matrix[$i][$j])>$maxlen){
				$maxlen = len($matrix[$i][$j]);
			}
		}
	}
	for($i="1";$i<=$m;$i=bcadd($i,"1")){
		if ($i == 1){
			print "[[";
		}else{
			print " [";
		}	
		for($j="1";$j<=$n;$j=bcadd($j,"1")){
			if($j != "1"){
				print " ";
			}
			$l = len($matrix[$i][$j]);
			if($l==0){
				$l = 1;
			}
			for($space=0;$space<($maxlen - $l);$space++){
				print " ";
			}
			if($matrix[$i][$j] >= 0 && $neg){
				print " ";
			}
			print $matrix[$i][$j] ;
		}
		if ($i == $m){
			print "]]\n";
		}else {
			print "]\n";
		}
	}
}

function printnparray($a,$start,$finish){
	$neg = 0;
	$maxlen = 0;
	for($i=$start;$i<$finish;$i=bcadd($i,"1")){
		if($a[$i] < 0){
			$neg = 1;
		}
		if(len($a[$i])>$maxlen){
			$maxlen = len($a[$i]);
		}
	}
	echo "[";
	for($i=$start;lt($i,$finish);$i=bcadd($i,"1")){
		if($i!=$start){
			print " ";
		}
		if($a[$i] >= 0 && $neg){
			print " ";
		}
		$l = len($a[$i]);
		if($l==0){
			$l = 1;
		}
		for($space=0;$space<($maxlen - $l);$space++){
			print " ";
		}
		print "$a[$i]";
	}
	print "]\n";
	return;
}


function print_all($m, $n, $m1, $n1) {
	global $B;
	global $L;
	global $A;
	global $D;
	print "A: \n";
	printnp($A, $m, $n);
	print "B: \n";
	printnp($B, $m, $m);
	print "L: \n";
	$neg = 0;
	for($i="2";$i<=$m;$i=bcadd($i,"1")){
		for($j="1";$j<$i;$j=bcadd($j,"1")){
			if($L[$i][$j] < 0){
				$neg = 1;
			}
		}
	}
	print "[[";
	for ($s="1";le($s, $m); $s=bcadd($s,"1")){
		if($s != "1"){
			print " ";
		}		
		if($neg){
			print " ";
		}
		print "0";
	}
	print "]\n";
	for($r="2";le($r,$m);$r=bcadd($r,"1")){
		print " [";
		for($s="1";lt($s,$r);$s=bcadd($s,"1")){
			if($s != "1"){
				print " ";
			}
			if($L[$r][$s] >= 0 && $neg){
				print " ";
			}			
			print $L[$r][$s];
		}
		for ($s=$r;le($s, $m); $s=bcadd($s,"1")){
			if($neg){
				print " ";
			}
			print " 0";
		}
		if ($r == $m){
			print "]]\n";
		}else {
			print "]\n";
		}
	}
	print "D: \n";
	printnparray($D, 0, $m);
}


global $col1;
global $col2;
global $nplus1;
global $B;
global $L;
global $A;
global $D;
global $hnf;
global $unimodular_matrix;
global $rank;
global $lcv;
global $choleskynum;
global $choleskyden;
global $G;

$arrays[0] = "-3 -2 -4 -3 -1 0 -3 0 1 3 3 -4 3 -1 3 -2 -4 -2 -1 0 2 1 0 -2 -4 3 3 -4 0 0 4 4 3 -4 2 4 1 0 -3 -2 1 2 2 1 -2 0 2 0 -3 -1 4 0 -2 -1 0 4 4 2 0 0 -4 1 -4 4 -4 0 -2 3 4 4";
$arrays[1] = "-1 0 0 -3 -3 -3 4 0 1 4 0 -2 -2 4 2 -4 0 -3 -4 2 -2 3 1 -4 2 -1 1 -4 0 1 4 -3 -2 2 -1 1 -4 -2 4 1 -3 -2 -1 -3 0 -4 1 -3 3 1 -4 -1 0 -3 0 0 3 3 -4 0 1 1 -1 -3 2 2 -3 3 2 2";
$arrays[2] = "-1 -2 3 0 4 0 -4 -3 4 -2 -3 2 -2 0 -4 3 3 2 0 -4 0 0 0 0 0 0 0 0 0 0 -1 -1 0 0 -1 -1 3 -3 4 2 3 -1 2 0 -4 -3 -1 -1 2 3 -4 0 -1 0 4 0 1 -4 -2 0 4 2 3 0 0 -2 2 -2 -4 1";
$arrays[3] = "-3 3 4 -1 0 -4 -1 -4 2 -2 1 2 3 -1 -3 3 -3 -2 1 -2 -4 2 2 -2 -3 -1 -2 -4 0 2 0 -4 -3 -3 1 2 0 -3 1 -1 -1 -1 3 1 1 4 -3 -3 0 2 0 1 -4 1 -3 0 -1 0 1 0 0 0 -2 -2 4 0 4 1 2 0";
$arrays[4] = "4 -3 0 0 0 3 4 -4 0 -3 -4 4 -3 0 -3 -3 2 0 -1 -1 0 3 4 0 2 -2 2 2 0 3 -3 1 0 0 2 0 0 -3 1 1 0 -4 -3 0 0 1 -3 -1 1 0 4 3 2 0 1 -1 0 -2 2 -2 0 0 0 0 0 0 0 0 0 0";

$xs[0] = "-99 3 4 1 0 -3 0 1 4 -2";
$xs[1] = "-99 0 4 4 3 2 4 0 1 4";
$xs[2] = "-99 0 -4 0 1 -4 4 4 -2 3";
$xs[3] = "-99 1 2 -3 3 -4 1 -3 -3 -4";
$xs[4] = "-99 0 -4 -1 -2 -4 0 4 3 -4";

$offset=0;
if($offset>0){
	$end = $offset + 1;
}else{
	$end = 5;
}
// $end = 1;

for($iii = $offset; lt ( $iii, $end ); $iii = bcadd ( $iii, "1" )) {
	$a = split ( '[ ]+', $arrays[$iii] );
	$x = split ( '[ ]+', $xs[$iii] );
	$rows = 7;
	$cols = count ( $a ) / 7;
	$ii = "0";
	for($i = "1"; le ( $i, $rows ); $i = bcadd ( $i, "1" )) {
		for($j = "1"; le ( $j, $cols ); $j = bcadd ( $j, "1" )) {
			$k = bcadd ( $ii, $j );
			$k = bcsub ( $k, "1" );
			$mat[$i][$j] = $a[$k];
		}
		$ii = bcadd ( $cols, $ii );
	}
	
	$m = bcsub ( $cols, "1" );
	$t = test_zeromat ( $mat, $rows, $m );
	$n = $rows;
	
	// echo "Augmented matrix [A|B]=";
	// printmat1 ( $mat, $rows, $cols );
	// echo "<br>\n";
// 	printmatrix($mat, $rows, $cols);
// 	for($i = "1"; le ( $i, $rows ); $i = bcadd ( $i, "1" )) {
// 		for($j = "1"; le ( $j, $cols ); $j = bcadd ( $j, "1" )) {
// 			$transposed[$j][$i] = $mat[$i][$j];
// 			print $transposed[$j][$i];			
// 		}
// 	}
	print "\n\n-------- $iii ---------\n";
 	$transposed = transpose ( $mat, $rows, $cols );
	axb_header ( $transposed, $m, $n, $m1, $n1 );
	$mplus1 = $m + 1;
	$nplus1 = $n + 1;
	lllhermite($G, $mplus1, $nplus1, $m1, $n1);	
// 	print "llhermite: " . $hnf . ", " . $unimodular_matrix . ", " . $rank . "\n";
	print_all($mplus1, $nplus1, $m1, $n1);
// 	print "flagcol(transposed, m, n): " . flagcol($A, $m, $n) . "\n";
	$k = 4;
	$i = $k - 1;
	$j = 4;
// 	swap2(4, $m, $n);
// 	print "swap2($k, $m, $n): " . "\n";
// 	print_all($m, $n, $m1, $n1);
	
// 	reduce2($k, $i, $mplus1, $nplus1, $D);
// 	print "reduce2($k, $i, $mplus1, $nplus1, D): " . $col1 . ", " . $col2 . "\n";
// 	print_all($mplus1, $nplus1, $m1, $n1);	
// 	minus($j, $n, $A);
// 	print "minus($j, $m, L):\n";
// 	print_all($mplus1, $nplus1, $m1, $n1);	
// 	print "zero_row_test(A, $n, $k): " . zero_row_test($A, $n, $k) . "\n";
// 	$X = gram($A, $mplus1, $nplus1);
// 	print "gram(A, $mplus1, $nplus1):\n";
// 	printnp($X, $mplus1, $mplus1);
 	lcasvector($A, $x, $m, $n);
	print "lcasvector(A, x, m, nplus1): ";
	printnparray($lcv, 1, $nplus1);
// 	cholesky($A, $m);
// 	print "cholesky($A, $m): ";
//  shortest_distance($A, $m, $n);
// 	print "shortest_distance($A, $m, $n): ";
}


// exit ();

// echo "<br>\n";
// print "<p>\n";
// flush ();
// print "<a href=\"./axb.html\">Return to main page</a><br>\n";
flush ();
// ?>
<!-- </body> -->
<!-- </html> -->
