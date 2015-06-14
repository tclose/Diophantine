<!-- <html> -->
<!-- <body bgcolor="FFFFCCC"> -->
<?php
include ("check_input.php");
include ("library.php");
include ("lllhermite_.php");
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
printmat1 ( $mat, $rows, $cols );
// echo "<br>\n";
$transposed = transpose ( $mat, $rows, $cols );
$m1 = "1";
$n1 = "1";
$n = $rows;

// axb ( $transposed, $m, $rows, $m1, $n1 );
// print flagcol($transposed, $m, $n);
print reduce2($k, $i, $m, $n, $D)
// minus($j, $m, $L)
// swap2($k, $m, $n)
// zero_row_test($matrix, $n, $i)
// shortest_distance($A, $m, $n)
// cholesky($A, $m)
// gram($A, $m, $n)
// introot($a, $b, $c, $d)
// egcd($p, $q)
// lnearint($a, $b)
// ratior($a, $b, $c, $d)
// multr($a, $b, $c, $d)
// subr($a, $b, $c, $d)
// addr($a, $b, $c, $d)
// comparer($a, $b, $c, $d)
// lcasvector($A, $X, $m, $n)
// exit ();

// echo "<br>\n";
// print "<p>\n";
// flush ();
// print "<a href=\"./axb.html\">Return to main page</a><br>\n";
flush ();
// ?>
<!-- </body> -->
<!-- </html> -->
