<html>
<body bgcolor="FFFFCCC">
<?php
include("check_input.php");
include("library.php");
include("lllhermite_.php");
global $transposed;

$rows=7;
$cols=13;
$matrix="0 1 0 0 0 -1 0 -1 -1 1 1 0 0 1 0 2 -2 0 -3 -4 1 -4 -2 3 2 0 0  4 1 -3 0 0 0 4 0 3 3 -3 -3 0 0  -3 0 -1 1 1 0 2 0 2 2 -2 -2 0 0  -1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0";
$a=split('[ ]+',$matrix);
$t=count($a);
if(le($t,"1")){
	print "number of entries is less than or equal to 1<br>\n";
        print "<a href=\"./axb.html\">Return to main page</a><br>\n";
	flush();
        exit;
}
$cols=bcadd($cols,"1");/* Here $cols is the number of columns of the augmented matrix */
$size=bcmul($rows,$cols);
if(lt($t,$size)){
	print "number $t of entries is less than m &times; n = $size<br>\n";
        print "<a href=\"./axb.html\">Return to main page</a><br>\n";
	flush();
        exit;
}
if(gt($t,$size)){
	print "number $t of entries is greater than m &times; n = $size<br>\n";
        print "<a href=\"./axb.html\">Return to main page</a><br>\n";
	flush();
        exit;
}

	$flag="0";
	for($i="0";lt($i,$t);$i=bcadd($i,"1")){
		$check=check_decimal($a[$i]);
		if (ezero($check)){
			print "<p>\n";
			flush();
			$flag="1";
			break;
		}
	}
	if(ezero($flag)){
                $ii="0";
		for($i="1";le($i,$rows);$i=bcadd($i,"1")){
		    for($j="1";le($j,$cols);$j=bcadd($j,"1")){
                        $k=bcadd($ii,$j);
                        $k=bcsub($k,"1");
                        $mat[$i][$j]=$a[$k];
                    }
                    $ii=bcadd($cols,$ii);
		}
                $m=bcsub($cols,"1");
                $t=test_zeromat($mat,$rows,$m);
                if(eq($t,"1")){
                   echo "Coeffficient matrix is the zero matrix<br>\n";
                }else{
                   echo "Augmented matrix [A|B]=";
                   printmat1($mat,$rows,$cols);
                   echo "<br>\n";
                   $transposed=transpose($mat,$rows,$cols);
                   $m1="1";
                   $n1="1";
                   axb($transposed,$m,$rows,$m1,$n1);
                   echo "<br>\n";
                }
	}
print "<p>\n";
flush();
print "<a href=\"./axb.html\">Return to main page</a><br>\n";
flush();
?>
</body>
</html>
