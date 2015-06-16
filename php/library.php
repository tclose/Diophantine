<?php
/* file library.php 
 * This file contains the following functions:
 * mod($a,$b)
 * abs_mod($a,$b)
 * int($a,$b)
 * mpower($a,$b,$c)
 * sign($a)
 * bcmul3($a,$b,$c)
 * bcabs($a)
 * gcd($m,$n)
 * gcda($m,$n)
 * egcd($p,$q)
 * minimum($x,$y)
 * maximum($x,$y)
 * printpoly($a,$n)
 * len($a)
 * ceiling($a,$b)
 * lcm($a,$b)
 * lcma($array,$n)
 * exponential($a,$b)
 * cong($m,$p,$n)
 * cong1($m,$p,$n)
 * chinese2($a,$b,$m,$n)
 * chinesea($a,$m,$n)
 * inverse($a,$m)
 * bcminus($a)
 * powerdd($a,$b,$dd,$n)
 * le($a,$b)
 * ge($a,$b)
 * lt($a,$b)
 * gt($a,$b)
 * neq($a,$b)
 * eq($a,$b)
 * neqzero($a)
 * ezero($a)
 * gezero($a)
 * lezero($a)
 * ltzero($a)
 * gtzero($a)
 * gcd3($a,$b,$c)
 * bezout($a,$b)
 * bezout1($a,$b)
 * parity($a)
 * lnearint($a,$b)
 * bcadd3($a,$b,$c)
 * abpluscd($a,$b,$c,$d)
 * abminuscd($a,$b,$c,$d)
 * aplusbc($a,$b,$c)
 * aminusbc($a,$b,$c)
 * printarray($a,$n)
 * printmatrix($matrix,$m,$n)
 * print_matrix($a,$b,$c,$d) // should change this to printmatrix2x2
 * printmat1($matrix,$m,$n) prints a matrix as a table with entries right justified
 * printmatrix($matrix,$m,$n) 
 * unit_matrix($m)
 * transpose($A,$m,$n)
 * transpose1(&$A,&$m,&$n)
 * row_submatrix($A,$p,$q)
 * col_submatrix($A,$rows,$p,$q)
 * submatrix($A,$rows,$p1,$q1,$p2,$q2)
 * rowiminusqrowj(&$A,$n,$i,$q,$j)
 * coliminusqcolj(&$A,$m,$i,$q,$j)
 * delete_row(&$B,$i,&$m)
 * delete_col(&$B,$j,$m,&$n)
 * swap_rows(&$P,$j,$k)
 * swap_cols(&$P,$m,$j,$k)
 * dotproduct($a,$b,$n)
 * multmat($A,$B,$m,$n,$p)
 * matrixperm(&$A,$a,$m)
 * equalmat($A,$B,$rowsA,$colsA,$rowsB,$colsB)
 * printbinaryform(a,b,c,x,y)
 */

/* the least non-negative remainder when an integer a is divided by a positive
 * integer b.
 * bcmod(a,b)=moda(a,b) if a>=0 or a<0 and b divides a;
 * bcmod(a,b)=mod(a,b)-b if a<0, b>0, a not divisible by b.
 */

function mod($a,$b){
	$c=bcmod($a,$b);
	if($a>="0"){
		 return($c);
	}
	if($c=="0"){
		return("0");
	}
	$temp=bcadd($c,$b);
        return($temp);
}

/* This returns $r=mod($a,$b) if r <= $b/2, otherwise $r-$b.*/
function abs_mod($a,$b){
   $r=mod($a,$b);
   $temp=bcmul("2",$r);
   if(gt($temp,$b)){
      $r=bcsub($r,$b);
   }
   return($r);
}

function int($a,$b){
	if($b<"0"){
	     $a=bcsub(0,$a);
	     $b=bcsub(0,$b);
	}
	$c=bcdiv($a,$b);
	$d=bcmod($a,$b);
	if($d=="0" || $a>"0"){
		return($c);
	}else{
		return(bcsub($c,1));
	}
}

function mpower($a,$b,$c){
	$x=mod($a,$c);
	$y=$b;
	$z="1";
	while($y){
		while(bcmod($y,"2")==0){
			$y=bcdiv($y,"2");
			$x=bcmod(bcmul($x,$x),$c);
		}
		$y=bcsub($y,"1");
		$z=bcmod(bcmul($z,$x),$c);
	}
	return($z);
}

/* sign of an integer a */
/* sign(a)=1,-1,0, according as a>0,a<0,a=0 */

function sign($a){
	if($a>"0"){
		return("1");
	}
	if($a<"0"){
		return("-1");
	}
	return("0");
}

/* signn of an integer a */
/* signn(a)=1,-1, according as a>=0,a<0 */

function signn($a){
	if($a>="0"){
		return("1");
	}else{
		return("-1");
	}
}

/* absolute value */
function bcabs($a){
	if($a>="0"){
		return($a);
	}else{
		$h=bcsub(0,$a);
		return($h);
	}
}

/*  $b=gcd($m,$n) for any integers m and n */
/* Euclid's division algorithm is used. */
/* We use gcd(m,n)=gcd(|m|,|n|) */

function gcd($m,$n){
	$a=bcabs($m);         /* a=r[0] */ 
	if($n=="0"){
	     return($a);
	}
        $b=bcabs($n);         /* b=r[1] */ 
        $c=mod($a,$b);        /* c=r[2]=r[0] mod(r[1]) */
        while($c){
		$a=$b;
                $b=$c;
                $c=mod($a,$b);    /* c=r[j]=r[j-2] mod(r[j-1]) */
        }
	return($b);
}   

function egcd($p,$q){
global $multiplier1;
global $multiplier2;

	if($q=="0"){
		if($p!="0"){
			$s=sign($p);
			if($s=="1"){
				$multiplier1="1";
			}else{
				$multiplier1="-1";
			}
			$multiplier2="0";
			return(bcabs($p));
		}else{
			$multiplier1="0";
			$multiplier2="0";
			return("0");
		}
	}
	$a=$p;
	$b=bcabs($q);
	$c=mod($a,$b);
	$s=sign($q);
	if($c=="0"){
		if($s=="1"){
			$multiplier2="1";
		}else{
			$multiplier2="-1";
		}
		$multiplier1="0";
		return($b);
	}
	$l1="1";
	$k1="0";
	$l2="0";
	$k2="1";
	while($c!="0"){
		$q=int($a,$b);
		$a=$b;
		$b=$c;
		$c=bcmod($a,$b);
		$temp1=bcmul($q,$k1);
		$temp2=bcmul($q,$k2);
		$h1=bcsub($l1,$temp1);
		$h2=bcsub($l2,$temp2);
		$l1=$k1;
		$l2=$k2;
		$k1=$h1;
		$k2=$h2;
	}
	$multiplier1=$k1;
	if($s=="-1"){
		$k2=bcsub("0",$k2);
	}
	$multiplier2=$k2;
	return($b);
}

/* min(x,y) */

function minimum($x,$y){
	if($y<$x){
		return($y);
	}else{
		return($x);
	}
}

/* max(x,y) */

function maximum($x,$y){
	if($y>$x){
		return($y);
	}else{
		return($x);
	}
}

/* php program printpoly.php */
function printpoly($a,$n){
	if(ezero($n)){
		print "$a[0]";
	}else{
	if(neq($a[$n],"1")){
	    if(eq($a[$n],"-1")){
	       print "-";
	    }else{
		print "$a[$n]";
	    }
	}
	if(gt($n,"1")){
	   print "x<sup>$n</sup>";
	}else if(eq($n,"1")){
           print "x";
	}
	$d=bcsub($n,"1");
  	for($i=$d;gezero($i);$i=bcsub($i,"1")){
	    if(neqzero($a[$i])){
		if(gt($a[$i],"1")){
		       print "+$a[$i]";
	        }
		if(eq($a[$i],"1")){
		   if(gtzero($i)){
		       print "+";
		   }
		   if(ezero($i)){
		       	  print "+1";
		   }
	        }
		if(eq($a[$i],"-1")){
		   if(gtzero($i)){
		       print "-";
		   }else{
			print "-1";
		   }
	        }
		if(lt($a[$i],"-1")){
		       print "$a[$i]";
	        }
		if(gt($i,"1")){
	            print "x<sup>$i</sup>";
		}
		if(eq($i,"1")){
		    print "x";
		}
            }
	}
	}
}

/* php program printpolylambda.php */
function printpolylambda($a,$n){
	if(ezero($n)){
		print "$a[0]";
	}else{
	if(neq($a[$n],"1")){
	    if(eq($a[$n],"-1")){
	       print "-";
	    }else{
		print "$a[$n]";
	    }
	}
	if(gt($n,"1")){
	   print "&lambda;<sup>$n</sup>";
	}else if(eq($n,"1")){
           print "&lambda;";
	}
	$d=bcsub($n,"1");
  	for($i=$d;gezero($i);$i=bcsub($i,"1")){
	    if(neqzero($a[$i])){
		if(gt($a[$i],"1")){
		       print "+$a[$i]";
	        }
		if(eq($a[$i],"1")){
		   if(gtzero($i)){
		       print "+";
		   }
		   if(ezero($i)){
		       	  print "+1";
		   }
	        }
		if(eq($a[$i],"-1")){
		   if(gtzero($i)){
		       print "-";
		   }else{
			print "-1";
		   }
	        }
		if(lt($a[$i],"-1")){
		       print "$a[$i]";
	        }
		if(gt($i,"1")){
	            print "&lambda;<sup>$i</sup>";
		}
		if(eq($i,"1")){
		    print "&lambda;";
		}
            }
	}
	}
}

/* If $n > 0, len($n) returns the number of base 10 digits of $n */

function len($n){
	$i="0";
	$x=bcabs($n);
	while($x!="0"){
		$x=int($x,10);
		$i=bcadd($i,1);
	}
	return($i);
}

/* ceiling function */
function ceiling($a,$b){
	$x=int($a,$b);
	if(bccomp(bcmul($b,$x),$a)==0){
           return($x);
        }
	else{
           return(bcadd($x,"1"));
        }
}

/* lcm(a,b) */
function lcm($a,$b){
	$g=gcd($a,$b);
	$h=bcmul($a,$b);
	$k=bcdiv($h,$g);
	return($k);
}

/* lcm(array[0],array[1],...,array[n-1]) */
function lcma($array,$n){
	for($i="0";bccomp($i,$n)<0;$i=bcadd($i,"1")){
          $b[$i]=$array[$i];
	}
	for($i="1";bccomp($i,$n)<0;$i=bcadd($i,"1")){
		$j=bcsub($i,"1");
		$b[$i]=lcm($b[$i],$b[$j]);
	}
	$j=bcsub($i,"1");
	return($b[$j]);
}

/* gcda(array[0],array[1],...,array[n-1]) */
function gcda($array,$n){
	for($i="0";$i<$n;$i=bcadd($i,"1")){
          $b[$i]=$array[$i];
	}
	for($i="1";$i<$n;$i=bcadd($i,"1")){
		$j=bcsub($i,"1");
		$b[$i]=gcd($b[$i],$b[$j]);
	}
	$j=bcsub($i,"1");
	return($b[$j]);
}

/* The bth power of a, where a is an integer, b a positive integer.
 * This performs the same function as bcpow($a,$b).
 */
function exponential($a,$b){
	$x=$a;
	$y=$b;
	$z="1";
	while(bccomp($y,"0")>0){
		while(bccomp(bcmod($y,"2"),"0")==0){
			$y=bcdiv($y,"2");
			$x=bcmul($x,$x);
		}
		$y=bcsub($y,"1");
		$z=bcmul($z,$x);
	}
	return($z);
}

/*  the congruence mx=p(mod n) */

function cong($m,$p,$n){
global $solution;
global $modulus;
global $multiplier1;
	$a=egcd($m,$n);
	$temp=bcmod($p,$a);
	if(bccomp($temp,"0")!=0){
		return("0");
	}
	$b=$multiplier1;
	$y=bcdiv($n,$a);
	$p=int($p,$a);
	$temp1=bcmul($b,$p);
	$solution=mod($temp1,$y);
	$modulus=$y;
	/*for(t=0;t<a;t++){print " ",z+t*y,","}
	print " mod ",n,"\n"*/
	return("1");
}

/*  the congruence mx=p(mod n) slightly modified version of cong($m,$p,$n) */
function cong1($m,$p,$n){
global $modulus;
global $multiplier1;
	$a=egcd($m,$n);
	$b=$multiplier1;
	$y=bcdiv($n,$a);
	$p=int($p,$a);
	$temp1=bcmul($b,$p);
	$solution=mod($temp1,$y);
	$modulus=$y;
	return($solution);
}

/* the Chinese remainder theorem for the congruences x=a(mod m)
 * and x=b(mod n), m>0, n>0, a and b arbitrary integers.
 * The construction of O. Ore, American Mathematical Monthly,
 * vol.59,pp.365-370,1952, is implemented.
 */

function chinese2($a,$b,$m,$n){
global $chinese_modulus;
global $chinese_solution;
global $multiplier1;
global $multiplier2;

	$d = egcd($m,$n);
	if(mod(bcsub($a,$b),$d)!="0"){
		return("0");
	}
	$x= bcdiv($m,$d);$y=bcdiv($n,$d);
	$z=bcdiv(bcmul($m,$n),$d);
	$temp1=bcmul($b,$multiplier1);
	$temp1=bcmul($temp1,$x);
	$temp2=bcmul($a,$multiplier2);
	$temp2=bcmul($temp2,$y);
	$c=mod(bcadd($temp1,$temp2),$z);
	$chinese_modulus=$z;
	$chinese_solution=$c;
	return("1");
}

function chinesea($a,$m,$n){
global $chinese_solution;
global $chinese_modulus;
        $chinese_modulus=$m[0];
        $chinese_solution=$a[0];
        for($i="1";$i<$n;$i=bcadd($i,"1")){
                $y=chinese2($a[$i],$chinese_solution,$m[$i],$chinese_modulus);
                if($y=="0"){
                        return("0");
                }
        }
        return("1");
}

/* Inverse of a (mod m) */
function inverse($a,$m){
	$t=cong1($a,"1",$m);
	return($t);
}

/* lezero($a) returns 1 if $a<="0", "0" otherwise. */
function lezero($a){
   $t=bccomp($a,"0");
   if($t<=0){
      return("1");
   }else{
      return("0");
   }
}

/* gtzero($a) returns 1 if $a>"0", "0" otherwise. */
function gtzero($a){
   $t=bccomp($a,"0");
   if($t>0){
      return("1");
   }else{
      return("0");
   }
}

/* ltzero($a) returns 1 if $a<"0", "0" otherwise. */
function ltzero($a){
   $t=bccomp($a,"0");
   if($t<0){
      return("1");
   }else{
      return("0");
   }
}

/* gezero($a) returns 1 if $a>="0", "0" otherwise. */
function gezero($a){
   $t=bccomp($a,"0");
   if($t>=0){
      return("1");
   }else{
      return("0");
   }
}
/* ezero($a) returns 1 if $a="0", "0" otherwise. */
function ezero($a){
   $t=bccomp($a,"0");
   if($t==0){
      return("1");
   }else{
      return("0");
   }
}

/* neqzero($a) returns 1 if $a != "0", "0" otherwise. */
function neqzero($a){
   $t=bccomp($a,"0");
   if($t!=0){
      return("1");
   }else{
      return("0");
   }
}

/* eq($a,$b) returns 1 if $a=$b, "0" otherwise. */
function eq($a,$b){
   $t=bccomp($a,$b);
   if($t==0){
      return("1");
   }else{
      return("0");
   }
}

/* neq($a,$b) returns 1 if $a != $b, "0" otherwise. */
function neq($a,$b){
   $t=bccomp($a,$b);
   if($t==0){
      return("0");
   }else{
      return("1");
   }
}

/* gt($a,$b) returns 1 if $a > $b, "0" otherwise. */
function gt($a,$b){
   $t=bccomp($a,$b);
   if($t>0){
      return("1");
   }else{
      return("0");
   }
}

/* lt($a,$b) returns 1 if $a < $b, "0" otherwise. */
function lt($a,$b){
   $t=bccomp($a,$b);
   if($t<0){
      return("1");
   }else{
      return("0");
   }
}

/* ge($a,$b) returns 1 if $a >= $b, "0" otherwise. */
function ge($a,$b){
   $t=bccomp($a,$b);
   if($t>=0){
      return("1");
   }else{
      return("0");
   }
}
/* le($a,$b) returns 1 if $a <= $b, "0" otherwise. */
function le($a,$b){
   $t=bccomp($a,$b);
   if($t<="0"){
      return("1");
   }else{
      return("0");
   }
}

function powerdd($a,$b,$dd,$n){
/* (a+bsqrt{dd})^n=zed1+zed2*sqrt{dd} */
global $zed1;
global $zed2;

        $x1=$a;
        $x2=$b;
        $y=$n;
        $zed1="1";
        $zed2="0";
        while(gtzero($y)){
                while(ezero(bcmod($y,"2"))){
                        $y=bcdiv($y,"2");
                        $temp=$x1;
                        $temp1=bcmul($x2,$x2);
                        $temp2=bcmul($x1,$x1);
                        $temp3=bcmul($dd,$temp1);
                        $x1=bcadd($temp2,$temp3);
                        $temp4=bcmul($temp,$x2);
                        $x2=bcmul("2",$temp4);
                }
                $y=bcsub($y,"1");
                $temp=$zed1;
                $temp1=bcmul($zed2,$x2);
                $temp2=bcmul($zed1,$x1);
                $temp3=bcmul($dd,$temp1);
                $zed1=bcadd($temp2,$temp3);
                $temp4=bcmul($temp,$x2);
                $temp5=bcmul($zed2,$x1);
                $zed2=bcadd($temp4,$temp5);
        }
       /* print "(zed1,zed2)=($zed1,$zed2)<br>\n";*/
        return;
}

function bcminus($a){

   $t=bcsub("0",$a);
    return($t);
}

function gcd3($a,$b,$c){

  $t=gcd($a,$b);
  $t=gcd($t,$c);
  return($t);
}

function falling_factorial($m,$n){
    $product="1";
    for($i=$m;le($i,$n);$i=bcadd($i,"1")){
         $product=bcmul($product,$i);
    }
    return($product);
}
function print_matrix($a,$b,$c,$d){
 print "<TABLE BORDER=\"1\" ALIGN=\"CENTER\"\n";
 print "<TR><TD>$a</TD><TD>$b</TD></TR>\n";
 print "<TR><TD>$c</TD><TD>$d</TD></TR>";
 print "</TABLE>";
 return;
}

function sort_array($a,$n){
global $sorted_array;
   $t=bcsub($n,"1");
   for($i="0";lt($i,$t);$i=bcadd($i,"1")){
      $temp1=bcadd($i,"1");
      for($j=$temp1;lt($j,$n);$j=bcadd($j,"1")){
         if(gt($a[$i],$a[$j])){
            $temp=$a[$i];
            $a[$i]=$a[$j];
            $a[$j]=$temp;
         }
      }
   }
   for($i="0";lt($i,$n);$i=bcadd($i,"1")){
      $sorted_array[$i]=$a[$i];
   }
}

/* From Henri Cohen' book, Alg. 1.3.6  13/07/2011
 * This assumes a >=0 and b >= 0.
 * returns d= gcd(a,b) and global variables globalu and globalv,
 * where d = globalu.a + globalv.b.
 */
function bezout($a,$b){
global $globalu;
global $globalv;
   $globalu="1";
   $d=$a;
   if(ezero($b)){
      $globalv="0";
      return($a);
   }else{
      $v1="0";
      $v3=$b;
   }
   while(gtzero($v3)){
      $q=bcdiv($d,$v3);
      $t3=bcmod($d,$v3);
      $temp=bcmul($q,$v1);
      $t1=bcsub($globalu,$temp);
      $globalu=$v1;
      $d=$v3;
      $v1=$t1;
      $v3=$t3;
   }
   $temp=bcmul($a,$globalu);
   $temp=bcsub($d,$temp);
   $globalv=bcdiv($temp,$b);
   return($d);
}

/* Here a and b are any integers.
 * returns d= gcd(a,b) and global variables globalu and globalv,
 * where d = globalu.a + globalv.b.
 */
function bezout1($a,$b){
global $globalu;
global $globalv;

   if(ltzero($a)){
     $absa=bcminus($a);
   }else{
     $absa=$a;
   }
   if(ltzero($b)){
     $absb=bcminus($b);
   }else{
     $absb=$b;
   }
   $d=bezout($absa,$absb);
   $ta=sign($a);
   $tb=sign($b);
   $globalu=bcmul($globalu,$ta);
   $globalv=bcmul($globalv,$tb);
   return($d);
}

function parity($a){
  $r=bcmod($a,"2");
  if(ezero($r)){
    return("0");
  }else{
    return("1");
  }
}

/* left nearest integer 
 * returns y+1/2 if a/b=y+1/2, y integral.
 */
function lnearint($a,$b){
	$y=int($a,$b);
        if(ltzero($b)){
          $a=bcminus($a);
          $b=bcminus($b);
        }
        $x=bcmul($b,$y);
        $z=bcsub($a,$x);
        $z=bcmul("2",$z);
        if(gt($z,$b)){
          $y=bcadd($y,"1");
        }
        return($y);
}

/* lmodd($m,$n) returns $r, $m=$q*$n+$r, -$n/2 < $r <= $n/2 */
function lmodd($m,$n){
	$t=lnearint($m,$n);
	$s=bcmul($n,$t);
	$r=bcsub($m,$s);
	return($r);
}

function mina($a,$n){
    $x=$a[1];
    for($i="2";le($i,$n);$i=bcadd($i,"1")){
         if(lt($a[$i],$x)){
            $x=$a[$i];
         }
    }
    return($x);
 
}

function bcadd3($a,$b,$c){
    $sum=bcadd($a,$b);
    $sum=bcadd($sum,$c);
    return($sum);
}

/* prints a matrix as a table with entries right justified */
function printmat1($matrix,$m,$n){
   print "<TABLE BORDER=\"1\" CELLSPACING=\"0\">\n";
   for($i="1";le($i,$m);$i=bcadd($i,"1")){
       print "<TR>";
       for($j="1";le($j,$n);$j=bcadd($j,"1")){
           $k=$matrix[$i][$j];
          print "<TD ALIGN=\"RIGHT\">$k</TD>";
       }
       print "</TR>\n";
   }
   print "</TABLE>\n";
}

function printmatrix($matrix,$m,$n){
    for($i="1";$i<=$m;$i=bcadd($i,"1")){
       for($j="1";$j<=$n;$j=bcadd($j,"1")){
       print $matrix[$i][$j] ;print " ";
       }
       print "<br>\n";
    }
}

function printmatrix2($matrix,$m,$n){
    for($i="0";lt($i,$m);$i=bcadd($i,"1")){
       for($j="0";lt($j,$n);$j=bcadd($j,"1")){
          print $matrix[$i][$j] ;print " ";
       }
       print "<br>\n";
    }
}
 function printarray($a,$n){
    echo "(";
    for($i="1";lt($i,$n);$i=bcadd($i,"1")){
       print "$a[$i], ";
    }
    print "$a[$n]) ";
    return;
}
 function printarray1($a,$n){
    for($i="0";lt($i,$n);$i=bcadd($i,"1")){
       print "$a[$i], ";
    }
    print "$a[$n]";
    return;
}
 function unit_matrix($m){
   for($i="1";le($i,$m);$i=bcadd($i,"1")){
       for($j="1";le($j,$m);$j=bcadd($j,"1")){
           if(eq($i,$j)){
              $P[$i][$j]="1";
           }else{
              $P[$i][$j]="0";
           }
        }
   }
   return($P);
 }

function transpose($A,$m,$n){
//  global $transposed;
     for($j="1";le($j,$n);$j=bcadd($j,"1")){
         for($i="1";le($i,$m);$i=bcadd($i,"1")){
             $transposed[$j][$i]=$A[$i][$j];
         }
     }
     return($transposed);
}

function transpose1(&$A,&$m,&$n){
     for($j="1";le($j,$n);$j=bcadd($j,"1")){
         for($i="1";le($i,$m);$i=bcadd($i,"1")){
             $transposed[$j][$i]=$A[$i][$j];
         }
     }
     $A=$transposed;
     $temp=$m;
     $m=$n;
     $n=$temp;
}

// creates the submatrix from rows $p to $q.
function row_submatrix($A,$p,$q){
global $new_row_size;
    $r=bcsub($p,"1");
    $s=bcsub($q,$p);
    $s=bcadd($s,"1");
    $new_row_size=$s;
    for($i="1";le($i,$s);$i=bcadd($i,"1")){
        $z=bcadd($i,$r);
        $B[$i]=$A[$z];
    }
    return($B);
}

// creates the submatrix from columns $p to $q.
function col_submatrix($A,$rows,$p,$q){
global $new_col_size;
    $r=bcsub($p,"1");
    $s=bcsub($q,$p);
    $s=bcadd($s,"1");
    $new_col_size=$s;
    for($j="1";le($j,$s);$j=bcadd($j,"1")){
        $z=bcadd($j,$r);
        for($i="1";le($i,$rows);$i=bcadd($i,"1")){
            $B[$i][$j]=$A[$i][$z];
        }
    }
    return($B);
}

// creates the submatrix from rows $p1 to $q1, columns $p2 to $q2.
function submatrix($A,$rows,$p1,$q1,$p2,$q2){
global $new_row_size;
global $new_col_size;
      $B=row_submatrix($A,$p1,$q1);
      $C=col_submatrix($B,$rows,$p2,$q2);
      return($C);
}

// replaces row $i of $A by $q times row $j, updating $A
function rowiminusqrowj(&$A,$n,$i,$q,$j){
    for($k="1";le($k,$n);$k=bcadd($k,"1")){
       $t=bcmul($A[$j][$k],$q);
       $A[$i][$k]=bcsub($A[$i][$k],$t);
    }
    return($A);
}

// replaces column $i of $A by $q times column $j, updating $A
function coliminusqcolj(&$A,$m,$i,$q,$j){
    for($k="1";le($k,$m);$k=bcadd($k,"1")){
       $t=bcmul($A[$k][$j],$q);
       $A[$k][$i]=bcsub($A[$k][$i],$t);
    }
}

function delete_row(&$B,$i,&$m){
   for($l=$i;lt($l,$m);$l=bcadd($l,"1")){
       $lplus1=bcadd($l,"1");
       $temp=$B[$lplus1];
       $B[$l]=$temp;
   }
   $m=bcsub($m,"1");
   return;
}

function delete_col(&$B,$j,$m,&$n){
     transpose1($B,$m,$n);
     delete_row($B,$j,$m);
     transpose1($B,$m,$n);
}

function swap_rows(&$P,$j,$k){
       $temp=$P[$k];
       $P[$k]=$P[$j];
       $P[$j]=$temp;
}

function swap_cols(&$P,$m,$j,$k){
   for($i="1";lt($i,$m);$i=bcadd($i,"1")){
       $temp=$P[$i][$j];
       $P[$i][$j]=$P[$i][$k];
       $P[$i][$k]=$temp;
   }
}

function dotproduct($a,$b,$n){
   $sum="0";
   for($j="1";le($j,$n);$j=bcadd($j,"1")){
       $temp=bcmul($a[$j],$b[$j]);
       $sum=bcadd($sum,$temp);
   }
   return($sum);
}

function multmat($A,$B,$m,$n,$p){
   for($i="1";le($i,$m);$i=bcadd($i,"1")){
       for($k="1";le($k,$p);$k=bcadd($k,"1")){
           $sum="0";
           for($j="1";le($j,$n);$j=bcadd($j,"1")){
               $t=bcmul($A[$i][$j],$B[$j][$k]);
               $sum=bcadd($sum,$t);
           }
           $C[$i][$k]=$sum;
       }
   }
   return($C);
}

// i->a[i] is a permutation of 1,..,m. A is m x m. The rows of A arr permuted.
function matrixperm(&$A,$a,$m){
   for($i="1";le($i,$m);$i=bcadd($i,"1")){
       for($j="1";le($j,$m);$j=bcadd($j,"1")){
           $B[$i]=$A[$a[$i]];
       }
   }
   $A=$B;
}

// outputs 1 or 0 according as $A=$B.
function equalmat($A,$B,$rowsA,$colsA,$rowsB,$colsB){
   if(neq($rowsA,$rowsB) || neq($colsA,$colsB)){
      return("0");
   }
   for($i="1";le($i,$rowsA);$i=bcadd($i,"1")){
       for($j="1";le($j,$colsA);$j=bcadd($j,"1")){
           if(neq($A[$i][$j],$B[$i][$j])){
              return("0");
           }
       }
   }
   return("1");
}
function abpluscd($a,$b,$c,$d){
   $s=bcmul($a,$b);
   $t=bcmul($c,$d);
   $u=bcadd($s,$t);
   return($u);
}
 function abminuscd($a,$b,$c,$d){
   $s=bcmul($a,$b);
   $t=bcmul($c,$d);
   $u=bcsub($s,$t);
   return($u);
}

function aplusbc($a,$b,$c){
   $s=bcmul($b,$c);
   $t=bcadd($a,$s);
   return($t);
}

function aminusbc($a,$b,$c){
   $s=bcmul($b,$c);
   $t=bcsub($a,$s);
   return($t);
}

// returns (a/b)/(c/d)
function ratior($a,$b,$c,$d){
  global $rationum;
  global $ratioden;
  $r=bcmul($a,$d);
  $s=bcmul($b,$c);
  $g=gcd($r,$s);
  if(ltzero($s)){
     $r=bcminus($r);
     $s=bcminus($s);
  }
  $rationum=bcdiv($r,$g);
  $ratioden=bcdiv($s,$g);
  return;
}

// returns (a/b)*(c/d)
function multr($a,$b,$c,$d){
  global $multnum;
  global $multden;
  $r=bcmul($a,$c);
  $s=bcmul($b,$d);
  $g=gcd($r,$s);
  $multnum=bcdiv($r,$g);
  $multden=bcdiv($s,$g);
  return;
}

function subr($a,$b,$c,$d){
  global $subnum;
  global $subden;
  $r=bcmul($a,$d);
  $s=bcmul($b,$c);
  $t=bcsub($r,$s);
  $u=bcmul($b,$d);
  $g=gcd($t,$u);
  $subnum=bcdiv($t,$g);
  $subden=bcdiv($u,$g);
  return;
}

function addr($a,$b,$c,$d){
  global $addnum;
  global $addden;
  $r=bcmul($a,$d);
  $s=bcmul($b,$c);
  $t=bcadd($r,$s);
  $u=bcmul($b,$d);
  $g=gcd($t,$u);
  $addnum=bcdiv($t,$g);
  $addden=bcdiv($u,$g);
  return;
}

// Assumes b>0 and d>0.  Returns -1, 0 or 1 according as a/b <,=,> c/d.
function comparer($a,$b,$c,$d){
  $t=abminuscd($a,$d,$b,$c);
  if(ltzero($t)){
     return("-1");
  }
  if(gtzero($t)){
     return("1");
  }
  return("0");
}

function printlc($A,$X,$m){
 $flag="0";
 $s=zeroarray($X,$m);
 if(ezero($s)){
    return("0");
 }
 for($i="1";le($i,$m);$i=bcadd($i,"1")){
     $t=$X[$i];
     if(ezero($t)){
        continue;
     }
     if(ezero($flag)){
       if(neq($t,"1") && neq($t,"-1")){
          echo "$t";echo "b[$i]";
       }
       if(eq($t,"-1")){
          echo "-b[$i]";
       }
       if(eq($t,"1")){
          echo "b[$i]";
       }
       $flag="1";
     }else{
       if(gtzero($t)){
          if(eq($t,"1")){
             echo "+b[$i]";
          }else{
             echo "+$t";echo "b[$i]";
          }
       }
       if(ltzero($t)){
          if(eq($t,"-1")){
             echo "-b[$i]";
          }else{
             echo "$t";echo "b[$i]";
          }
       }
     }
 }
 return("1");
}

// lcv[j]=X[1]A[1][j]=...+X[m]A[m][j], 1 <= j <= n.
function lcasvector($A,$X,$m,$n){
global $lcv;
// 		printarray($X, $m);
// 		print "\n";
// 		printmatrix($A, $m, $n);
// 		print "\n";
   for($j="1";le($j,$n);$j=bcadd($j,"1")){
      $sum="0";
      for($i="1";le($i,$m);$i=bcadd($i,"1")){
         $t=bcmul($X[$i],$A[$i][$j]);
         $sum=bcadd($sum,$t);
      }
      $lcv[$j]=$sum;
   }
   return;
}
 function minusa(&$a,$n){
   for($j="1";le($j,$n);$j=bcadd($j,"1")){
       $a[$j]=bcminus($a[$j]);
   }
   return;
}

// This returns 1 if $A is the zero matrix, otherwise returns 0.
 function test_zeromat($A,$m,$n){
      for($i="1";le($i,$m);$i=bcadd($i,"1")){
          for($j="1";le($j,$n);$j=bcadd($j,"1")){
              if(neqzero($A[$i][$j])){
                 return("0");
              }
          }
      }
      return("1");
}

function bcmul3($a,$b,$c){
 $temp=bcmul($a,$b);
 $temp=bcmul($temp,$c);
 return($temp);
}

function pparity($e){
   $t=bcmod($e,"2");
   if(eq($t,"1")){
     return("-1");
   }else{
     return("1");
   }
}

function printbinaryform($a,$b,$c,$x,$y){
	if(gt($a,"1") || lt($a,"-1")){
           print"$a&#8203;$x<sup>2</sup>";
        }
        if(eq($a,"1")){
           print"$x<sup>2</sup>";
        }
        if(eq($a,"-1")){
           print"-$x<sup>2</sup>";
        }
	if(neqzero($b)){
	   if(gt($b,"1")){
              print"+$b&#8203;$x$y";
           }
	   if(eq($b,"1")){
              print"+$x$y";
           }
	   if(eq($b,"-1")){
              print"-$x$y";
           }
           if(lt($b,"-1")){
              print"$b&#8203;$x$y";
           }
        }
	if(neqzero($c)){
	   if(gt($c,"1")){
              print"+$c&#8203;$y<sup>2</sup>";
           }
	   if(lt($c,"-1")){
              print"$c&#8203;$y<sup>2</sup>";
           }
           if(eq($c,"1")){
              print"+$y<sup>2</sup>";
           }
           if(eq($c,"-1")){
              print"-$y<sup>2</sup>";
           }
        }
}

function bcmul4($a,$b,$c,$d){
   $temp1=bcmul($a,$b);
   $temp2=bcmul($c,$d);
   $temp=bcmul($temp1,$temp2);
   return($temp);
}

function bcmul5($a,$b,$c,$d,$e){
   $temp1=bcmul4($a,$b,$c,$d);
   $temp=bcmul($temp1,$e);
   return($temp);
}
 
function bcadd4($a,$b,$c,$d){
   $temp1=bcadd3($a,$b,$c);
   $temp=bcadd($temp1,$d);
   return($temp);
}

function bcadd5($a,$b,$c,$d,$e){
   $temp1=bcadd4($a,$b,$c,$d);
   $temp=bcadd($temp1,$e);
   return($temp);
}

function bcadd6($a,$b,$c,$d,$e,$f){
   $temp1=bcadd5($a,$b,$c,$d,$e);
   $temp=bcadd($temp1,$f);
   return($temp);
}

function printaxplusby($a,$x,$b,$y){
   if(eq($a,"1")){
      print "$x";
   }else if(eq($a,"-1")){
      print "-$x";
   }else{
      print "$a$x";
   }
   if(gtzero($b)){
      if(eq($b,"1")){
         print " + $y";
      }else{
         print " + $b$y";
      }
   }
   if(ltzero($b)){
      if(eq($b,"-1")){
         print " - $y";
      }else{
         $minusb=bcminus($b);
         print " - $minusb$y";
      }
   }
   return;
}

function printxplusa($x,$a){
   print $x;
   if(gtzero($a)){
      print " + $a";
   }
   if(ltzero($a)){
      $minusa=bcminus($a);
      print " - $minusa";
   }
   return;
}

function printaxplusbyplusc($a,$x,$b,$y,$c){
   $null=printaxplusby($a,$x,$b,$y);
      if(gtzero($c)){
         print " + $c";
      }
      if(ltzero($c)){
         $minusc=bcminus($c);
         print " - $minusc";
      }
      return;
}

function gcd4($a,$b,$c,$d){

  $t=gcd($a,$b);
  $t=gcd($t,$c);
  $t=gcd($t,$d);
  return($t);
}

function printdXequalsaxplusbyplusc($d,$X,$a,$x,$b,$y,$c){
   $t=gcd4($d,$a,$b,$c);
   $d=bcdiv($d,$t);
   $a=bcdiv($a,$t);
   $b=bcdiv($b,$t);
   $c=bcdiv($c,$t);
   print "$d$X = ";
   $null=printaxplusbyplusc($a,$x,$b,$y,$c);
   return;
}
?>
