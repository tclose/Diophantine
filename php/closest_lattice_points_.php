<?php
function cholesky($A,$m){ // A is positive definite mxm
  global $rationum;
  global $ratioden;
  global $multnum;
  global $multden;
  global $subnum;
  global $subden;
  global $choleskynum;
  global $choleskyden;

  $Qnum=$A;
  for($i=1;le($i,$m);$i=bcadd($i,"1")){
      for($j=1;le($j,$m);$j=bcadd($j,"1")){
          $Qden[$i][$j]="1";
      }
  }
  for($i=1;lt($i,$m);$i=bcadd($i,"1")){
     $iplus1=bcadd($i,"1");
     for($j=$iplus1;le($j,$m);$j=bcadd($j,"1")){
         $Qnum[$j][$i]=$Qnum[$i][$j];
         $Qden[$j][$i]=$Qden[$i][$j];
         ratior($Qnum[$i][$j],$Qden[$i][$j],$Qnum[$i][$i],$Qden[$i][$i]);
         $Qnum[$i][$j]=$rationum;
         $Qden[$i][$j]=$ratioden;
     }
     for($k=$iplus1;le($k,$m);$k=bcadd($k,"1")){
         for($l=$k;le($l,$m);$l=bcadd($l,"1")){
               multr($Qnum[$k][$i],$Qden[$k][$i],$Qnum[$i][$l],$Qden[$i][$l]);
               $t2num=$multnum;
               $t2den=$multden;
               subr($Qnum[$k][$l],$Qden[$k][$l],$t2num,$t2den);
               $Qnum[$k][$l]=$subnum;
               $Qden[$k][$l]=$subden;
         }
     }
  }
  for($i=1;le($i,$m);$i=bcadd($i,"1")){
      for($j=1;le($j,$m);$j=bcadd($j,"1")){
          $choleskynum[$i][$j]=$Qnum[$i][$j];
          $choleskyden[$i][$j]=$Qden[$i][$j];
      }
  }
  return;
}

function gram($A,$m,$n){
   for($i="1";le($i,$m);$i=bcadd($i,"1")){
       for($j="1";le($j,$m);$j=bcadd($j,"1")){
           $B[$i][$j]=dotproduct($A[$i],$A[$j],$n);
       }
   }
   return($B);
}

// With Z=$a/$b, U=$c/$d, returns [sqrt($a/$b)+$c/$d]. First ANSWER = [sqrt(Z)] + [U]. One then 
// tests if Z < ([sqrt(Z)] + 1 -{U})^2. If this does not hold, ANSWER++.
// For use in fincke_pohst().
function introot($a,$b,$c,$d){
global $subnum;
global $subden;
global $addnum;
global $addden;
global $multnum;
global $multden;
   $y=int($c,$d);
   if(ezero($a)){
      return($y);
   }
   $x=bcdiv($a,$b);
   $x=bcsqrt($x);
   $answer=bcadd($x,$y);
   subr($c,$d,$y,"1");
   subr("1","1",$subnum,$subden);
   addr($x,"1",$subnum,$subden);
   multr($addnum,$addden,$addnum,$addden);
   $t=comparer($multnum,$multden,$a,$b);
   if(le($t,"0")){
      $answer=bcadd($answer,"1");
   }
   return($answer);
}

function shortest_distance($A,$m,$n){
global $choleskynum;
global $choleskyden;
global $multnum;
global $multden;
global $addnum;
global $addden;
global $subnum;
global $subden;
global $rationum;
global $ratioden;
global $lcv;

    $count="0";
    $min_count="0";

    echo "matrix A:";
    printmat1($A,$m,$n);
    echo "<br>\n";
    $nplus1=bcadd($n,"1");
    //for($j="1";le($j,$n);$j=bcadd($j,"1")){
    // $Am[$j]=$A[$m][$j];
  //}
    print "P = A[$m] =  ";printarray($A[$m],$n);print "<br>\n";
   // $lengthj=dotproduct($Am,$Am,$m);

    $mplus1=$m;
    $mminus1=bcsub($m,"1");
    if(gt($mminus1,"1")){
       print "&#8466; is the lattice spanned by the first $mminus1 rows of A<br>\n";
    }else{
       print "&#8466; is the lattice spanned by the first row of A<br>\n";
    }
    for($i="1";le($i,$mminus1);$i=bcadd($i,"1")){  // AA consists of the first m-1 rows of A
        for($j="1";le($j,$n);$j=bcadd($j,"1")){
            $AA[$i][$j]=$A[$i][$j];
        }
    }
    $G=gram($A,$m,$n); 
    $lengthj=$G[$m][$m];
    cholesky($G,$m);
    $Qnum=$choleskynum;
    $Qden=$choleskyden;
    $QQnum=transpose($snum,$m,$m);
    $QQden=transpose($sden,$m,$m);
    $m=bcsub($m,"1");
    for($i="1";le($i,$m);$i=bcadd($i,"1")){// the N vector
        $Nnum[$i]=$Qnum[$i][$mplus1];
        $Nden[$i]=$Qden[$i][$mplus1];
    }

    $Cnum="0";
    $Cden="1";
    for($i="1";le($i,$m);$i=bcadd($i,"1")){
        multr($Nnum[$i],$Nden[$i],$Nnum[$i],$Nden[$i]);
        multr($multnum,$multden,$Qnum[$i][$i],$Qden[$i][$i]);
        addr($Cnum,$Cden,$multnum,$multden);
        $Cnum=$addnum;
        $Cden=$addden;
    }
    $i=$m;
    $Tnum[$m]=$Cnum;
    $Tden[$m]=$Cden;
    $Unum[$m]="0";
    $Uden[$m]="1";
    while(1){
       ratior($Tnum[$i],$Tden[$i],$Qnum[$i][$i],$Qden[$i][$i]);
       $Znum=$rationum;
       $Zden=$ratioden;
       subr($Nnum[$i],$Nden[$i],$Unum[$i],$Uden[$i]);
       $UB[$i]=introot($Znum,$Zden,$subnum,$subden);
       subr($Unum[$i],$Uden[$i],$Nnum[$i],$Nden[$i]);
       $temp2=introot($Znum,$Zden,$subnum,$subden);
       $temp3=bcminus($temp2);
       $x[$i]=bcsub($temp3,"1");
       while("1"){
          $x[$i]=bcadd($x[$i],"1");
          if(le($x[$i],$UB[$i])){
              if(eq($i,"1")){
                   //$s=printlc($A,$x,$m);
                   lcasvector($AA,$x,$m,$n);
                   $count=bcadd($count,"1");
               //  print "X[$count]=";printarray($x,$m);
                   $lcva[$count]=$lcv;
               //  print "lcv[$count]=";printarray($lcv,$n);
               //  print "<br>\n";
                   $coord[$count]=$x;
                   for($k="1";le($k,$n);$k=bcadd($k,"1")){
                       $temp=$A[$mplus1][$k];
                       $multiplier_vector[$count][$k]=bcsub($temp,$lcv[$k]);
                   }
                   $l=lengthsquared($multiplier_vector[$count],$n);
                   $multiplier_vector[$count][$nplus1]=$l;
                // print "P-X[$count]=";printarray($multiplier_vector[$count],$n);print": $l<br>\n";
                   $lengtharray[$count]=$l;
                   continue;
              }else{
                $i=bcsub($i,"1");
                // now update $U[$i];
                $sumnum="0";
                $sumden="1";
                $iplus1=bcadd($i,"1");
                for($j=$iplus1;le($j,$m);$j=bcadd($j,"1")){
                    multr($Qnum[$i][$j],$Qden[$i][$j],$x[$j],"1");
                    addr($sumnum,$sumden,$multnum,$multden);
                    $sumnum=$addnum;
                    $sumden=$addden;
                }
                $Unum[$i]=$sumnum;
                $Uden[$i]=$sumden;
                // now update $T[$i];
                addr($x[$iplus1],"1",$Unum[$iplus1],$Uden[$iplus1]);
                subr($addnum,$addden,$Nnum[$iplus1],$Nden[$iplus1]);
                multr($subnum,$subden,$subnum,$subden);
                multr($Qnum[$iplus1][$iplus1],$Qden[$iplus1][$iplus1],$multnum,$multden);
                subr($Tnum[$iplus1],$Tden[$iplus1],$multnum,$multden);
                $Tnum[$i]=$subnum;
                $Tden[$i]=$subden;
                break;
              }
          }else{
             $i=bcadd($i,"1");
             if(gt($i,$m)){
                    echo "Here are the X[k] &isin; &#8466;, P - X[k], ||P-X[k]||<sup>2</sup> such that ||P-X[k]||<sup>2</sup> &le; $lengthj<br>\n";
                    print "<TABLE BORDER=\"1\" CELLSPACING=\"0\">\n";
                    for($k="1";le($k,$count);$k=bcadd($k,"1")){
                           print "<TR>";
                  //       print "<TD ALIGN=\"RIGHT\">";
                  //       minusa($coord[$k],$m); 
                  //       $s=printlc($AA,$coord[$k],$m);
                  //       if(ezero($s)){
                  //          print "b[$mplus1]";
                  //       }else{
                  //          print "+b[$mplus1]";
                  //       }
                  //       print "=";
                  //       print "</TD>";
                           print "<TD ALIGN=\"RIGHT\">";
                           printarray($lcva[$k],$n);
                           print "</TD>";
                           print "<TD ALIGN=\"RIGHT\">";
                           printarray($multiplier_vector[$k],$n);
                           print "</TD>";
                           print "<TD>";
                           print " $lengtharray[$k]";
                           print "</TD>";
                           print "</TR>\n";
                    }
                    print "</TABLE>\n";
                    print "Also<br>\n";
                    $min_length=mina($lengtharray,$count);
                    print "<TABLE BORDER=\"0\" CELLSPACING=\"0\">\n";
                    for($k="1";le($k,$count);$k=bcadd($k,"1")){
                        if(eq($multiplier_vector[$k][$nplus1],$min_length)){
                           print "<TR>";
                           print "<TD ALIGN=\"RIGHT\">";
                           printarray($lcva[$k],$n);
                           print "</TD>";
                           print "</TR>\n";
                           $min_count=bcadd($min_count,"1");
                        }
                    }
                    print "</TABLE>\n";
                    if(eq($min_count,"1")){
                       print " is the closest vector of &#8466; to P, with shortest distance squared $min_length<br>\n";
                    }else{
                       print " are the closest vectors of &#8466; to P, with shortest distance squared $min_length<br>\n";
                    }
                    return;
             }
             continue;
          }
       }
   }
}

function lengthsquared($a,$n){
   $sum="0";
   for($i="1";le($i,$n);$i=bcadd($i,"1")){
      $temp=bcmul($a[$i],$a[$i]);
      $sum=bcadd($sum,$temp);
   }
   return($sum);
}

// returns 0 if all of $a[$i] are zero, otherwise 1.
function zeroarray($a,$n){
  $flag="0";
  for($i="1";le($i,$n);$i=bcadd($i,"1")){
      if(neqzero($a[$i])){
        $flag="1";
        break;
      }
  }
  return($flag);
}

?>