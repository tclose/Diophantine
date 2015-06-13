<?php
/* php program check_input.php */
function check_decimal($Mvalue){
if ( preg_match( '/^0$/', $Mvalue) ){
	$value="1";
	return ($value);
} 

if ( preg_match( '/^(-[1-9])[0-9]*$/', $Mvalue) ) {
	$value="1";
	return($value);
}
elseif ( preg_match( '/^([1-9])[0-9]*$/', $Mvalue) ) {
	$value="1";
	return($value);
} else {
	print "The value \"$Mvalue\" is badly formatted!<br>\n";
	$value="0";
	return($value);
}
}
?>
