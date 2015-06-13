
"""  file library.php 
  This file contains the following functions:
  mod(a,b)
  abs_mod(a,b)
  int(a,b)
  mpower(a,b,c)
  sign(a)
  bcmul3(a,b,c)
  bcabs(a)
  gcd(m,n)
  gcda(m,n)
  egcd(p,q)
  minimum(x,y)
  maximum(x,y)
  printpoly(a,n)
  len(a)
  ceiling(a,b)
  lcm(a,b)
  lcma(array,n)
  exponential(a,b)
  cong(m,p,n)
  cong1(m,p,n)
  chinese2(a,b,m,n)
  chinesea(a,m,n)
  inverse(a,m)
  bcminus(a)
  powerdd(a,b,dd,n)
  le(a,b)
  ge(a,b)
  lt(a,b)
  gt(a,b)
  neq(a,b)
  eq(a,b)
  neqzero(a)
  ezero(a)
  gezero(a)
  lezero(a)
  ltzero(a)
  gtzero(a)
  gcd3(a,b,c)
  bezout(a,b)
  bezout1(a,b)
  parity(a)
  lnearint(a,b)
  bcadd3(a,b,c)
  abpluscd(a,b,c,d)
  abminuscd(a,b,c,d)
  aplusbc(a,b,c)
  aminusbc(a,b,c)
  print[a,n]
  printmatrix(matrix,m,n)
  print_matrix(a,b,c,d) # should change this to printmatrix2x2
  printmat1(matrix,m,n) prints a matrix as a table with entries right justified
  printmatrix(matrix,m,n) 
  unit_matrix(m)
  transpose(A,m,n)
  transpose1(&A,&m,&n)
  row_submatrix(A,p,q)
  col_submatrix(A,rows,p,q)
  submatrix(A,rows,p1,q1,p2,q2)
  rowiminusqrowj(&A,n,i,q,j)
  coliminusqcolj(&A,m,i,q,j)
  delete_row(&B,i,&m)
  delete_col(&B,j,m,&n)
  swap_rows(&P,j,k)
  swap_cols(&P,m,j,k)
  dotproduct(a,b,n)
  multmat(A,B,m,n,p)
  matrixperm(&A,a,m)
  equalmat(A,B,rowsA,colsA,rowsB,colsB)
  printbinaryform(a,b,c,x,y)
  """

"""  the least non-negative remainder when an integer a is divided by a positive
  integer b+ 
  bcmod(a,b)=moda(a,b) if a>=0 or a<0 and b divides a
  bcmod(a,b)=mod(a,b)-b if a<0, b>0, a not divisible by b+ 
  """

def mod(a,b):
	c=bcmod(a,b)
	if a>="0":
		 return(c)
	
	if c=="0":
		return("0")
	
	temp=bcadd(c,b)
        return(temp)


def  abs_mod(a,b):
"""  This returns r=mod(a,b) if r <= b/2, otherwise r-b. """
   r=mod(a,b)
   temp=bcmul("2",r)
   if temp > b:
      r=bcsub(r,b)
   
   return(r)


def int(a,b):
	if b<"0":
	     a=bcsub(0,a)
	     b=bcsub(0,b)
	
	c=bcdiv(a,b)
	d=bcmod(a,b)
	if d=="0" || a>"0":
		return(c)
	else:
		return(bcsub(c,1))
	


def mpower(a,b,c):
	x=mod(a,c)
	y=b
	z="1"
	while(y)
		while(bcmod(y,"2")==0)
			y=bcdiv(y,"2")
			x=bcmod(bcmul(x,x),c)
		
		y=bcsub(y,"1")
		z=bcmod(bcmul(z,x),c)
	
	return(z)


"""  sign of an integer a  """
"""  sign(a)=1,-1,0, according as a>0,a<0,a=0  """

def sign(a):
	if a>"0":
		return("1")
	
	if a<"0":
		return("-1")
	
	return("0")


"""  signn of an integer a  """
"""  signn(a)=1,-1, according as a>=0,a<0  """

def signn(a):
	if a>="0":
		return("1")
	else:
		return("-1")
	


def  bcabs(a):
"""  absolute value  """
	if a>="0":
		return(a)
	else:
		h=bcsub(0,a)
		return(h)
	


"""   b=gcd(m,n) for any integers m and n  """
"""  Euclid's division algorithm is used.  """
"""  We use gcd(m,n)=gcd(|m|,|n|)  """

def gcd(m,n):
	a=bcabs(m)         """  a=r[0]  """ 
	if n=="0":
	     return(a)
	
        b=bcabs(n)         """  b=r[1]  """ 
        c=mod(a,b)        """  c=r[2]=r[0] mod(r[1])  """
        while(c)
		a=b
                b=c
                c=mod(a,b)    """  c=r[j]=r[j-2] mod(r[j-1])  """
        
	return(b)
   

def egcd(p,q):
global multiplier1
global multiplier2

	if q=="0":
		if p!="0":
			s=sign(p)
			if s=="1":
				multiplier1="1"
			else:
				multiplier1="-1"
			
			multiplier2="0"
			return(bcabs(p))
		else:
			multiplier1="0"
			multiplier2="0"
			return("0")
		
	
	a=p
	b=bcabs(q)
	c=mod(a,b)
	s=sign(q)
	if c=="0":
		if s=="1":
			multiplier2="1"
		else:
			multiplier2="-1"
		
		multiplier1="0"
		return(b)
	
	l1="1"
	k1="0"
	l2="0"
	k2="1"
	while(c!="0")
		q=int(a,b)
		a=b
		b=c
		c=bcmod(a,b)
		temp1=bcmul(q,k1)
		temp2=bcmul(q,k2)
		h1=bcsub(l1,temp1)
		h2=bcsub(l2,temp2)
		l1=k1
		l2=k2
		k1=h1
		k2=h2
	
	multiplier1=k1
	if s=="-1":
		k2=bcsub("0",k2)
	
	multiplier2=k2
	return(b)


"""  min(x,y)  """

def minimum(x,y):
	if y<x:
		return(y)
	else:
		return(x)
	


"""  max(x,y)  """

def maximum(x,y):
	if y>x:
		return(y)
	else:
		return(x)
	


def  printpoly(a,n):
"""  php program printpoly.php  """
	if ezero(n):
		print "a[0]"
	else:
	if neq(a[n],"1"):
	    if eq(a[n],"-1"):
	       print "-"
	    else:
		print "a[n]"
	    
	
	if n > "1":
	   print "x<sup>n</sup>"
	else:
           print "x"
	
	d=bcsub(n,"1")
  	for(i=dgezero(i)i=bcsub(i,"1"))
	    if neqzero(a[i]):
		if a[i] > "1":
		       print "+a[i]"
	        
		if eq(a[i],"1"):
		   if gtzero(i):
		       print "+"
		   
		   if ezero(i):
		       	  print "+1"
		   
	        
		if eq(a[i],"-1"):
		   if gtzero(i):
		       print "-"
		   else:
			print "-1"
		   
	        
		if a[i] < "-1":
		       print "a[i]"
	        
		if i > "1":
	            print "x<sup>i</sup>"
		
		if eq(i,"1"):
		    print "x"
		
            
	
	


def  printpolylambda(a,n):
"""  php program printpolylambda.php  """
	if ezero(n):
		print "a[0]"
	else:
	if neq(a[n],"1"):
	    if eq(a[n],"-1"):
	       print "-"
	    else:
		print "a[n]"
	    
	
	if n > "1":
	   print "&lambda<sup>n</sup>"
	else:
           print "&lambda"
	
	d=bcsub(n,"1")
  	for(i=dgezero(i)i=bcsub(i,"1"))
	    if neqzero(a[i]):
		if a[i] > "1":
		       print "+a[i]"
	        
		if eq(a[i],"1"):
		   if gtzero(i):
		       print "+"
		   
		   if ezero(i):
		       	  print "+1"
		   
	        
		if eq(a[i],"-1"):
		   if gtzero(i):
		       print "-"
		   else:
			print "-1"
		   
	        
		if a[i] < "-1":
		       print "a[i]"
	        
		if i > "1":
	            print "&lambda<sup>i</sup>"
		
		if eq(i,"1"):
		    print "&lambda"
		
            
	
	


"""  If n > 0, len(n) returns the number of base 10 digits of n  """

def len(n):
	i="0"
	x=bcabs(n)
	while(x!="0")
		x=int(x,10)
		i=bcadd(i,1)
	
	return(i)


"""  ceiling def  """
function ceiling(a,b):
	x=int(a,b)
	if bccomp(bcmul(b,x),a)==0:
           return(x)
        
	else:
           return(bcadd(x,"1"))
        


def  lcm(a,b):
"""  lcm(a,b)  """
	g=gcd(a,b)
	h=bcmul(a,b)
	k=bcdiv(h,g)
	return(k)


def  lcma(array,n):
"""  lcm(array[0],array[1],...,array[n-1])  """
	for(i="0"bccomp(i,n)<0i=bcadd(i,"1"))
          b[i]=array[i]
	
	for(i="1"bccomp(i,n)<0i=bcadd(i,"1"))
		j=bcsub(i,"1")
		b[i]=lcm(b[i],b[j])
	
	j=bcsub(i,"1")
	return(b[j])


def  gcda(array,n):
"""  gcda(array[0],array[1],...,array[n-1])  """
	for i in xrange(n):
          b[i]=array[i]
	
	for(i="1"i<ni=bcadd(i,"1"))
		j=bcsub(i,"1")
		b[i]=gcd(b[i],b[j])
	
	j=bcsub(i,"1")
	return(b[j])


"""  The bth power of a, where a is an integer, b a positive integer+ 
  This performs the same def as bcpow(a,b):
	x=a
	y=b
	z="1"
	while(bccomp(y,"0")>0)
		while(bccomp(bcmod(y,"2"),"0")==0)
			y=bcdiv(y,"2")
			x=bcmul(x,x)
		
		y=bcsub(y,"1")
		z=bcmul(z,x)
	
	return(z)


"""   the congruence mx=p(mod n)  """

def cong(m,p,n):
global solution
global modulus
global multiplier1
	a=egcd(m,n)
	temp=bcmod(p,a)
	if bccomp(temp,"0")!=0:
		return("0")
	
	b=multiplier1
	y=bcdiv(n,a)
	p=int(p,a)
	temp1=bcmul(b,p)
	solution=mod(temp1,y)
	modulus=y
	""" for(t=0t<at += 1)print " ",z+ty,","
	print " mod ",n,"\n" """
	return("1")


def  cong1(m,p,n):
"""   the congruence mx=p(mod n) slightly modified version of cong(m,p,n)  """
global modulus
global multiplier1
	a=egcd(m,n)
	b=multiplier1
	y=bcdiv(n,a)
	p=int(p,a)
	temp1=bcmul(b,p)
	solution=mod(temp1,y)
	modulus=y
	return(solution)


"""  the Chinese remainder theorem for the congruences x=a(mod m)
  and x=b(mod n), m>0, n>0, a and b arbitrary integers+ 
  The construction of O. Ore, American Mathematical Monthly,
  vol.59,pp.365-370,1952, is implemented+ 
  """

def chinese2(a,b,m,n):
global chinese_modulus
global chinese_solution
global multiplier1
global multiplier2

	d = egcd(m,n)
	if mod(bcsub(a,b),d)!="0":
		return("0")
	
	x= bcdiv(m,d)y=bcdiv(n,d)
	z=bcdiv(bcmul(m,n),d)
	temp1=bcmul(b,multiplier1)
	temp1=bcmul(temp1,x)
	temp2=bcmul(a,multiplier2)
	temp2=bcmul(temp2,y)
	c=mod(bcadd(temp1,temp2),z)
	chinese_modulus=z
	chinese_solution=c
	return("1")


def chinesea(a,m,n):
global chinese_solution
global chinese_modulus
        chinese_modulus=m[0]
        chinese_solution=a[0]
        for(i="1"i<ni=bcadd(i,"1"))
                y=chinese2(a[i],chinese_solution,m[i],chinese_modulus)
                if y=="0":
                        return("0")
                
        
        return("1")


def  inverse(a,m):
"""  Inverse of a (mod m)  """
	t=cong1(a,"1",m)
	return(t)


"""  lezero(a) returns 1 if a<="0", "0" otherwise.  """
def lezero(a):
   t=bccomp(a,"0")
   if t<=0:
      return("1")
   else:
      return("0")
   


"""  gtzero(a) returns 1 if a>"0", "0" otherwise.  """
def gtzero(a):
   t=bccomp(a,"0")
   if t>0:
      return("1")
   else:
      return("0")
   


"""  ltzero(a) returns 1 if a<"0", "0" otherwise.  """
def ltzero(a):
   t=bccomp(a,"0")
   if t<0:
      return("1")
   else:
      return("0")
   


"""  gezero(a) returns 1 if a>="0", "0" otherwise.  """
def gezero(a):
   t=bccomp(a,"0")
   if t>=0:
      return("1")
   else:
      return("0")
   

"""  ezero(a) returns 1 if a="0", "0" otherwise.  """
def ezero(a):
   t=bccomp(a,"0")
   if t==0:
      return("1")
   else:
      return("0")
   


"""  neqzero(a) returns 1 if a != "0", "0" otherwise.  """
def neqzero(a):
   t=bccomp(a,"0")
   if t!=0:
      return("1")
   else:
      return("0")
   


"""  eq(a,b) returns 1 if a=b, "0" otherwise.  """
def eq(a,b):
   t=bccomp(a,b)
   if t==0:
      return("1")
   else:
      return("0")
   


"""  neq(a,b) returns 1 if a != b, "0" otherwise.  """
def neq(a,b):
   t=bccomp(a,b)
   if t==0:
      return("0")
   else:
      return("1")
   


"""  gt(a,b) returns 1 if a > b, "0" otherwise.  """
def gt(a,b):
   t=bccomp(a,b)
   if t>0:
      return("1")
   else:
      return("0")
   


"""  lt(a,b) returns 1 if a < b, "0" otherwise.  """
def lt(a,b):
   t=bccomp(a,b)
   if t<0:
      return("1")
   else:
      return("0")
   


"""  ge(a,b) returns 1 if a >= b, "0" otherwise.  """
def ge(a,b):
   t=bccomp(a,b)
   if t>=0:
      return("1")
   else:
      return("0")
   

"""  le(a,b) returns 1 if a <= b, "0" otherwise.  """
def le(a,b):
   t=bccomp(a,b)
   if t<="0":
      return("1")
   else:
      return("0")
   


def powerdd(a,b,dd,n):
"""  (a+bsqrtdd)^n=zed1+zed2sqrtdd  """
global zed1
global zed2

        x1=a
        x2=b
        y=n
        zed1="1"
        zed2="0"
        while(gtzero(y))
                while(ezero(bcmod(y,"2")))
                        y=bcdiv(y,"2")
                        temp=x1
                        temp1=bcmul(x2,x2)
                        temp2=bcmul(x1,x1)
                        temp3=bcmul(dd,temp1)
                        x1=bcadd(temp2,temp3)
                        temp4=bcmul(temp,x2)
                        x2=bcmul("2",temp4)
                
                y=bcsub(y,"1")
                temp=zed1
                temp1=bcmul(zed2,x2)
                temp2=bcmul(zed1,x1)
                temp3=bcmul(dd,temp1)
                zed1=bcadd(temp2,temp3)
                temp4=bcmul(temp,x2)
                temp5=bcmul(zed2,x1)
                zed2=bcadd(temp4,temp5)
        
       """  print "(zed1,zed2)=(zed1,zed2)<br>\n" """
        return


def bcminus(a):

   t=bcsub("0",a)
    return(t)


def gcd3(a,b,c):

  t=gcd(a,b)
  t=gcd(t,c)
  return(t)


def falling_factorial(m,n):
    product="1"
    for(i=mle(i,n)i=bcadd(i,"1"))
         product=bcmul(product,i)
    
    return(product)

def print_matrix(a,b,c,d):
 print "<TABLE BORDER=\"1\" ALIGN=\"CENTER\"\n"
 print "<TR><TD>a</TD><TD>b</TD></TR>\n"
 print "<TR><TD>c</TD><TD>d</TD></TR>"
 print "</TABLE>"
 return


def sort_[a,n]:
global sorted_array
   t=bcsub(n,"1")
   for i in xrange(t):
      temp1=bcadd(i,"1")
      for j in xrange(temp1, n):
         if a[i] > a[j]:
            temp=a[i]
            a[i]=a[j]
            a[j]=temp
         
      
   
   for i in xrange(n):
      sorted_array[i]=a[i]
   


def  bezout(a,b):
"""  From Henri Cohen' book, Alg. 1.3.6  13/07/2011
  This assumes a >=0 and b >= 0+ 
  returns d= gcd(a,b) and global variables globalu and globalv,
  where d = globalu.a + globalv.b+ 
  """
global globalu
global globalv
   globalu="1"
   d=a
   if ezero(b):
      globalv="0"
      return(a)
   else:
      v1="0"
      v3=b
   
   while(gtzero(v3))
      q=bcdiv(d,v3)
      t3=bcmod(d,v3)
      temp=bcmul(q,v1)
      t1=bcsub(globalu,temp)
      globalu=v1
      d=v3
      v1=t1
      v3=t3
   
   temp=bcmul(a,globalu)
   temp=bcsub(d,temp)
   globalv=bcdiv(temp,b)
   return(d)


def  bezout1(a,b):
"""  Here a and b are any integers+ 
  returns d= gcd(a,b) and global variables globalu and globalv,
  where d = globalu.a + globalv.b+ 
  """
global globalu
global globalv

   if ltzero(a):
     absa=bcminus(a)
   else:
     absa=a
   
   if ltzero(b):
     absb=bcminus(b)
   else:
     absb=b
   
   d=bezout(absa,absb)
   ta=sign(a)
   tb=sign(b)
   globalu=bcmul(globalu,ta)
   globalv=bcmul(globalv,tb)
   return(d)


def parity(a):
  r=bcmod(a,"2")
  if ezero(r):
    return("0")
  else:
    return("1")
  


def  lnearint(a,b):
"""  left nearest integer 
  returns y+1/2 if a/b=y+1/2, y integral+ 
  """
	y=int(a,b)
        if ltzero(b):
          a=bcminus(a)
          b=bcminus(b)
        
        x=bcmul(b,y)
        z=bcsub(a,x)
        z=bcmul("2",z)
        if z > b:
          y=bcadd(y,"1")
        
        return(y)


def  lmodd(m,n):
"""  lmodd(m,n) returns r, m=qn+r, -n/2 < r <= n/2  """
	t=lnearint(m,n)
	s=bcmul(n,t)
	r=bcsub(m,s)
	return(r)


def mina(a,n):
    x=a[1]
    for(i="2"le(i,n)i=bcadd(i,"1"))
         if a[i] < x:
            x=a[i]
         
    
    return(x)
 


def bcadd3(a,b,c):
    sum=bcadd(a,b)
    sum=bcadd(sum,c)
    return(sum)


def  printmat1(matrix,m,n):
"""  prints a matrix as a table with entries right justified  """
   print "<TABLE BORDER=\"1\" CELLSPACING=\"0\">\n"
   for i in xrange(m):
       print "<TR>"
       for j in xrange(n):
           k=matrix[i][j]
          print "<TD ALIGN=\"RIGHT\">k</TD>"
       
       print "</TR>\n"
   
   print "</TABLE>\n"


def printmatrix(matrix,m,n):
    for(i="1"i<=mi=bcadd(i,"1"))
       for(j="1"j<=nj=bcadd(j,"1"))
       print matrix[i][j] print " "
       
       print "<br>\n"
    


def printmatrix2(matrix,m,n):
    for i in xrange(m):
       for j in xrange(n):
          print matrix[i][j] print " "
       
       print "<br>\n"
    

 def print[a,n]:
    print "("
    for i in xrange(n - 1):
       print "a[i], "
    
    print "a[n]) "
    return

 def printarray1(a,n):
    for i in xrange(n):
       print "a[i], "
    
    print "a[n]"
    return

 def unit_matrix(m):
   for i in xrange(m):
       for j in xrange(m):
           if eq(i,j):
              P[i][j]="1"
           else:
              P[i][j]="0"
           
        
   
   return(P)
 

def transpose(A,m,n):
#global transposed
     for j in xrange(n):
         for i in xrange(m):
             transposed[j][i]=A[i][j]
         
     
     return(transposed)


def transpose1(&A,&m,&n):
     for j in xrange(n):
         for i in xrange(m):
             transposed[j][i]=A[i][j]
         
     
     A=transposed
     temp=m
     m=n
     n=temp


# creates the submatrix from rows p to q+ 
def row_submatrix(A,p,q):
global new_row_size
    r=bcsub(p,"1")
    s=bcsub(q,p)
    s=bcadd(s,"1")
    new_row_size=s
    for i in xrange(s):
        z=bcadd(i,r)
        B[i]=A[z]
    
    return(B)


# creates the submatrix from columns p to q+ 
def col_submatrix(A,rows,p,q):
global new_col_size
    r=bcsub(p,"1")
    s=bcsub(q,p)
    s=bcadd(s,"1")
    new_col_size=s
    for j in xrange(s):
        z=bcadd(j,r)
        for i in xrange(rows):
            B[i][j]=A[i][z]
        
    
    return(B)


# creates the submatrix from rows p1 to q1, columns p2 to q2+ 
def submatrix(A,rows,p1,q1,p2,q2):
global new_row_size
global new_col_size
      B=row_submatrix(A,p1,q1)
      C=col_submatrix(B,rows,p2,q2)
      return(C)


# replaces row i of A by q times row j, updating A
def rowiminusqrowj(&A,n,i,q,j):
    for k in xrange(n):
       t=bcmul(A[j][k],q)
       A[i][k]=bcsub(A[i][k],t)
    
    return(A)


# replaces column i of A by q times column j, updating A
def coliminusqcolj(&A,m,i,q,j):
    for k in xrange(m):
       t=bcmul(A[k][j],q)
       A[k][i]=bcsub(A[k][i],t)
    


def delete_row(&B,i,&m):
   for l in xrange(i, m):
       lplus1=bcadd(l,"1")
       temp=B[lplus1]
       B[l]=temp
   
   m=bcsub(m,"1")
   return


def delete_col(&B,j,m,&n):
     transpose1(B,m,n)
     delete_row(B,j,m)
     transpose1(B,m,n)


def swap_rows(&P,j,k):
       temp=P[k]
       P[k]=P[j]
       P[j]=temp


def swap_cols(&P,m,j,k):
   for i in xrange(m - 1):
       temp=P[i][j]
       P[i][j]=P[i][k]
       P[i][k]=temp
   


def dotproduct(a,b,n):
   sum="0"
   for j in xrange(n):
       temp=bcmul(a[j],b[j])
       sum=bcadd(sum,temp)
   
   return(sum)


def multmat(A,B,m,n,p):
   for i in xrange(m):
       for k in xrange(p):
           sum="0"
           for j in xrange(n):
               t=bcmul(A[i][j],B[j][k])
               sum=bcadd(sum,t)
           
           C[i][k]=sum
       
   
   return(C)


# i.a[i] is a permutation of 1,..,m. A is m x m. The rows of A arr permuted+ 
def matrixperm(&A,a,m):
   for i in xrange(m):
       for j in xrange(m):
           B[i]=A[a[i]]
       
   
   A=B


# outputs 1 or 0 according as A=B+ 
def equalmat(A,B,rowsA,colsA,rowsB,colsB):
   if neq(rowsA,rowsB) || neq(colsA,colsB):
      return("0")
   
   for i in xrange(rowsA):
       for j in xrange(colsA):
           if neq(A[i][j],B[i][j]):
              return("0")
           
       
   
   return("1")

def abpluscd(a,b,c,d):
   s=bcmul(a,b)
   t=bcmul(c,d)
   u=bcadd(s,t)
   return(u)

 def abminuscd(a,b,c,d):
   s=bcmul(a,b)
   t=bcmul(c,d)
   u=bcsub(s,t)
   return(u)


def aplusbc(a,b,c):
   s=bcmul(b,c)
   t=bcadd(a,s)
   return(t)


def aminusbc(a,b,c):
   s=bcmul(b,c)
   t=bcsub(a,s)
   return(t)


# returns (a/b)/(c/d)
def ratior(a,b,c,d):
  global rationum
  global ratioden
  r=bcmul(a,d)
  s=bcmul(b,c)
  g=gcd(r,s)
  if ltzero(s):
     r=bcminus(r)
     s=bcminus(s)
  
  rationum=bcdiv(r,g)
  ratioden=bcdiv(s,g)
  return


# returns (a/b)(c/d)
def multr(a,b,c,d):
  global multnum
  global multden
  r=bcmul(a,c)
  s=bcmul(b,d)
  g=gcd(r,s)
  multnum=bcdiv(r,g)
  multden=bcdiv(s,g)
  return


def subr(a,b,c,d):
  global subnum
  global subden
  r=bcmul(a,d)
  s=bcmul(b,c)
  t=bcsub(r,s)
  u=bcmul(b,d)
  g=gcd(t,u)
  subnum=bcdiv(t,g)
  subden=bcdiv(u,g)
  return


def addr(a,b,c,d):
  global addnum
  global addden
  r=bcmul(a,d)
  s=bcmul(b,c)
  t=bcadd(r,s)
  u=bcmul(b,d)
  g=gcd(t,u)
  addnum=bcdiv(t,g)
  addden=bcdiv(u,g)
  return


# Assumes b>0 and d>0.  Returns -1, 0 or 1 according as a/b <,=,> c/d+ 
def comparer(a,b,c,d):
  t=abminuscd(a,d,b,c)
  if ltzero(t):
     return("-1")
  
  if gtzero(t):
     return("1")
  
  return("0")


def printlc(A,X,m):
 flag="0"
 s=zero[X,m]
 if ezero(s):
    return("0")
 
 for i in xrange(m):
     t=X[i]
     if ezero(t):
        continue
     
     if ezero(flag):
       if neq(t,"1") and neq(t,"-1"):
          print "t"print "b[i]"
       
       if eq(t,"-1"):
          print "-b[i]"
       
       if eq(t,"1"):
          print "b[i]"
       
       flag="1"
     else:
       if gtzero(t):
          if eq(t,"1"):
             print "+b[i]"
          else:
             print "+t"print "b[i]"
          
       
       if ltzero(t):
          if eq(t,"-1"):
             print "-b[i]"
          else:
             print "t"print "b[i]"
          
       
     
 
 return("1")


# lcv[j]=X[1]A[1][j]=...+X[m]A[m][j], 1 <= j <= n+ 
def lcasvector(A,X,m,n):
global lcv
   for j in xrange(n):
      sum="0"
      for i in xrange(m):
         t=bcmul(X[i],A[i][j])
         sum=bcadd(sum,t)
      
      lcv[j]=sum
   
   return

 def minusa(&a,n):
   for j in xrange(n):
       a[j]=bcminus(a[j])
   
   return


# This returns 1 if A is the zero matrix, otherwise returns 0+ 
 def test_zeromat(A,m,n):
      for i in xrange(m):
          for j in xrange(n):
              if neqzero(A[i][j]):
                 return("0")
              
          
      
      return("1")


def bcmul3(a,b,c):
 temp=bcmul(a,b)
 temp=bcmul(temp,c)
 return(temp)


def pparity(e):
   t=bcmod(e,"2")
   if eq(t,"1"):
     return("-1")
   else:
     return("1")
   


def printbinaryform(a,b,c,x,y):
	if gt(a,"1") || lt(a,"-1"):
           print"a&#8203x<sup>2</sup>"
        
        if eq(a,"1"):
           print"x<sup>2</sup>"
        
        if eq(a,"-1"):
           print"-x<sup>2</sup>"
        
	if neqzero(b):
	   if b > "1":
              print"+b&#8203xy"
           
	   if eq(b,"1"):
              print"+xy"
           
	   if eq(b,"-1"):
              print"-xy"
           
           if b < "-1":
              print"b&#8203xy"
           
        
	if neqzero(c):
	   if c > "1":
              print"+c&#8203y<sup>2</sup>"
           
	   if c < "-1":
              print"c&#8203y<sup>2</sup>"
           
           if eq(c,"1"):
              print"+y<sup>2</sup>"
           
           if eq(c,"-1"):
              print"-y<sup>2</sup>"
           
        


def bcmul4(a,b,c,d):
   temp1=bcmul(a,b)
   temp2=bcmul(c,d)
   temp=bcmul(temp1,temp2)
   return(temp)


def bcmul5(a,b,c,d,e):
   temp1=bcmul4(a,b,c,d)
   temp=bcmul(temp1,e)
   return(temp)

 
def bcadd4(a,b,c,d):
   temp1=bcadd3(a,b,c)
   temp=bcadd(temp1,d)
   return(temp)


def bcadd5(a,b,c,d,e):
   temp1=bcadd4(a,b,c,d)
   temp=bcadd(temp1,e)
   return(temp)


def bcadd6(a,b,c,d,e,f):
   temp1=bcadd5(a,b,c,d,e)
   temp=bcadd(temp1,f)
   return(temp)


def printaxplusby(a,x,b,y):
   if eq(a,"1"):
      print "x"
   else:
      print "-x"
   else:
      print "ax"
   
   if gtzero(b):
      if eq(b,"1"):
         print " + y"
      else:
         print " + by"
      
   
   if ltzero(b):
      if eq(b,"-1"):
         print " - y"
      else:
         minusb=bcminus(b)
         print " - minusby"
      
   
   return


def printxplusa(x,a):
   print x
   if gtzero(a):
      print " + a"
   
   if ltzero(a):
      minusa=bcminus(a)
      print " - minusa"
   
   return


def printaxplusbyplusc(a,x,b,y,c):
   null=printaxplusby(a,x,b,y)
      if gtzero(c):
         print " + c"
      
      if ltzero(c):
         minusc=bcminus(c)
         print " - minusc"
      
      return


def gcd4(a,b,c,d):

  t=gcd(a,b)
  t=gcd(t,c)
  t=gcd(t,d)
  return(t)


def printdXequalsaxplusbyplusc(d,X,a,x,b,y,c):
   t=gcd4(d,a,b,c)
   d=bcdiv(d,t)
   a=bcdiv(a,t)
   b=bcdiv(b,t)
   c=bcdiv(c,t)
   print "dX = "
   null=printaxplusbyplusc(a,x,b,y,c)
   return


