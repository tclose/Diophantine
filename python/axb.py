include("check_input.php")
include("library.php")
include("lllhermite_.php")
global transposed

#print "<pre>" print str(_POST) print "</pre>"

#rows=_POST['rows']
#cols=_POST['cols']
rows=trim(rows)
cols=trim(cols)
check1=check_decimal(rows)
check2=check_decimal(cols)
if check1=="0" || check2=="0":
        print "<p>\n"
        exit

if rows > "50":
  print "row dimension is &gt 50<br>\n"
  print "<a href=\"./axb.html\">Return to main page</a><br>\n"
  exit

if cols > "50":
  print "column dimension is &gt 50<br>\n"
  print "<a href=\"./axb.html\">Return to main page</a><br>\n"
  exit

if rows < "1":
  print "row dimension is &lt 1<br>\n"
  print "<a href=\"./axb.html\">Return to main page</a><br>\n"
  exit

if cols < "1":
  print "column dimension is &lt 1<br>\n"
  print "<a href=\"./axb.html\">Return to main page</a><br>\n"
  exit

matrix=trim(matrix)
a=split('[ ]+',matrix)
t=len(a)
if(le(t,"1"))
	print "number of entries is less than or equal to 1<br>\n"
        print "<a href=\"./axb.html\">Return to main page</a><br>\n"
	flush()
        exit

cols=bcadd(cols,"1")"""  Here cols is the number of columns of the augmented matrix  """
size=bcmul(rows,cols)
if t < size:
	print "number t of entries is less than m &times n = size<br>\n"
        print "<a href=\"./axb.html\">Return to main page</a><br>\n"
	flush()
        exit

if t > size:
	print "number t of entries is greater than m &times n = size<br>\n"
        print "<a href=\"./axb.html\">Return to main page</a><br>\n"
	flush()
        exit


	flag="0"
	for(i="0"lt(i,t)i=bcadd(i,"1"))
		check=check_decimal(a[i])
		if ezero(check:)
			print "<p>\n"
			flush()
			flag="1"
			break
		
	
	if(ezero(flag))
                ii="0"
		for(i="1"le(i,rows)i=bcadd(i,"1"))
		    for(j="1"le(j,cols)j=bcadd(j,"1"))
                        k=bcadd(ii,j)
                        k=bcsub(k,"1")
                        mat[i][j]=a[k]
                    
                    ii=bcadd(cols,ii)
		
                m=bcsub(cols,"1")
                t=test_zeromat(mat,rows,m)
                if(eq(t,"1"))
                   print "Coeffficient matrix is the zero matrix<br>\n"
                else:
                   print "Augmented matrix [A|B]="
                   printmat1(mat,rows,cols)
                   print "<br>\n"
                   transposed=transpose(mat,rows,cols)
                   m1="1"
                   n1="1"
                   axb(transposed,m,rows,m1,n1)
                   print "<br>\n"
                
	
print "<p>\n"
flush()
print "<a href=\"./axb.html\">Return to main page</a><br>\n"
flush()

</body>
</html>
