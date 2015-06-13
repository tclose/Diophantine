
def cholesky(A,m): # A is positive definite mxm
  global rationum
  global ratioden
  global multnum
  global multden
  global subnum
  global subden
  global choleskynum
  global choleskyden

  Qnum=A
  for i in xrange(m):
      for j in xrange(m):
          Qden[i][j]=1
      
  
  for i in xrange(1, m):
     iplus1=i + 1
     for j in xrange(i, m):
         Qnum[j][i]=Qnum[i][j]
         Qden[j][i]=Qden[i][j]
         ratior(Qnum[i][j],Qden[i][j],Qnum[i][i],Qden[i][i])
         Qnum[i][j]=rationum
         Qden[i][j]=ratioden
     
     for k in xrange(i, m):
         for l in xrange(k - 1, m):
               multr(Qnum[k][i],Qden[k][i],Qnum[i][l],Qden[i][l])
               t2num=multnum
               t2den=multden
               subr(Qnum[k][l],Qden[k][l],t2num,t2den)
               Qnum[k][l]=subnum
               Qden[k][l]=subden
         
     
  
  for i in xrange(m):
      for j in xrange(m):
          choleskynum[i][j]=Qnum[i][j]
          choleskyden[i][j]=Qden[i][j]
      
  
  return


def gram(A,m,n):
   for i in xrange(m):
       for j in xrange(m):
           B[i][j]=dotproduct(A[i],A[j],n)
       
   
   return(B)


# With Z=a/b, U=c/d, returns [sqrt(a/b)+c/d]. First ANSWER = [sqrt(Z)] + [U]. One then 
# tests if Z < ([sqrt(Z)] + 1 -U)^2. If this does not hold, ANSWER += 1+ 
# For use in fincke_pohst()+ 
def introot(a,b,c,d):
global subnum
global subden
global addnum
global addden
global multnum
global multden
   y=int(c,d)
   if ezero(a):
      return(y)
   
   x=a / b
   x=bcsqrt(x)
   answer=x + y
   subr(c,d,y,1)
   subr(1,1,subnum,subden)
   addr(x,1,subnum,subden)
   multr(addnum,addden,addnum,addden)
   t=comparer(multnum,multden,a,b)
   if le(t,0):
      answer=answer + 1
   
   return(answer)


def shortest_distance(A,m,n):
global choleskynum
global choleskyden
global multnum
global multden
global addnum
global addden
global subnum
global subden
global rationum
global ratioden
global lcv

    count=0
    min_count=0

    print "matrix A:"
    printmat1(A,m,n)
    print "<br>\n"
    nplus1=n + 1
    #for j in xrange(n):
    # Am[j]=A[m][j]
  #
    print "P = A[m] =  "print[A[m],n]print "<br>\n"
   # lengthj=dotproduct(Am,Am,m)

    mplus1=m
    mminus1=m - 1
    if mminus1 > 1:
       print "&#8466 is the lattice spanned by the first mminus1 rows of A<br>\n"
    else:
       print "&#8466 is the lattice spanned by the first row of A<br>\n"
    
    for i in xrange(mminus1):  # AA consists of the first m-1 rows of A
        for j in xrange(n):
            AA[i][j]=A[i][j]
        
    
    G=gram(A,m,n) 
    lengthj=G[m][m]
    cholesky(G,m)
    Qnum=choleskynum
    Qden=choleskyden
    QQnum=transpose(snum,m,m)
    QQden=transpose(sden,m,m)
    m=m - 1
    for i in xrange(m):# the N vector
        Nnum[i]=Qnum[i][mplus1]
        Nden[i]=Qden[i][mplus1]
    

    Cnum=0
    Cden=1
    for i in xrange(m):
        multr(Nnum[i],Nden[i],Nnum[i],Nden[i])
        multr(multnum,multden,Qnum[i][i],Qden[i][i])
        addr(Cnum,Cden,multnum,multden)
        Cnum=addnum
        Cden=addden
    
    i=m
    Tnum[m]=Cnum
    Tden[m]=Cden
    Unum[m]=0
    Uden[m]=1
    while 1:
       ratior(Tnum[i],Tden[i],Qnum[i][i],Qden[i][i])
       Znum=rationum
       Zden=ratioden
       subr(Nnum[i],Nden[i],Unum[i],Uden[i])
       UB[i]=introot(Znum,Zden,subnum,subden)
       subr(Unum[i],Uden[i],Nnum[i],Nden[i])
       temp2=introot(Znum,Zden,subnum,subden)
       temp3=-temp2
       x[i]=temp3 - 1
       while 1:
          x[i]=x[i] + 1
          if le(x[i],UB[i]):
              if i == 1:
                   #s=printlc(A,x,m)
                   lcasvector(AA,x,m,n)
                   count=count + 1
               #  print "X[count]="print[x,m]
                   lcva[count]=lcv
               #  print "lcv[count]="print[lcv,n]
               #  print "<br>\n"
                   coord[count]=x
                   for k in xrange(n):
                       temp=A[mplus1][k]
                       multiplier_vector[count][k]=temp - lcv[k]
                   
                   l=lengthsquared(multiplier_vector[count],n)
                   multiplier_vector[count][nplus1]=l
                # print "P-X[count]="print[multiplier_vector[count],n]print": l<br>\n"
                   lengtharray[count]=l
                   continue
              else:
                i=i - 1
                # now update U[i]
                sumnum=0
                sumden=1
                iplus1=i + 1
                for j in xrange(i, m):
                    multr(Qnum[i][j],Qden[i][j],x[j],1)
                    addr(sumnum,sumden,multnum,multden)
                    sumnum=addnum
                    sumden=addden
                
                Unum[i]=sumnum
                Uden[i]=sumden
                # now update T[i]
                addr(x[iplus1],1,Unum[iplus1],Uden[iplus1])
                subr(addnum,addden,Nnum[iplus1],Nden[iplus1])
                multr(subnum,subden,subnum,subden)
                multr(Qnum[iplus1][iplus1],Qden[iplus1][iplus1],multnum,multden)
                subr(Tnum[iplus1],Tden[iplus1],multnum,multden)
                Tnum[i]=subnum
                Tden[i]=subden
                break
              
          else:
             i=i + 1
             if i > m:
                    print "Here are the X[k] &isin &#8466, P - X[k], ||P-X[k]||<sup>2</sup> such that ||P-X[k]||<sup>2</sup> &le lengthj<br>\n"
                    print "<TABLE BORDER=\1\ CELLSPACING=\0\>\n"
                    for k in xrange(count):
                           print "<TR>"
                  #       print "<TD ALIGN=\"RIGHT\">"
                  #       minusa(coord[k],m) 
                  #       s=printlc(AA,coord[k],m)
                  #       if ezero(s):
                  #          print "b[mplus1]"
                  #       else:
                  #          print "+b[mplus1]"
                  #       
                  #       print "="
                  #       print "</TD>"
                           print "<TD ALIGN=\"RIGHT\">"
                           print[lcva[k],n]
                           print "</TD>"
                           print "<TD ALIGN=\"RIGHT\">"
                           print[multiplier_vector[k],n]
                           print "</TD>"
                           print "<TD>"
                           print " lengtharray[k]"
                           print "</TD>"
                           print "</TR>\n"
                    
                    print "</TABLE>\n"
                    print "Also<br>\n"
                    min_length=mina(lengtharray,count)
                    print "<TABLE BORDER=\0\ CELLSPACING=\0\>\n"
                    for k in xrange(count):
                        if multiplier_vector[k][nplus1] == min_length:
                           print "<TR>"
                           print "<TD ALIGN=\"RIGHT\">"
                           print[lcva[k],n]
                           print "</TD>"
                           print "</TR>\n"
                           min_count=min_count + 1
                        
                    
                    print "</TABLE>\n"
                    if min_count == 1:
                       print " is the closest vector of &#8466 to P, with shortest distance squared min_length<br>\n"
                    else:
                       print " are the closest vectors of &#8466 to P, with shortest distance squared min_length<br>\n"
                    
                    return
             
             continue
          
       
   


def lengthsquared(a,n):
   sum=0
   for i in xrange(n):
      temp=a[i] * a[i]
      sum=sum + temp
   
   return(sum)


# returns 0 if all of a[i] are zero, otherwise 1+ 
def zero[a,n]:
  flag=0
  for i in xrange(n):
      if a[i] != 0:
        flag=1
        break
      
  
  return(flag)



