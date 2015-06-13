
"""  php program lllhermite_.php  """
"""  Input: integer mxn matrix A, nonzero, at least two rows+ 
  Output: small unimodular matrix B and HNF(A), such that BA=HNF(A)+ 
  The Havas, Majewski, Matthews LLL method is used+ 
  We usually take alpha=m1/n1, with (m1,n1)=(1,1) to get best results+ 
  """

include("closest_lattice_points_.php")

global col1
global col2
global nplus1
global B
global L
global A
global D
global hnf
global unimodular_matrix
global rank

def  lllhermite(G,m,n,m1,n1):
   """  G is a nonzero matrix with at least two rows.  """
   for(i="1"le(i,m)i=bcadd(i,"1"))
       for(j="1"le(j,m)j=bcadd(j,"1"))
           if(eq(i,j))
              B[i][j]="1"
           else:
              B[i][j]="0"
           
        
   
   for(r="2"le(r,m)r=bcadd(r,"1"))
       for(s="1"lt(s,r)s=bcadd(s,"1"))
           L[r][s]="0"
       
   
   for(i="0"le(i,m)i=bcadd(i,"1"))
        D[i]="1"
   
   for(i="1"le(i,m)i=bcadd(i,"1"))
       for(j="1"le(j,n)j=bcadd(j,"1"))
           A[i][j]=G[i][j]
       
   

   flag=flagcol(A,m,n)
   if(eq(flag,"1"))
      B[m][m]="-1"
      for(j="1"le(j,n)j=bcadd(j,"1"))
          A[m][j]=bcminus(A[m][j])
      
   
   k="2"
   nplus1=bcadd(n,"1")
   while(le(k,m))
         kminus1=bcsub(k,"1")
         reduce2(k,kminus1,m,n,D)
         kminus2=bcsub(k,"2")
         minim=minimum(col2,n)
         temp1=bcmul(D[kminus2],D[k])
         temp2=bcmul(L[k][kminus1],L[k][kminus1])
         temp3=bcadd(temp1,temp2)
         u=bcmul(n1,temp3)
         temp1=bcmul(D[kminus1],D[kminus1])
         v=bcmul(m1,temp1)
         if(le(col1,minim) || (eq(col1,col2) and eq(col1,nplus1) and lt(u,v)))
            swap2(k,m,n)
            if k > "2":
               k=kminus1
            
         else:
            for(i=kminus2ge(i,"1")i=bcsub(i,"1"))
                reduce2(k,i,m,n,D)
            
            k=bcadd(k,"1")
         
   
   for(i="1"le(i,m)i=bcadd(i,"1"))
       for(j="1"le(j,n)j=bcadd(j,"1"))
           hnf[i][j]=A[i][j]
       
   
   for(i="1"le(i,m)i=bcadd(i,"1"))
       for(j="1"le(j,m)j=bcadd(j,"1"))
           unimodular_matrix[i][j]=B[i][j]
       
   
   for(i=mge(i,"1")i=bcsub(i,"1"))
       test=zero_row_test(A,n,i)
       if(ezero(test))
          break
       
   
   rank=bcsub(m,i)
   #print "rank = rank<br>\n"
   mplus1=bcadd(m,"1")
   for(i="1"le(i,m)i=bcadd(i,"1"))
        for(j="1"le(j,n)j=bcadd(j,"1"))
           k=bcsub(mplus1,i)
           hnf[i][j]=A[k][j]
       
   
   for(i="1"le(i,m)i=bcadd(i,"1"))
       for(j="1"le(j,m)j=bcadd(j,"1"))
           k=bcsub(mplus1,i)
           unimodular_matrix[i][j]=B[k][j]
       
   
#   rankplus1=bcadd(rank,"1")
#   for(i=rankplus1le(i,m)i=bcadd(i,"1"))
 #      for(j="1"le(j,n)j=bcadd(j,"1"))
  #         hnf[i][j]="0"
   #    
 #  
 #  for(i=rankplus1le(i,m)i=bcadd(i,"1"))
  #     for(j="1"le(j,n)j=bcadd(j,"1"))
   #        hnf[i][j]="0"
    #   
 #  
   return


def  flagcol(A,m,n):
    """  returns 0 if the first nonzero column j of A contains more than one
      nonzero entry, or contains only one nonzero entry and which is positive+ 
      returns 1 if the first nonzero column j of A contains only one nonzero entry, which is negative+ 
      This assumes A is a nonzero matrix with at least two rows+ 
    """
    flag="0"
    for(j="1"le(j,n)j=bcadd(j,"1"))
       for(i="1"le(i,m)i=bcadd(i,"1"))
           if(neqzero(A[i][j]))#found the first column with a nonzero elt, which is in row i
              flag="1"
              break
           
       if(eq(flag,"1"))
          break
       
    iplus1=bcadd(i,"1")
    for(k=iplus1le(k,m)k=bcadd(k,"1"))
       if neqzero(A[k][j]:)
           return("0")

    if(gtzero(A[i][j]))# A[i][j] is the only elt in column j and is positive
      return("0")
    else:# A[i][j] is the only elt in column j and is negative
      return("1")

    #   found:
    #  if i < m:
    #    return("0")
    #else:
    # if(ltzero(A[m][j]))
     #   return("1")
    # else:
     #   return("0")
     #
     #

def reduce2(k,i,m,n,D):
global col1
global col2
global nplus1
global B
global L
global A
   col1=nplus1
   for(j="1"le(j,n)j=bcadd(j,"1"))
       if(neqzero(A[i][j]))
         col1=j
         if(ltzero(A[i][col1]))
            minus(i,m,L)
#print "Row i . -Row i<br>\n"
            for(jj="1"le(jj,n)jj=bcadd(jj,"1"))
                A[i][jj]=bcminus(A[i][jj])
            
            for(jj="1"le(jj,m)jj=bcadd(jj,"1"))
                B[i][jj]=bcminus(B[i][jj])
            
         
         break
       
   
   col2=nplus1
   for(j="1"le(j,n)j=bcadd(j,"1"))
       if(neqzero(A[k][j]))
         col2=j
         break
       
   
   if(le(col1,n))
      q=int(A[k][col1],A[i][col1])
   else:
      t=bcabs(L[k][i])
      t=bcmul("2",t)
      if t > D[i]:
        q=lnearint(L[k][i],D[i])
      else:
        q="0"
      
   
   if(neqzero(q))
#print "Row k . Row k - q &times Row i<br>\n"
      for(j="1"le(j,n)j=bcadd(j,"1"))
          temp=bcmul(q,A[i][j])
          A[k][j]=bcsub(A[k][j],temp)
      
      for(j="1"le(j,m)j=bcadd(j,"1"))
          temp=bcmul(q,B[i][j])
          B[k][j]=bcsub(B[k][j],temp)
      
      temp=bcmul(q,D[i])
      L[k][i]=bcsub(L[k][i],temp)
      for(j="1"lt(j,i)j=bcadd(j,"1"))
          temp=bcmul(q,L[i][j])
          L[k][j]=bcsub(L[k][j],temp)
      
   


def minus(j,m,&L):
   for(r="2"le(r,m)r=bcadd(r,"1"))
       for(s="1"lt(s,r)s=bcadd(s,"1"))
           if(eq(r,j) || eq(s,j))
             L[r][s]=bcminus(L[r][s])
           
       
   


def swap2(k,m,n):
global B
global L
global A
global D
   
   kminus1=bcsub(k,"1")
   #print "Row k <. Row kminus1<br>\n"
   for(j="1"le(j,n)j=bcadd(j,"1"))
       temp=A[k][j]
       A[k][j]=A[kminus1][j]
       A[kminus1][j]=temp
   
   for(j="1"le(j,m)j=bcadd(j,"1"))
       temp=B[k][j]
       B[k][j]=B[kminus1][j]
       B[kminus1][j]=temp
   
   kminus2=bcsub(k,"2")
   for(j="1"le(j,kminus2)j=bcadd(j,"1"))
       temp=L[k][j]
       L[k][j]=L[kminus1][j]
       L[kminus1][j]=temp
   
   kplus1=bcadd(k,"1")
   for(i=kplus1le(i,m)i=bcadd(i,"1"))
       temp1=bcmul(L[i][kminus1],D[k])
       temp2=bcmul(L[i][k],L[k][kminus1])
       t=bcsub(temp1,temp2)
       temp1=bcmul(L[i][kminus1],L[k][kminus1])
       temp2=bcmul(L[i][k],D[kminus2])
       temp3=bcadd(temp1,temp2)
       L[i][kminus1]=bcdiv(temp3,D[kminus1])
       L[i][k]=bcdiv(t,D[kminus1])
   

   temp1=bcmul(D[kminus2],D[k])
   temp2=bcmul(L[k][kminus1],L[k][kminus1])
   t=bcadd(temp1,temp2)
#var_dump(t)
   D[kminus1]=bcdiv(t,D[kminus1])
   return


def  zero_row_test(matrix,n,i):
"""  This tests the i-th row of matrix to see if there is a nonzero
  entry. If there is one and the first occurs in column j, then j
  is returned. Otherwise 0 is returned+ 
  """
    for(j="1"le(j,n)j=bcadd(j,"1"))
       if(neqzero(matrix[i][j]))
         return(j)
       
    
    return("0")


def  axb(Ab,m,n,m1,n1):
"""  A is m x n, b is m x 1, solving AX=b, X is n x 1+ 
  Ab is the (n+1) x m transposed augmented matrix. G=[A^t|0]
                                                     [b^t]1]
  """
global hnf
global unimodular_matrix
global rank
    mplus1=bcadd(m,"1")
    for(i="1"le(i,mplus1)i=bcadd(i,"1"))
        for(j="1"le(j,n)j=bcadd(j,"1"))
               G[i][j]=Ab[i][j]
        
    
    nplus1=bcadd(n,"1")
    for(i="1"le(i,m)i=bcadd(i,"1"))
        G[i][nplus1]="0"
    
    G[mplus1][nplus1]="1"
    print "G="
    printmat1(G,mplus1,nplus1)
    print "<br>\n"
    lllhermite(G,mplus1,nplus1,m1,n1)
    print "HNF(G)="
    printmat1(hnf,mplus1,nplus1)
    print "<br>\n"
    print "P ="
    printmat1(unimodular_matrix,mplus1,mplus1)
    print "is a unimodular matrix such that PG = HNF(G)"
    print "<br>\n"
    flag="0"
    for(i="1"lt(i,rank)i=bcadd(i,"1"))
        if(neqzero(hnf[i][nplus1]))
           flag="1"
           break
        
    
    flag1="0"
    for(j="1"le(j,n)j=bcadd(j,"1"))
        if(neqzero(hnf[rank][j]))
           flag1="1"
           break
        
    
    #t=hnf[rank][nplus1] this was erroneous - fixed 25th October 2011 thanks to an example of Mostafa
    #Khorramizadeh, Int, J, Computing math. 86, issue 5,2009, 883-896
    if(eq(flag,"0") and eq(hnf[rank][nplus1],"1") and eq(flag1,"0"))
        print "<img align=\"middle\" src=\"../jpgs/matrixP.png\"><br>\n"
        for(j="1"le(j,m)j=bcadd(j,"1"))
          y[j]=bcminus(unimodular_matrix[rank][j])
        
        print "AX=B has a solution: Y = "
        print[y,m]
        print "<br>\n"
        nullity=bcsub(mplus1,rank)
        if(ezero(nullity))
           print "AX=B has a unique solution in integers<br>\n"
           return
        else:
           lim=bcsub(mplus1,rank)
           for(i="1"le(i,lim)i=bcadd(i,"1"))
               rankplusi=bcadd(rank,i)
               for(j="1"le(j,m)j=bcadd(j,"1"))
                   basis[i][j]=unimodular_matrix[rankplusi][j]
               
           
           if(eq(nullity,"1"))
              print "the row: "
           else:
              print "the rows: "
           
           printmat1(basis,lim,m)
           if(eq(nullity,"1"))
              print "of submatrix R of P forms a Z-basis for the lattice AX=0<br>\n"
           else:
              print "of submatrix R of P form a Z-basis for the lattice AX=0<br>\n"
           
        
    else:
       print "AX=B has no solution in integers<br>\n"
       return
    
    # joining basis and y
    limplus1=bcadd(lim,"1")
    for(j="1"le(j,m)j=bcadd(j,"1"))
        basis[limplus1][j]=y[j]
    
    shortest_distance_axb(basis,limplus1,m)
    return


def shortest_distance_axb(A,m,n):
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

    count="0"
    min_count="0"

    #print "matrix A:"
    #printmat1(A,m,n)
    #print "<br>\n"
    nplus1=bcadd(n,"1")
    #for(j="1"le(j,n)j=bcadd(j,"1"))
    # Am[j]=A[m][j]
  #
    #print "P = A[m] =  "print[A[m],n]print "<br>\n"
   # lengthj=dotproduct(Am,Am,m)

    mplus1=m
    mminus1=bcsub(m,"1")
    #if mminus1 > "1":
       #print "&#8466 is the lattice spanned by the first mminus1 rows of A<br>\n"
    #else:
       #print "&#8466 is the lattice spanned by the first row of A<br>\n"
    #
    for(i="1"le(i,mminus1)i=bcadd(i,"1"))  # AA consists of the first m-1 rows of A
        for(j="1"le(j,n)j=bcadd(j,"1"))
            AA[i][j]=A[i][j]
        
    
    G=gram(A,m,n) 
    lengthj=G[m][m]
    cholesky(G,m)
    Qnum=choleskynum
    Qden=choleskyden
    QQnum=transpose(Qnum,m,m)
    QQden=transpose(Qden,m,m)
    m=bcsub(m,"1")
    for(i="1"le(i,m)i=bcadd(i,"1"))# the N vector
        Nnum[i]=Qnum[i][mplus1]
        Nden[i]=Qden[i][mplus1]
    

    Cnum="0"
    Cden="1"
    for(i="1"le(i,m)i=bcadd(i,"1"))
        multr(Nnum[i],Nden[i],Nnum[i],Nden[i])
        multr(multnum,multden,Qnum[i][i],Qden[i][i])
        addr(Cnum,Cden,multnum,multden)
        Cnum=addnum
        Cden=addden
    
    i=m
    Tnum[m]=Cnum
    Tden[m]=Cden
    Unum[m]="0"
    Uden[m]="1"
    while(1)
       ratior(Tnum[i],Tden[i],Qnum[i][i],Qden[i][i])
       Znum=rationum
       Zden=ratioden
       subr(Nnum[i],Nden[i],Unum[i],Uden[i])
       UB[i]=introot(Znum,Zden,subnum,subden)
       subr(Unum[i],Uden[i],Nnum[i],Nden[i])
       temp2=introot(Znum,Zden,subnum,subden)
       temp3=bcminus(temp2)
       x[i]=bcsub(temp3,"1")
       while("1")
          x[i]=bcadd(x[i],"1")
          if(le(x[i],UB[i]))
              if(eq(i,"1"))
                   #s=printlc(A,x,m)
                   lcasvector(AA,x,m,n)
                   count=bcadd(count,"1")
               #  print "X[count]="print[x,m]
                   lcva[count]=lcv
               #  print "lcv[count]="print[lcv,n]
               #  print "<br>\n"
                   coord[count]=x
                   for(k="1"le(k,n)k=bcadd(k,"1"))
                       temp=A[mplus1][k]
                       multiplier_vector[count][k]=bcsub(temp,lcv[k])
                   
                   l=lengthsquared(multiplier_vector[count],n)
                   multiplier_vector[count][nplus1]=l
                # print "P-X[count]="print[multiplier_vector[count],n]print": l<br>\n"
                   lengtharray[count]=l
                   continue
              else:
                i=bcsub(i,"1")
                # now update U[i]
                sumnum="0"
                sumden="1"
                iplus1=bcadd(i,"1")
                for(j=iplus1le(j,m)j=bcadd(j,"1"))
                    multr(Qnum[i][j],Qden[i][j],x[j],"1")
                    addr(sumnum,sumden,multnum,multden)
                    sumnum=addnum
                    sumden=addden
                
                Unum[i]=sumnum
                Uden[i]=sumden
                # now update T[i]
                addr(x[iplus1],"1",Unum[iplus1],Uden[iplus1])
                subr(addnum,addden,Nnum[iplus1],Nden[iplus1])
                multr(subnum,subden,subnum,subden)
                multr(Qnum[iplus1][iplus1],Qden[iplus1][iplus1],multnum,multden)
                subr(Tnum[iplus1],Tden[iplus1],multnum,multden)
                Tnum[i]=subnum
                Tden[i]=subden
                break
              
          else:
             i=bcadd(i,"1")
             if i > m:
                    print "Here are the solution vectors with length squared &le lengthj<br>\n"
                    print "<TABLE BORDER=\"1\" CELLSPACING=\"0\">\n"
                    for(k="1"le(k,count)k=bcadd(k,"1"))
                           print "<TR>"
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
                    print "<TABLE BORDER=\"0\" CELLSPACING=\"0\">\n"
                    for(k="1"le(k,count)k=bcadd(k,"1"))
                        if(eq(multiplier_vector[k][nplus1],min_length))
                           print "<TR>"
                           print "<TD ALIGN=\"RIGHT\">"
                           print[multiplier_vector[k],n]
                           print "</TD>"
                           print "</TR>\n"
                           min_count=bcadd(min_count,"1")
                        
                    
                    print "</TABLE>\n"
                    if(eq(min_count,"1"))
                       print " is the shortest solution vector, length squared min_length<br>\n"
                    else:
                       print " are the shortest solution vectors, length squared min_length<br>\n"
                    
                    return
             
             continue

