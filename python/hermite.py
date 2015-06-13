"""  lllhermite
  Input: integer mxn matrix A, nonzero, at least two rows+ 
  Output: small unimodular matrix B and HNF(A), such that BA=HNF(A)+ 
  The Havas, Majewski, Matthews LLL method is used+ 
  We usually take alpha=m1/n1, with (m1,n1)=(1,1) to get best results+ 
"""
#global col1
#global col2
#global n + 1
#global B
#global L
#global A
#global D
#global hnf
#global unimodular_matrix
#global rank


def  axb(Ab,m,n,m1,n1):
    """  A is m x n, b is m x 1, solving AX=b, X is n x 1+ 
      Ab is the (n+1) x m transposed augmented matrix. G=[A^t|0]
                                                     [b^t]1]
    """
    #global hnf
    #global unimodular_matrix
    #global rank
    for i in xrange(m + 1):
        for j in xrange(n):
            G[i][j]=Ab[i][j]
    for i in xrange(m):
        G[i][n + 1]=0
    G[m + 1][n + 1]=1
#    print "G="
#    printmat1(G,m + 1,n + 1)
#    print "<br>\n"
    lllhermite(G,m + 1,n + 1,m1,n1)
#    print "HNF(G)="
#    printmat1(hnf,m + 1,n + 1)
#    print "<br>\n"
#    print "P ="
#    printmat1(unimodular_matrix,m + 1,m + 1)
#    print "is a unimodular matrix such that PG = HNF(G)"
#    print "<br>\n"
    flag=0
    for i in xrange(rank - 1):
        if hnf[i][n + 1] != 0:
            flag=1
            break
    flag1=0
    for j in xrange(n):
        if hnf[rank][j] != 0:
            flag1=1
            break
    if flag == 0 and hnf[rank][n + 1] == 1 and flag1 == 0:
        #print "<img align=\"middle\" src=\"../jpgs/matrixP.png\"><br>\n"
        for j in xrange(m):
            y[j]=-unimodular_matrix[rank][j]
        #print "AX=B has a solution: Y = "
        #print[y,m]
        #print "<br>\n"
        nullity=m + 1 - rank
        if nullity == 0:
             print "AX=B has a unique solution in integers<br>\n"
             return
        else:
             lim=m + 1 - rank
             for i in xrange(lim):
                 rankplusi=rank + i
                 for j in xrange(m):
                     basis[i][j]=unimodular_matrix[rankplusi][j]
              if nullity == 1:
                  print "the row: "
              else:
                  print "the rows: "
              printmat1(basis,lim,m)
              if nullity == 1:
                  print "of submatrix R of P forms a Z-basis for the lattice AX=0<br>\n"
              else:
                  print "of submatrix R of P form a Z-basis for the lattice AX=0<br>\n"
    else:
         print "AX=B has no solution in integers<br>\n"
         return
    # joining basis and y
    for j in xrange(m):
         basis[lim + 1][j]=y[j]
    shortest_distance_axb(basis,lim + 1,m)
    return


def  lllhermite(G,m,n,m1,n1):
    """  G is a nonzero matrix with at least two rows.  """
    for i in xrange(m):
        for j in xrange(m):
            if i == j:
                 B[i][j]=1
            else:
                 B[i][j]=0
    for r in xrange(1, m):
        for s in xrange(r - 1):
            L[r][s]=0
    for i in xrange(m + 1):
        D[i]=1
    for i in xrange(m):
        for j in xrange(n):
            A[i][j]=G[i][j]
    flag=flagcol(A,m,n)
    if flag == 1:
        B[m][m]=-1
        for j in xrange(n):
            A[m][j]=-A[m][j]
    k=2
    while k <= m:
        reduce2(k,k - 1,m,n,D)
        minim=min(col2,n)
        temp1=D[k - 2] * D[k]
        temp2=L[k][k - 1] * L[k][k - 1]
        temp3=temp1 + temp2
        u=n1 * temp3
        temp1=D[k - 1] * D[k - 1]
        v=m1 * temp1
        if col1 <= minim || (col1 == col2 and col1 == n + 1 and u < v):
            swap2(k,m,n)
            if k > 2:
                k=k - 1
        else:
            for i in xrange(k - 2 - 1, 0, -1):
                reduce2(k,i,m,n,D)
            k=k + 1
    for i in xrange(m):
        for j in xrange(n):
            hnf[i][j]=A[i][j]
    for i in xrange(m):
        for j in xrange(m):
            unimodular_matrix[i][j]=B[i][j]
    for i in xrange(m - 1, 0, -1):
        test=zero_row_test(A,n,i)
        if test == 0:
            break
    rank=m - i
    for i in xrange(m):
        for j in xrange(n):
            k=m + 1 - i
            hnf[i][j]=A[k][j]
    for i in xrange(m):
        for j in xrange(m):
            k=m + 1 - i
            unimodular_matrix[i][j]=B[k][j]
    return


def  flagcol(A,m,n):
    """
      returns 0 if the first nonzero column j of A contains more than one
      nonzero entry, or contains only one nonzero entry and which is positive+ 
      returns 1 if the first nonzero column j of A contains only one nonzero entry,
      which is negative+ This assumes A is a nonzero matrix with at least two rows+ 
    """
    flag=0
    for j in xrange(n):
        for i in xrange(m):
            if(A[i][j] != 0)#found the first column with a nonzero elt, which is in row i
                flag=1
                break
        if flag == 1:
            break
    for k in xrange(i, m):
        if A[k][j]: != 0
            return 0
    if(A[i][j] > 0)# A[i][j] is the only elt in column j and is positive
        return 0
    else:# A[i][j] is the only elt in column j and is negative
        return 1

def reduce2(k,i,m,n,D):
    #global col1
    #global col2
    #global n + 1
    #global B
    #global L
    #global A
    col1=n + 1
    for j in xrange(n):
        if A[i][j] != 0:
            col1=j
         if A[i][col1] < 0:
             minus(i,m,L)
             for jj in xrange(n):
                 A[i][jj]=-A[i][jj]
             for jj in xrange(m):
                 B[i][jj]=-B[i][jj]
         break
    col2=n + 1
    for j in xrange(n):
         if A[k][j] != 0:
             col2=j
             break
    if col1 <= n:
         q=int(A[k][col1],A[i][col1])
    else:
         t=abs(L[k][i])
         t=2 * t
         if t > D[i]:
             q=lnearint(L[k][i],D[i])
         else:
             q=0
    if q != 0:
         for j in xrange(n):
             temp=q * A[i][j]
             A[k][j]=A[k][j] - temp
         for j in xrange(m):
             temp=q * B[i][j]
             B[k][j]=B[k][j] - temp
         temp=q * D[i]
         L[k][i]=L[k][i] - temp
         for j in xrange(i - 1):
             temp=q * L[i][j]
             L[k][j]=L[k][j] - temp


def minus(j,m,&L):
    for r in xrange(1, m):
        for s in xrange(r - 1):
            if r == j || s == j:
                L[r][s]=-L[r][s]
           
       
def swap2(k,m,n):
    #global B
    #global L
    #global A
    #global D
    #print "Row k <. Row k - 1<br>\n"
    for j in xrange(n):
         temp=A[k][j]
         A[k][j]=A[k - 1][j]
         A[k - 1][j]=temp
    for j in xrange(m):
         temp=B[k][j]
         B[k][j]=B[k - 1][j]
         B[k - 1][j]=temp
    for j in xrange(k - 2):
         temp=L[k][j]
         L[k][j]=L[k - 1][j]
         L[k - 1][j]=temp
   for i in xrange(k, m):
         temp1=L[i][k - 1] * D[k]
         temp2=L[i][k] * L[k][k - 1]
         t=temp1 - temp2
         temp1=L[i][k - 1] * L[k][k - 1]
         temp2=L[i][k] * D[k - 2]
         temp3=temp1 + temp2
         L[i][k - 1]=temp3 / D[k - 1]
         L[i][k]=t / D[k - 1]
   temp1=D[k - 2] * D[k]
   temp2=L[k][k - 1] * L[k][k - 1]
   t=temp1 + temp2
   D[k - 1]=t / D[k - 1]
   return


def  zero_row_test(matrix,n,i):
    """  This tests the i-th row of matrix to see if there is a nonzero
      entry. If there is one and the first occurs in column j, then j
      is returned. Otherwise 0 is returned+ 
    """
    for j in xrange(n):
        if matrix[i][j] != 0:
            return j
    return 0


def shortest_distance_axb(A,m,n):
    #global choleskynum
    #global choleskyden
    #global multnum
    #global multden
    #global addnum
    #global addden
    #global subnum
    #global subden
    #global rationum
    #global ratioden
    #global lcv
    count=0
    min_count=0
    #print "matrix A:"
    #printmat1(A,m,n)
    #print "<br>\n"
    #for j in xrange(n):
    # Am[j]=A[m][j]
    #
    #print "P = A[m] =  "print[A[m],n]print "<br>\n"
    # lengthj=dotproduct(Am,Am,m)
    m + 1=m
    #if m - 1 > 1:
       #print "&#8466 is the lattice spanned by the first m - 1 rows of A<br>\n"
    #else:
       #print "&#8466 is the lattice spanned by the first row of A<br>\n"
    #
    for i in xrange(m - 1):  # AA consists of the first m-1 rows of A
        for j in xrange(n):
            AA[i][j]=A[i][j]
    G=gram(A,m,n) 
    lengthj=G[m][m]
    cholesky(G,m)
    Qnum=choleskynum
    Qden=choleskyden
    QQnum=transpose(Qnum,m,m)
    QQden=transpose(Qden,m,m)
    m=m - 1
    for i in xrange(m):# the N vector
        Nnum[i]=Qnum[i][m + 1]
        Nden[i]=Qden[i][m + 1]
    Cnum=0
    Cden=1
    for i in xrange(m):
        multnum, multden = multr(Nnum[i],Nden[i],Nnum[i],Nden[i])
        multnum, multden = multr(multnum,multden,Qnum[i][i],Qden[i][i])
        addnum, addden = addr(Cnum,Cden,multnum,multden)
        Cnum=addnum
        Cden=addden
    i=m
    Tnum[m]=Cnum
    Tden[m]=Cden
    Unum[m]=0
    Uden[m]=1
    while 1:
        rationum, ratioden = ratior(Tnum[i],Tden[i],Qnum[i][i],Qden[i][i])
        Znum=rationum
        Zden=ratioden
        subnum, subden = subr(Nnum[i],Nden[i],Unum[i],Uden[i])
        UB[i]=introot(Znum,Zden,subnum,subden)
        subnum, subden = subr(Unum[i],Uden[i],Nnum[i],Nden[i])
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
                        temp=A[m + 1][k]
                        multiplier_vector[count][k]=temp - lcv[k]
                    l=lengthsquared(multiplier_vector[count],n)
                    multiplier_vector[count][n + 1]=l
                  # print "P-X[count]="print[multiplier_vector[count],n]print": l<br>\n"
                    lengtharray[count]=l
                    continue
                else:
                    i=i - 1
                    # now update U[i]
                    sumnum=0
                    sumden=1
                    for j in xrange(i, m):
                        multnum, multden = multr(Qnum[i][j],Qden[i][j],x[j],1)
                        addnum, addden = addr(sumnum,sumden,multnum,multden)
                        sumnum=addnum
                        sumden=addden
                    Unum[i]=sumnum
                    Uden[i]=sumden
                    # now update T[i]
                    addnum, addden = addr(x[i + 1],1,Unum[i + 1],Uden[i + 1])
                    subnum, subden = subr(addnum,addden,Nnum[i + 1],Nden[i + 1])
                    multnum, multden = multr(subnum,subden,subnum,subden)
                    multnum, multden = multr(Qnum[i + 1][i + 1],Qden[i + 1][i + 1],multnum,multden)
                    subnum, subden = subr(Tnum[i + 1],Tden[i + 1],multnum,multden)
                    Tnum[i]=subnum
                    Tden[i]=subden
                    break
            else:
                i=i + 1
                if i > m:
                    print "Here are the solution vectors with length squared &le lengthj<br>\n"
                    print "<TABLE BORDER=\1\ CELLSPACING=\0\>\n"
                    for k in xrange(count):
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
                    print "<TABLE BORDER=\0\ CELLSPACING=\0\>\n"
                    for k in xrange(count):
                        if multiplier_vector[k][n + 1] == min_length:
                            print "<TR>"
                            print "<TD ALIGN=\"RIGHT\">"
                            print[multiplier_vector[k],n]
                            print "</TD>"
                            print "</TR>\n"
                            min_count=min_count + 1
                    print "</TABLE>\n"
                    if min_count == 1:
                        print " is the shortest solution vector, length squared min_length<br>\n"
                    else:
                        print " are the shortest solution vectors, length squared min_length<br>\n"
                    return
                continue

def gram(A,m,n):
   for i in xrange(m):
       for j in xrange(m):
           B[i][j]=dotproduct(A[i],A[j],n)
   return B


def introot(a,b,c,d):
    """
    With Z=a/b, U=c/d, returns [sqrt(a/b)+c/d]. First ANSWER = [sqrt(Z)] + [U]. One then 
    tests if Z < ([sqrt(Z)] + 1 -U)^2. If this does not hold, ANSWER += 1+ 
    For use in fincke_pohst()+ 
    """
    y=int(c,d)
    if a == 0:
        return y
    x=a / b
    x=sqrt(x)
    answer=x + y
    subnum, subden = subr(c,d,y,1)
    subnum, subden = subr(1,1,subnum,subden)
    addnum, addden = addr(x,1,subnum,subden)
    multnum, multden = multr(addnum,addden,addnum,addden)
    t=comparer(multnum,multden,a,b)
    if le(t,0):
        answer=answer + 1
    return answer

def lengthsquared(a,n):
    sm=0
    for i in xrange(n):
        temp=a[i] * a[i]
        sm=sm + temp
    return sm


def int(a,b):
    if b<0:
         a=0 - a
         b=0 - b
    c=a / b
    d=a % b
    if d==0 || a>0:
        return c
    else:
        return c - 1
    

def sign(a):
    """  sign of an integer a  """
    """  sign(a)=1,-1,0, according as a>0,a<0,a=0  """
    if a>0:
        return 1
    if a<0:
        return -1
    return 0


def gcd(m,n):
    """   b=gcd(m,n) for any integers m and n  """
    """  Euclid's division algorithm is used.  """
    """  We use gcd(m,n)=gcd(|m|,|n|)  """
    a=abs(m)  # a=r[0]
    if n==0:
         return a
    b=abs(n)  # b=r[1]
    c=mod(a,b)  # c=r[2]=r[0] mod(r[1])
    while c:
        a=b
        b=c
        c=mod(a,b)  # c=r[j]=r[j-2] mod(r[j-1]) 
    return b
   

def egcd(p,q):
    if q==0:
        if p!=0:
            s=sign(p)
            if s==1:
                multiplier1=1
            else:
                multiplier1=-1
            return abs(p), multiplier1, 0
        else:
            return 0, 0, 0
    a=p
    b=abs(q)
    c=a % b
    s=sign(q)
    if c==0:
        if s==1:
            multiplier2=1
        else:
            multiplier2=-1
        return b, 0, multiplier2
    l1=1
    k1=0
    l2=0
    k2=1
    while c!=0:
        q=int(a,b)
        a=b
        b=c
        c=a % b
        temp1=q * k1
        temp2=q * k2
        h1=l1 - temp1
        h2=l2 - temp2
        l1=k1
        l2=k2
        k1=h1
        k2=h2
    if s==-1:
        k2=0 - k2
    return b, k1, k2


def num_digits(n):
    """  If n > 0, len(n) returns the number of base 10 digits of n  """
    i=0
    x=abs(n)
    while x!=0:
        x=int(x,10)
        i=i + 1
    return i


def  lnearint(a,b):
    """
    left nearest integer 
    returns y+1/2 if a/b=y+1/2, y integral+ 
    """
    y=int(a,b)
    if b < 0:
        a=-a
        b=-b
    x=b * y
    z=a - x
    z=2 * z
    if z > b:
        y=y + 1
    return y


def mina(a,n):
    x=a[1]
    for i in xrange(1, n):
    if a[i] < x:
        x=a[i]
    return x


def dotproduct(a,b,n):
    sum=0
    for j in xrange(n):
        temp=a[j] * b[j]
        sum=sum + temp
    return sum


def abminuscd(a,b,c,d):
    s=a * b
    t=c * d
    u=s - t
    return u


def rationum, ratioden = ratior(a,b,c,d):
    """ returns (a/b)/(c/d)"""
    r=a * d
    s=b * c
    g=gcd(r,s)
    if s < 0:
        r=-r
        s=-s
    return r / g, s / g


def multnum, multden = multr(a,b,c,d):
    # returns (a/b)(c/d)
    r=a * c
    s=b * d
    g=gcd(r,s)
    return r / g, s / g


def subnum, subden = subr(a,b,c,d):
    r=a * d
    s=b * c
    t=r - s
    u=b * d
    g=gcd(t,u)
    return t / g, u / g


def addnum, addden = addr(a,b,c,d):
    r=a * d
    s=b * c
    t=r + s
    u=b * d
    g=gcd(t,u)
    return t / g, u / g


def comparer(a,b,c,d):
    """ Assumes b>0 and d>0.  Returns -1, 0 or 1 according as a/b <,=,> c/d+ """
    t=abminuscd(a,d,b,c)
    if t < 0:
        return -1
    if t > 0:
        return 1
    return 0


def lcasvector(A,X,m,n):
    """lcv[j]=X[1]A[1][j]=...+X[m]A[m][j], 1 <= j <= n+"""
    #global lcv
    for j in xrange(n):
        sum=0
        for i in xrange(m):
            t=X[i] * A[i][j]
            sum=sum + t
        lcv[j]=sum
    return lcv


def minusa(&a,n):
    for j in xrange(n):
        a[j]=-a[j]
    return


def test_zeromat(A,m,n):
    """ This returns 1 if A is the zero matrix, otherwise returns 0+ """
    for i in xrange(m):
        for j in xrange(n):
            if A[i][j] != 0:
                return 0
    return 1

