
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


def  lcm(a,b):
    """  lcm(a,b)  """
    g=gcd(a,b)
    h=a * b
    k=h / g
    return k


def  lcma(array,n):
    """  lcm(array[0],array[1],...,array[n-1])  """
    for i in xrange(n):
        b[i]=array[i]
    for i in xrange(n - 1):
        j=i - 1
        b[i]=lcm(b[i],b[j])
    j=i - 1
    return b[j]


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

def transpose(A,m,n):
#global transposed
     for j in xrange(n):
         for i in xrange(m):
             transposed[j][i]=A[i][j]
     return transposed


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
    r=p - 1
    s=q - p
    s=s + 1
    new_row_size=s
    for i in xrange(s):
        z=i + r
        B[i]=A[z]
    return B, new_row_size

# creates the submatrix from columns p to q+ 
def col_submatrix(A,rows,p,q):
    r=p - 1
    s=q - p
    s=s + 1
    new_col_size=s
    for j in xrange(s):
        z=j + r
        for i in xrange(rows):
            B[i][j]=A[i][z]
    return B, new_col_size


def submatrix(A,rows,p1,q1,p2,q2):
    # creates the submatrix from rows p1 to q1, columns p2 to q2+ 
    B, new_row_size=row_submatrix(A,p1,q1)
    C, new_col_size=col_submatrix(B,rows,p2,q2)
    return C, new_row_size, new_col_size


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
    # Assumes b>0 and d>0.  Returns -1, 0 or 1 according as a/b <,=,> c/d+ 
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


# This returns 1 if A is the zero matrix, otherwise returns 0+ 
 def test_zeromat(A,m,n):
      for i in xrange(m):
          for j in xrange(n):
              if A[i][j] != 0:
                 return 0
      return 1

