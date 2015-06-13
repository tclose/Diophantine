#
#def cholesky(A,m): # A is positive definite mxm
#  global rationum
#  global ratioden
#  global multnum
#  global multden
#  global subnum
#  global subden
#  global choleskynum
#  global choleskyden
#
#  Qnum=A
#  for i in xrange(m):
#      for j in xrange(m):
#          Qden[i][j]=1
#      
#  
#  for i in xrange(1, m):
#
#     for j in xrange(i, m):
#         Qnum[j][i]=Qnum[i][j]
#         Qden[j][i]=Qden[i][j]
#         rationum, ratioden = ratior(Qnum[i][j],Qden[i][j],Qnum[i][i],Qden[i][i])
#         Qnum[i][j]=rationum
#         Qden[i][j]=ratioden
#     
#     for k in xrange(i, m):
#         for l in xrange(k - 1, m):
#               multnum, multden = multr(Qnum[k][i],Qden[k][i],Qnum[i][l],Qden[i][l])
#               t2num=multnum
#               t2den=multden
#               subnum, subden = subr(Qnum[k][l],Qden[k][l],t2num,t2den)
#               Qnum[k][l]=subnum
#               Qden[k][l]=subden
#         
#     
#  
#  for i in xrange(m):
#      for j in xrange(m):
#          choleskynum[i][j]=Qnum[i][j]
#          choleskyden[i][j]=Qden[i][j]
#  return
#

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

