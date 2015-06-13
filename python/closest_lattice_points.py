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

