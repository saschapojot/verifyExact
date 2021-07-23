from s2Rk2Fcnt import *
from datetime import datetime

t0=datetime.now()
#calculate exact values
solutionExact=[]

for q in range(0,Q+1):
    tCurr=dt*q
    psiCurr=[psiN(n,tCurr) for n in range(0,N)]
    solutionExact.append(psiCurr)
t1=datetime.now()
print("Exact solution calculated, time = ",t1-t0)
#numerical solution
t2=datetime.now()
dataAll=[]
dataAll.append(psi0)
for q in range(0,Q):
    psiCurr=dataAll[-1]
    psiNext=oneStepS2(psiCurr)
    dataAll.append(psiNext)

t3=datetime.now()
print("numerical solution calculated, time = ", t3-t2)


tAll=[dt*q for q in range(0,Q+1)]
distAll=[l2Dist(solutionExact[q],dataAll[q]) for q in range(0,Q+1)]
# distAll=[lInfDist(solutionExact[q],dataAll[q]) for q in range(0,Q+1)]
plt.figure()
plt.plot(tAll,distAll,color="black")
plt.title("dt = "+str(dt))
plt.savefig("l2dist.png")
plt.close()


