#@author: Prathyush Sambaturu

from gurobipy import *
import sys
import time

start_time = time.time()

#inputs
R = [] # regions
T = [] # time periods
L = [] # powers of expansion
Neff = {} #effective population in region i at time t
theta = []
obs = []
vacc = []

fr = open('vacc_output.csv','w')

#set up all input data from the file
#=========================================
with open(sys.argv[1]) as fp:
    lines = fp.readlines()
    #print(len(lines))
    l1 = lines[0].split(",")
    noRegions = int(l1[0])
    timePeriods = int(l1[1])
    noLValues = int(l1[2])

    for i in range(0, noRegions):
        R.append(i)
    for i in range(0, timePeriods):
        T.append(i)
    for i in range(0, noLValues):
        L.append(i)
    b = int(l1[3])
    error = float(l1[4])
    print(b)
    print(R)
    print(T)
    print(L)
    print(error)

    #theta = []
    l2 = lines[4].strip("\n").split(",")
    k = 0
    print("no of values in Theta Matrix :"+str(len(l2)))
    for i in range(0, noRegions):
        theta.append([])
        for j in range(0, noRegions):
            theta[i].append(float(l2[k]))
            k = k+1
    #print(theta)

    l3 = lines[1].strip("\n").split(",")
    N = []
    for i in range(0, noRegions):
            N.append(float(l3[i]))
    #print(N)

    l4 = lines[2].strip("\n").split(",")
    alpha = float(l4[0])
    beta = float(l4[1])
    gamma = float(l4[2])

    l5 = lines[3].strip("\n").split(",")
    k = 0
    print("no of values in Obs Matrix :"+str(len(l5)))
    for i in range(0, noRegions):
        obs.append([])
        for t in range(0, timePeriods):
            obs[i].append(float(l5[k]))
            k = k+1
    #print(obs)
            
    #theta = []
    l6 = lines[5].strip("\n").split(",")
    k = 0
    print("no of values in vacc Matrix :"+str(len(l6)))
    for i in range(0, noRegions):
        vacc.append([])
        for t in range(0, timePeriods):
            vacc[i].append(int(l6[k]))
            k = k+1
    #print(vacc)

#Creating the model for calibration
m = Model("Calibration")

#variables depending on only time and region i.e., of type x[i,t]
xe = {} #number of exposed in region i at time t
xs = {} #number of susceptible in region i at time t
xr = {} #number of recovered in region i at time t
xi = {} #number of infected in region i at time t
xie = {} #number of effectively infected in region i at time t
delta = {} #difference between observed and actual infected in region i at time t
seed0 = {} #seed infected in region i at time 0

for i in R:
    for t in T:
        xe[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xv("+str(i)+","+str(t)+")")
        xs[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xs("+str(i)+","+str(t)+")")
        xr[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xr("+str(i)+","+str(t)+")")
        xi[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xi("+str(i)+","+str(t)+")")
        xie[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xie("+str(i)+","+str(t)+")")
        delta[i, t] = m.addVar(lb = 0, obj = 0.0, name ="delta("+str(i)+","+str(t)+")")

for i in R:
    seed0[i] = m.addVar(lb = 0, obj = 0.0, name = "seed0("+str(i)+")")

#variables depending on time, region, and l values
ys = {}
yie = {}

for i in R:
    for t in T:
        for l in L:
            ys[i, t, l] = m.addVar(vtype = GRB.BINARY, name = "ys("+str(i)+","+str(t)+","+str(l)+")")
            yie[i, t, l] = m.addVar(vtype = GRB.BINARY, name = "yie("+str(i)+","+str(t)+","+str(l)+")")


#constraints on expansions of xs and xie
for i in R:
   for t in T:
      m.addConstr( quicksum(ys[i,t,l]* pow(b, l) for l in L) >= xs[i,t] - 1, name = "RC1["+str(i)+","+str(t)+"]")
      m.addConstr( quicksum(ys[i,t,l]* pow(b, l) for l in L) <= xs[i,t], "RC2["+str(i)+","+str(t)+"]")
      m.addConstr( quicksum(yie[i,t,l]* pow(b, l) for l in L) >= xie[i,t] - 1, name = "RC3["+str(i)+","+str(t)+"]")
      m.addConstr( quicksum(yie[i,t,l]* pow(b, l) for l in L) <= xie[i,t], name = "RC4["+str(i)+","+str(t)+"]")


#variable z in the objective
z = {}

for i in R:
    for j in R:
       for t in T:
          for l in L:
             for ll in L:
                z[i,j,t,l,ll] = m.addVar(lb = 0, obj = 1.0, name = "z("+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+")")
                m.addConstr( z[i,j,t,l,ll] <= ys[i,t,l], name = "RC5["+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+"]")
                m.addConstr( z[i,j,t,l,ll] <= yie[i,t,ll], name = "RC6["+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+"]")


#print(len(T))
#difference in number of people in each compartment at time t+1 and t
for i in R:
   for t in range(0, len(T)-1):
       m.addConstr((xr[i, t+1] - xr[i, t] - gamma* xi[i, t]) == 0, name = "RC7["+str(i)+","+str(t)+"]")
       m.addConstr((xi[i, t+1] - xi[i, t] - alpha* xe[i, t] + gamma* xi[i,t]) == 0, name = "RC8["+str(i)+","+str(t)+"]")

#effective infected equation
for i in R:
    for t in T:
        m.addConstr( quicksum(theta[j][i]*xi[j,t] for j in R) == xie[i,t], name = "RC9["+str(i)+","+str(t)+"]")

#erro related constraints
for i in R:
    for t in T:
        m.addConstr( (delta[i, t] - xi[i, t] + obs[i][t]) >= 0, name = "RC10["+str(i)+","+str(t)+"]")
        m.addConstr( (delta[i, t] + xi[i, t] - obs[i][t]) >= 0, name = "RC11["+str(i)+","+str(t)+"]")

m.addConstr( quicksum(delta[i, t] for i in R for t in T) <= error, name = "RC12")

#seed constraints
for i in R:
    m.addConstr( xi[i, 0] == seed0[i], name = "RC13["+str(i)+"]")

#computing Neff in region i at time t
for i in R:
    for t in T:
        if t==0:
           Neff[i, t] = N[i]
        else:
           Neff[i, t] = 0

#print(Neff)
for j in R:
    for t in range(1, len(T)):
        #print(j, t)
        for i in R:
            Neff[j, t] += theta[i][j] * Neff[i, t-1]

#print(Neff)

#difference in susceptible and exposed constraints
for i in R:
    for t in range(0,len(T)-1):
        m.addConstr(quicksum( ((theta[i][j]*beta)/Neff[j, t]) * quicksum( z[i,j,t,l,ll] * pow(b, l+ll) for l in L for ll in L) for j in R) == xs[i,t] - xs[i,t+1] -vacc[i][t] , name = "RC14["+str(i)+","+str(t)+"]")
        m.addConstr(quicksum( ((theta[i][j] *beta)/Neff[j, t]) * quicksum( z[i,j,t,l,ll] * pow(b, l+ll) for l in L for ll in L) for j in R) == -xe[i,t] + xe[i,t+1] + alpha *xe[i,t], name = "RC15["+str(i)+","+str(t)+"]")

#setup objective and update the program
m.setObjective( quicksum(1.0 * z[i,j,t,l,ll] for i in R for j in R for t in T for l in L for ll in L) , GRB.MAXIMIZE)
m.update()
m.write('calibration.lp')
m.write('calibration.mps')


m.optimize()

elapsed_time = time.time() - start_time
print("#################################################")
print("Total time elapsed: "+str(elapsed_time))

