#Code: Prathyush Sambaturu.  Linear Program: Anil Vullikanti
from gurobipy import *
import sys
import time

start_time = time.time()
# your code
#elapsed_time = time.time() - start_time
#create the model for calibration
m = Model("Vaccination")

#R = [0, 1, 2] #three regions 
#T = [0, 1] #three weeks. Tmax = 3
#L = [0, 1, 2] # for now -> will be changed
#Seed = [1500, 1000, 2200] #Seed[i] gives the infected people in region i at time 0.
#b = 10 #base value
#theta = [[0.2, 0.6, 0.2], [0.6, 0.1, 0.3], [0.2, 0.3, 0.5]] #theta array nxn 
#Neff = [[20000, 23200], [10000, 12000], [24000, 20000]] #effective population in region i at time t

#guess = 1000000 #bound on total number of infections.
#gamma = 0.4 # recovery rate
#alpha = 0.3 # rate
#beta = 0.3 # rate

#set up all input data from the file
#=========================================
with open(sys.argv[1]) as fp:  
    lines = fp.readlines()
    #print(len(lines))
    l1 = lines[0].split(",")
    noRegions = int(l1[0])
    timePeriods = int(l1[1])
    noLValues = int(l1[2])
    #print(str(noRegions)+"adaf")
    R = []
    T = []
    L = []
    for i in range(0, noRegions):
        R.append(i)
    for i in range(0, timePeriods):
        T.append(i)
    for i in range(0, noLValues):
        L.append(i)
    b = int(l1[3])
    #print(b)
    #print(R)
    #print(T)
    #print(L)
    theta = []
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
    k = 0
    Neff = []
    for i in range(0, noRegions):
            Neff.append(float(l3[0]))
            k = k+1
    #print(Neff)    
    l4 = lines[5].strip("\n").split(",")
    k = 0 
    B = []
    for t in range(0, timePeriods):
        B.append(float(l4[0]))
        k = k +1 
    #print(B) 
    guess = int(lines[6].strip("\n"))
    #print(guess)
    l6 = lines[2].strip("\n").split(",")
    alpha = float(l6[0])
    beta = float(l6[1])
    gamma = float(l6[2])
    l7 = lines[3].strip("\n").split(",")
    Seed = []
    for i in range(0, noRegions):
       Seed.append(float(l7[0]))
    print(alpha, beta, gamma)
    #print(Seed)
    #print(theta[2][2])
#===i================================


#variables depending on only time and region i.e of type x[i,t]
xv = {} #xv[i,t] is variable denoting the # of vaccinated in region i at time t
xe = {} #xe -> exposed
xs = {} #xs ->susceptible
xr = {} #xr -> recovered
xi = {} #xi -> infected
xie = {} #xie -> effetively infected

for i in R:
   for t in T:
        xv[i, t] = m.addVar(lb = 0, obj = 0, name = "xv("+str(i)+","+str(t)+")")
        xe[i, t] = m.addVar(lb = 0, obj = 0, name = "xe("+str(i)+","+str(t)+")")
        xs[i, t] = m.addVar(lb = 0, obj = 0, name = "xs("+str(i)+","+str(t)+")")
        xr[i, t] = m.addVar(lb = 0, obj = 0, name = "xr("+str(i)+","+str(t)+")")
        xi[i, t] = m.addVar(lb = 0, obj = 0, name = "xi("+str(i)+","+str(t)+")")
        xie[i, t] = m.addVar(lb = 0, obj = 0, name = "xie("+str(i)+","+str(t)+")")

#add Budget Constraints: for every time t, total no. of vaccinated in all regions is less than or equal to Budget of that time period
for t in T:
        m.addConstr(quicksum(xv[i,t] for i in R) <= B[t],name = "RC1["+str(t)+"]")

#add Seed Infected Equality Constraints: for every region i, no. of infected people at time 0, xi[i][0] is equal to Seed[i].
for i in R:
    m.addConstr((xi[i, 0]) == Seed[i], name = "RC2["+str(i)+"]")
    m.addConstr((xr[i, 0]) == 0.0, name = "RC15["+str(i)+"]")
    m.addConstr((xe[i, 0]) == 0.0, name = "RC16["+str(i)+"]")
    #m.addConstr((xs[i,0]) == 4000.0, name = "RC17["+str(i)+"]")
    #m.addConstr((xv[i,0]) == 0.0, name = "RC18["+str(i)+"]")
    #m.addConstr((xie[i,0]) == 0.0, name = "RC19["+str(i)+"]")

#add Infected Constraint: sum of no. of infected in all regions over all time is less than Obj
m.addConstr( quicksum(xi[i, t] for i in R for t in T) <= guess, name = "RC3")

#add Difference in Recovered at t+1 and t Constraints
for i in R:
   for t in range(0,len(T)-1):
       m.addConstr((xr[i, t+1] - xr[i, t] - gamma* xi[i, t]) == 0, name = "RC4["+str(i)+","+str(t)+"]")
       
#add Difference in Infected at t+1 and t Constraints
for i in R:
   for t in range(0,len(T)-1):
        m.addConstr((xi[i, t+1] - xi[i, t] - alpha* xe[i, t] + gamma* xi[i,t]) == 0, name = "RC5["+str(i)+","+str(t)+"]")

ys = {}
yie = {}
for i in R:
    for t in T:
       for l in L:
       	   ys[i, t, l] = m.addVar(vtype=GRB.BINARY, name = "ys("+str(i)+","+str(t)+","+str(l)+")")
           yie[i, t, l] = m.addVar(vtype=GRB.BINARY, name = "yie("+str(i)+","+str(t)+","+str(l)+")")

for i in R:
   for t in T: 
      m.addConstr( quicksum(ys[i,t,l]* pow(b, l) for l in L) >= xs[i,t] - 1, name = "RC6["+str(i)+","+str(t)+"]")	
      m.addConstr( quicksum(ys[i,t,l]* pow(b, l) for l in L) <= xs[i,t], "RC7["+str(i)+","+str(t)+"]")
      m.addConstr( quicksum(yie[i,t,l]* pow(b, l) for l in L) >= xie[i,t]-1, name = "RC8["+str(i)+","+str(t)+"]")
      m.addConstr( quicksum(yie[i,t,l]* pow(b, l) for l in L) <= xie[i,t], name = "RC9["+str(i)+","+str(t)+"]")

z = {}

for i in R:
    for j in R:
       for t in T:
          for l in L:
             for ll in L:
                z[i,j,t,l,ll] = m.addVar(lb = 0, name = "z("+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+")")
                m.addConstr( z[i,j,t,l,ll] <= ys[i,t,l], name = "RC10["+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+"]")
                m.addConstr( z[i,j,t,l,ll] <= yie[i,t,ll], name = "RC11["+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+"]")
  

for i in R:
    for t in T:
        #print(i,t)
        m.addConstr( quicksum(theta[j][i]*xi[j,t] for j in R) >= xie[i,t], name = "RC12["+str(i)+","+str(t)+"]") 

for i in R:
    for t in range(0,len(T)-1):
        m.addConstr(quicksum( ((theta[i][j]*beta)/Neff[j]) * quicksum( z[i,j,t,l,ll] * pow(b, l+ll) for l in L for ll in L) for j in R) == xs[i,t] - xs[i,t+1] -xv[i,t] , name = "RC13["+str(i)+","+str(t)+"]")  
        m.addConstr(quicksum( ((theta[i][j] *beta)/Neff[j]) * quicksum( z[i,j,t,l,ll] * pow(b, l+ll) for l in L for ll in L) for j in R) == -xe[i,t] + xe[i,t+1] + alpha *xe[i,t], name = "RC14["+str(i)+","+str(t)+"]") 

m.setObjective( quicksum(1.0 * z[i,j,t,l,ll] for i in R for j in R for t in T for l in L for ll in L) , GRB.MAXIMIZE)

m.update()

m.write('vaccination.lp')
m.write('vaccination.mps')

m.optimize()
if m.status == GRB.Status.OPTIMAL:
        print('\nobjective value: %g' % m.objVal)
        #print("Infected:")
        #print(xi)
        #print("Recovered:")
        #print(xr)
        #print("Vaccinated:")
        #print(xv)
        infected = m.getAttr('X',xi)
        #print("Infected :"+str(infected))
        exposed = m.getAttr('X',xe)
        #print("Exposed :"+str(exposed))
        vaccinated = m.getAttr('X',xv)
        #print("Vaccinated :"+str(vaccinated))
        recovered = m.getAttr('X',xr)
        #print("Recovered :"+str(recovered))
        susceptible = m.getAttr('X', xs)
        #print("Susceptible :"+str(susceptible))
        #zees = m.getAttr('X', z)
        #print("Z values :"+str(zees))
else:
        print('No solution')

elapsed_time = time.time() - start_time
print("#############################################")
print("Total time elapsed: "+str(elapsed_time))
