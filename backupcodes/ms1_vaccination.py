#Code: Prathyush Sambaturu.  Linear Program: Anil Vullikanti
from gurobipy import *

#create the model for calibration
m = Model("Vaccination")

R = [0, 1, 2] #three regions 
T = [0, 1, 2] #three weeks. Tmax = 3
L = [0, 1, 2] # for now -> will be changed
Seed = [50, 100, 20] #Seed[i] gives the infected people in region i at time 0.
b = 10 #base value
theta = [[0.2, 0.6, 0.2], [0.6, 0.1, 0.3], [0.2, 0.3, 0.5]] #theta array nxn 
Neff = [[1000, 2000, 1800], [800, 1200, 1300], [2400, 2200, 2300]] #effective population in region i at time t


print(theta[1][2])

#Budget of a time period (week) is defined by B[week]
B = {}
for t in T:
   B[t] = 10000

guess = 3000 #bound on total number of infections.
gamma = 0.5 # rate
alpha = 0.3 # rate
beta = 0.4 # rate


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
        m.addConstr(quicksum(xv[i,t] for i in R) <= B[t], name = "RC1["+str(t)+"]")

#add Seed Infected Equality Constraints: for every region i, no. of infected people at time 0, xi[i][0] is equal to Seed[i].
for i in R:
	m.addConstr((xi[i, 0]) == Seed[i], name = "RC2["+str(i)+"]")

#add Infected Constraint: sum of no. of infected in all regions over all time is less than Obj
m.addConstr( quicksum(xi[i, t] for i in R for t in T) <= guess, name = "RC3")

#add Difference in Recovered at t+1 and t Constraints
for i in R:
   for t in [0,len(T)-2]:
       m.addConstr((xr[i, t+1] - xr[i, t] - gamma* xi[i, t]) == 0, name = "RC4["+str(i)+","+str(t)+"]")
       
#add Difference in Infected at t+1 and t Constraints
for i in R:
   for t in [0,len(T)-2]:
        m.addConstr((xi[i, t+1] - xi[i, t] - alpha* xe[i, t] + gamma* xr[i,t]) == 0, name = "RC5["+str(i)+","+str(t)+"]")

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
                z[i,j,t,l,ll] = m.addVar(lb = 0, obj = 1, name = "z("+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+")")
                m.addConstr( z[i,j,t,l,ll] <= ys[i,t,l], name = "RC10["+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+"]")
                m.addConstr( z[i,j,t,l,ll] <= yie[i,t,ll], name = "RC11["+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+"]")
  

for i in R:
    for t in T:
        m.addConstr((xie[i,t] - quicksum(theta[j][i]*xi[j,t] for j in R)) == 0, name = "RC12["+str(i)+","+str(t)+"]") 

for i in R:
    for t in [0,len(T)-2]:
        m.addConstr(xs[i,t+1] - xs[i, t] +xv[i,t] + quicksum( ((theta[i][j]*beta)/(Neff[j][t])) * quicksum( z[i,j,t,l,ll] * pow(b, l+ll) for l in L for ll in L) for j in R) == 0, name = "RC13["+str(i)+","+str(t)+"]")  
        m.addConstr(xe[i,t+1] - xe[i,t] + alpha * xe[i,t] - quicksum( (theta[i][j] *beta)/(Neff[j][t]) * quicksum( z[i,j,t,l,ll] * pow(b, l+ll) for l in L for ll in L) for j in R) == 0, name = "RC14["+str(i)+","+str(t)+"]") 

m.setObjective( quicksum(1.0 * z[i,j,t,l,ll] for i in R for j in R for t in T for l in L for ll in L) , GRB.MAXIMIZE)

m.update()

m.write('vaccination.lp')
m.write('vaccination.mps')

m.optimize()
if m.status == GRB.Status.OPTIMAL:
        print('\nobjective value: %g' % m.objVal)
else:
        print('No solution')

print(m)
