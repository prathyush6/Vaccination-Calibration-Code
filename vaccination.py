# @author: Prathyush Sambaturu.

from gurobipy import *
import sys
import time

start_time = time.time()

#create the model for vaccination
m = Model("Vaccination")

#R = [0, 1, 2] #Regions
#T = [0, 1, 2] #Time Periods
#L = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] #L values
#N = [20000, 15000, 18000] #Population of R regions at time 0
#B = [1000, 1000, 1000] # Budget at time t
#theta = [ [0.7, 0.1, 0.2], [0.2, 0.5, 0.3], [0.1, 0.3, 0.6] ] #the fraction of individuals from region i at region j
#b = 2
R = []
T = []
L = []
N = []
B = []
theta = []
noRegions = 0
timePeriods = 0
noLValues = 0
guess = 100000
Seed = [] #Seed infected in region i at time 0
alpha = 0.4
beta = 0.3
gamma = 0.5

with open(sys.argv[1]) as fp:
     lines = fp.readlines()
     l1 = lines[0].split(",")
     noRegions = int(l1[0])
     timePeriods = int(l1[1])
     noLValues = int(l1[2])
     b = int(l1[3])

     for i in range(0, noRegions):
         R.append(i)
     for i in range(0, timePeriods):
         T.append(i)
     for i in range(0, noLValues):
         L.append(i)
     print(R)
    
     l5 = lines[4].strip("\n").split(",")
     k = 0
     for i in range(0, noRegions):
         theta.append([])
         for j in range(0, noRegions):
             theta[i].append(float(l5[k]))
             k = k+1
     
     l2 = lines[1].strip("\n").split(",")
     N = []
     for i in range(0, noRegions):
         N.append(float(l2[i]))

     l6 = lines[5].strip("\n").split(",")

     for t in range(0, timePeriods):
         B.append(float(l6[t]))
    
     guess = int(lines[6].strip("\n"))
     print(guess)
   
     l3 = lines[2].strip("\n").split(",")
     alpha = float(l3[0]) 
     beta = float(l3[1])
     gamma = float(l3[2])
  
     l4 = lines[3].strip("\n").split(",")
          
     for i in range(0, noRegions):
         Seed.append(float(l4[i]))
   



#computing Neff
Neff = []
for j in R:
    eff_pop = 0
    for i in R:
        eff_pop = eff_pop + theta[i][j] * N[i]
    Neff.append(eff_pop)


#variables depending only on time and region i.e of type x[i, t]
xv = {} #vaccined in region i at time t
xe = {} #exposed in region i at time t
xs = {} #suscceptible in region i at time t
xi = {} #infected in region i at time t
xr = {} #recovered in region i at time t
xie = {} #effective infected in region i at time t

for i in R:
    for t in T:
        xv[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xv("+str(i)+","+str(t)+")")
        xe[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xe("+str(i)+","+str(t)+")")
        xs[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xs("+str(i)+","+str(t)+")")
        xr[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xr("+str(i)+","+str(t)+")")
        xi[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xi("+str(i)+","+str(t)+")")
        xie[i, t] = m.addVar(lb = 0, obj = 0.0, name = "xie("+str(i)+","+str(t)+")")

#budget constraints
for t in T:
    m.addConstr( quicksum(xv[i, t] for i in R) <= B[t], name = "RC1["+str(t)+"]")

#seed infected constraints 
for i in R:
    m.addConstr( (xi[i, 0]) == Seed[i], name = "RC2["+str(i)+"]")
    
# sum of infected in all regions over all time is less than guess
m.addConstr( quicksum(xi[i, t] for i in R for t in T) <= guess, name = "RC3")

#recovered constraints
for i in R:
    for t in range(0, len(T)-1):
        m.addConstr( (xr[i, t+1] - xr[i, t] - gamma*xi[i,t] ) == 0, name = "RC4["+str(i)+","+str(t)+"]")

#infected constriants
for i in R:
    for t in range(0, len(T)-1):
        m.addConstr( (xi[i, t+1] - xi[i, t] - alpha* xe[i, t] + gamma* xi[i, t] ) == 0, name = "RC5["+str(i)+","+str(t)+"]")

ys = {}
yie = {}

for i in R:
    for t in T:
        for l in L:
            ys[i, t, l] = m.addVar(vtype = GRB.BINARY, name = "ys("+str(i)+","+str(t)+","+str(l)+")")
            yie[i, t, l] = m.addVar(vtype = GRB.BINARY, name = "yie("+str(i)+","+str(t)+","+str(l)+")")

 
for i in R:
    for t in T:   
        m.addConstr( quicksum(ys[i, t, l]* pow(b, l) for l in L) >= xs[i, t]-1, name = "RC6["+str(i)+","+str(t)+"]")
        m.addConstr( quicksum(ys[i, t, l]* pow(b, l) for l in L) <= xs[i, t], name ="RC7["+str(i)+","+str(t)+"]")
        m.addConstr( quicksum(yie[i, t, l]* pow(b, l) for l in L) >= xie[i, t] - 1, name = "RC8["+str(i)+","+str(t)+"]")
        m.addConstr( quicksum(yie[i, t, l]* pow(b, l) for l in L) <= xie[i, t], name = "RC9["+str(i)+","+str(t)+"]")

z = {}

for i in R:
    for j in R:
        for t in T:
            for l in L:
                for ll in L:
                    z[i, j, t, l, ll] = m.addVar(lb = 0, name = "z("+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+")")
                    m.addConstr( z[i, j, t, l, ll] <= ys[i, t, l], name= "RC10["+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+")")
                    m.addConstr( z[i, j, t, l, ll] <= yie[i, t, ll], name = "RC11["+str(i)+","+str(j)+","+str(t)+","+str(l)+","+str(ll)+")")

for i in R:
    for t in T:
        m.addConstr( quicksum(theta[j][i] * xi[j, t] for j in R) <= xie[i, t], name = "RC12["+str(i)+","+str(t)+"]")

for i in R:
    for t in range(0, len(T)-1):
        m.addConstr( quicksum( ((theta[i][j]*beta)/Neff[j]) * quicksum( z[i, j, t, l, ll] * pow(b, l+ll) for l in L for ll in L) for j in R) == xs[i, t] - xs[i,t+1] - xv[i,t], name = "RC13["+str(i)+","+str(t)+"]")
        m.addConstr( quicksum( ((theta[i][j]*beta)/Neff[j]) * quicksum( z[i, j, t, l, ll] * pow(b, l+ll) for l in L for ll in L) for j in R) == -xe[i, t] + xe[i, t+1] + alpha * xe[i, t], name = "RC14["+str(i)+","+str(t)+"]")

m.setObjective( quicksum(1.0 * z[i, j, t, l, ll] for i in R for j in R for t in T for l in L for ll in L), GRB.MAXIMIZE)

gap = 5
m.setParam(GRB.Param.MIPGap, gap)

m.update()
m.write('vaccination.lp')


m.optimize()



elapsed_time = time.time() - start_time
print("Total time elapsed: "+str(elapsed_time)) 

       










