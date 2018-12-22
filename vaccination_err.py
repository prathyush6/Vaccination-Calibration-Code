#@Author: Prathyush Sambaturu.
from gurobipy import *
import sys
import time

start_time = time.time()

#create the model for calibration
m = Model("Vaccination")

R = [] #Regions
T = [] #Time Periods
L = [] #L values
N = [] #Population of R regions at time 0
B = [] #Budget at time t
noRegions = 0
timePeriods = 0
noLValues = 0
guess = 0
Seed = [] #Seed Infected in each region

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
    l2 = lines[1].strip("\n").split(",")
    N = []
    for i in range(0, noRegions):
            N.append(float(l2[i]))
    #print(N)  
    
    l3 = lines[5].strip("\n").split(",")
    
    for t in range(0, timePeriods):
        B.append(float(l3[t]))
    #print(B) 

    guess = int(lines[6].strip("\n"))
    #print(guess)

    l4 = lines[2].strip("\n").split(",")
    alpha = float(l4[0])
    beta = float(l4[1])
    gamma = float(l4[2])

    l7 = lines[3].strip("\n").split(",")
    Seed = []
    for i in range(0, noRegions):
       Seed.append(float(l7[i]))
    #print(alpha, beta, gamma)
    #print(Seed)
    #print(theta[2][2])
#===i================================
#set up Neff
Neff = []
for j in range(0, noRegions):
    eff_population = 0
    for i in range(0, noRegions):
        eff_population = eff_population+ theta[i][j] * N[i]
    Neff.append(eff_population)
#print(Neff)

#Seed_eff = []
#for j in range(0, noRegions):
#    eff_infected = 0
#    for i in range(0, noRegions):
#        eff_infected = eff_infected + theta[i][j] * Seed[i]
#    Seed_eff.append(eff_infected)
#print(Seed_eff)

#setup list K: aggregated regions
noK = int(sys.argv[2])
K = []
for i in range(0, noK):
    K.append(i)  
print(K)

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
    #m.addConstr((xie[i,0]) == Seed_eff[i], name = "RC17["+str(i)+"]")

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

lot = []
f = int(noRegions/noK)
for i in range(0,f):
    lot.append(i)
print(lot)

ys = {}
yie = {}
for k in K:
    for t in T:
       for l in L:
       	   ys[k, t, l] = m.addVar(lb=0, name = "ys("+str(k)+","+str(t)+","+str(l)+")")
           yie[k, t, l] = m.addVar(ub=0, name = "yie("+str(k)+","+str(t)+","+str(l)+")")

for k in K:
   for t in T: 
       m.addConstr( quicksum(ys[k,t,l]* pow(b, l) for l in L) - quicksum(xs[f*k+i,t] for i in lot) >= - 10000, name = "RC6["+str(k)+","+str(t)+"]")	
       m.addConstr( quicksum(ys[k,t,l]* pow(b, l) for l in L) <= quicksum(xs[f*k+i,t] for i in lot), "RC7["+str(k)+","+str(t)+"]")
       m.addConstr( quicksum(yie[k,t,l]* pow(b, l) for l in L) - quicksum(xie[f*k+i,t] for i in lot) >=  -10000, name = "RC8["+str(k)+","+str(t)+"]")
       m.addConstr( quicksum(yie[k,t,l]* pow(b, l) for l in L) <= quicksum(xie[f*k+i,t] for i in lot), name = "RC9["+str(k)+","+str(t)+"]")

z = {}

for k in K:
    for kk in K:
       for t in T:
          for l in L:
             for ll in L:
                z[k,kk,t,l,ll] = m.addVar(lb = 0, name = "z("+str(k)+","+str(kk)+","+str(t)+","+str(l)+","+str(ll)+")")
                m.addConstr( z[k,kk,t,l,ll] <= ys[k,t,l], name = "RC10["+str(k)+","+str(kk)+","+str(t)+","+str(l)+","+str(ll)+"]")
                m.addConstr( z[k,kk,t,l,ll] <= yie[kk,t,ll], name = "RC11["+str(k)+","+str(kk)+","+str(t)+","+str(l)+","+str(ll)+"]")
  

for i in R:
    for t in T:
        #print(i,t)
        m.addConstr( quicksum(theta[j][i]*xi[j,t] for j in R) >= xie[i,t], name = "RC12["+str(i)+","+str(t)+"]") 

for i in R:
    for t in range(0,len(T)-1):
        k = int(i/f)
        m.addConstr(quicksum( ((theta[i][j]*beta)/Neff[j]) * quicksum( z[k,int(j/f),t,l,ll] * pow(b, l+ll) for l in L for ll in L) for j in R)/100.0 == xs[i,t] - xs[i,t+1] -xv[i,t] , name = "RC13["+str(i)+","+str(t)+"]")  
        m.addConstr(quicksum( ((theta[i][j] *beta)/Neff[j]) * quicksum( z[k,int(j/f),t,l,ll] * pow(b, l+ll) for l in L for ll in L) for j in R)/100.0 == -xe[i,t] + xe[i,t+1] + alpha *xe[i,t], name = "RC14["+str(i)+","+str(t)+"]") 

m.setObjective( quicksum(1.0 * z[k,kk,t,l,ll] for k in K for kk in K for t in T for l in L for ll in L) , GRB.MAXIMIZE)
#m.Params.IterationLimit = 100000
#m.Params.NodeLimit = 20
#m.Params.BarConvTol = 0.5
#m.Params.Cutoff = 42000

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
        fr.write("Infected :"+str(infected)+"\n\n")
        exposed = m.getAttr('X',xe)
        fr.write("Exposed :"+str(exposed)+"\n\n")
        vaccinated = m.getAttr('X',xv)
        fr.write("Vaccinated :"+str(vaccinated)+"\n\n")
        recovered = m.getAttr('X',xr)
        fr.write("Recovered :"+str(recovered)+"\n\n")
        vaccinated = m.getAttr('X', xv)
        #print("Vaccinated :"+str(vaccinated))
        #zees = m.getAttr('X', z)
        #print("Z values :"+str(zees))
else:
        #infected = m.getAttr('X', xi)
        #fr.write("Infected :"+str(infected)+"\n\n")
        #vaccinated = m.getAttr('X',xv)
        #fr.write("Vaccinated :"+str(vaccinated)+"\n\n")
        fr.write('No solution')


fp = open('vacc.csv','w') 
for i in R:
    for t in T:
        if (t == len(T)-1):
           fp.write(str(vaccinated[i,t])+"\n")
        else:
           fp.write(str(vaccinated[i,t])+",")
fp.close()

fp = open('infected.csv', 'w')
for i in R:
    for t in T:
        if(t == len(T)-1):
           fp.write(str(infected[i,t])+"\n")
        else:
           fp.write(str(infected[i,t])+",")
fp.close()

#zee = m.getAttr('X',z)
#obj_current = 0.0
#for i in R:
#    for j in R:
#        for t in T:
#            for l in L:
#                for ll in L:
#                    obj_current = obj_current + zee[i,j,t,l,ll]


#print("Objective Value: "+str(obj_current))

elapsed_time = time.time() - start_time
print("#############################################")
print("Total time elapsed: "+str(elapsed_time))
fr.write("#############################################")
fr.write("Total time elapsed: "+str(elapsed_time))
