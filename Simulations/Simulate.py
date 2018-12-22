import sys

R = []
T = []
N = []

fr = open('sim_output.txt','w')


#set up all input data from the file
#=========================================
with open(sys.argv[1]) as fp:  
    lines = fp.readlines()
    
    
    l1 = lines[0].split(",")
    noRegions = int(l1[0])
    timePeriods = int(l1[1])
        
    for i in range(0, noRegions):
        R.append(i)
    for i in range(0, timePeriods):
        T.append(i)
    print("R: "+str(R))
    print("T: "+str(T))
    
    l2 = lines[1].strip("\n").split(",")
    k = 0

    N = []
    for i in range(0, noRegions):
        N.append(int(l2[i]))
        k = k+1
    #print("N: "+str(N))

    l3 = lines[2].strip("\n").split(",")
    alpha = float(l3[0])
    beta = float(l3[1])
    gamma = float(l3[2])
    #print("a,b,g: "+str(alpha)+","+str(beta)+","+str(gamma))

    l4 = lines[3].strip("\n").split(",")
    I = []
    for i in range(0, noRegions):
        I.append([])
        I[i].append(int(l4[i]))
        for t in range(1, timePeriods):
            I[i].append(0)
      
    #print("I: "+str(I))

    #l5 = lines[4].strip("\n").split(",")
    #V = []
    #k = 0
    #for i in range(0, noRegions):
    #    V.append([])
    #    for t in range(0,timePeriods):
    #        V[i].append(float(l5[k]))
    #        k = k+1
    #print("V: "+str(V))

    theta = []
    l6 = lines[4].strip("\n").split(",")
    k = 0
    for i in range(0, noRegions):
        theta.append([])
        for j in range(0, noRegions):
            theta[i].append(float(l6[k]))
            k = k+1
    #print("theta: "+str(theta))

fq = open(sys.argv[2],'r')

lines = fq.readlines()

V = []
for i in range(0,noRegions):
    l = lines[i].strip("\n").split(",")
    V.append([])
    for t in range(0,timePeriods):
        V[i].append(float(l[t]))
#print("V:"+str(V))

        


#intialize seed values for V, E, R, and S
#=======================================
E = []
R = []
S = []

for i in range(0, noRegions):
    E.append([])
    R.append([])
    S.append([])
    #V.append([])
    for t in range(0, timePeriods):
        E[i].append(0)
        R[i].append(0)
        S[i].append(0)
        #V[i].append(0.1 * I[i][0])

#print(E)
#print(V)

#initializing susceptibles as population - infected
for i in range(0, noRegions):
    S[i][0] = N[i] - I[i][0]

#print(S)     

#Neff and Ieff Calculations
#=====================================
Neff = []
for i in range(0, noRegions):
    Neff.append(0)
    for j in range(0, noRegions):
        Neff[i] = Neff[i] + N[j] * theta[j][i]

#print("Neff:\t"+str(Neff))


Ieff = []
for i in range(0, noRegions):
    Ieff.append([])
    for t in range(0, timePeriods):
        Ieff[i].append(0)
        Ieff[i][t] = 0
        for j in range(0, noRegions):
            Ieff[i][t] = Ieff[i][t] + I[j][t] * theta[j][i]

#print("Ieff:\t"+str(Ieff)) 

print("Simulation Begins")

#equations for simulation rounds
#======================================
for i in range(0, noRegions):
    for t in range(1, timePeriods):
        R[i][t] = R[i][t-1] + gamma * I[i][t-1]
        I[i][t] = I[i][t-1] + alpha * E[i][t-1] - gamma * I[i][t-1]
    print("x")            
for i in range(0, noRegions):
    for t in range(1, timePeriods):
        for j in range(0, noRegions):
            Ieff[i][t] = Ieff[i][t] + I[j][t] * theta[j][i]
    print("x")
for i in range(0, noRegions):
    for t in range(1, timePeriods):
        term = 0
        for j in range(0, noRegions):
            term = term + theta[i][j] * beta * (Ieff[j][t-1]/ Neff[j]) * S[i][t-1]
        E[i][t] = E[i][t-1] - alpha * E[i][t-1] + term
        S[i][t] = S[i][t-1] - V[i][t] - term
    print("x")
#fr.write("\n\nR:\t"+str(R)+"\n\n")
#fr.write("I:\t"+str(I)+"\n\n")

for i in range(0, noRegions):
    for t in range(0, timePeriods):
        if(t == timePeriods - 1):
           fr.write(str(I[i][t])+"\n")
        else:
           fr.write(str(I[i][t])+",")

#fr.write("Ieff:\t"+str(Ieff)+"\n\n")
#fr.write("E:\t"+str(E)+"\n\n")
#fr.write("S:\t"+str(S))
fr.close()
print("Simulation Ends")
	
