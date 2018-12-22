import sys
import random

theta = []

R = int(sys.argv[1])
T = int(sys.argv[2])
L = int(sys.argv[3])
b = int(sys.argv[4])
print(R,T)


N = []
for i in range(0, R):
    x = random.randint(3000,4000)
    N.append(x)

Seed = []
for i in range(0, R):
    x = 0
    toss = random.randint(0,2)
    if(toss == 1):
       x = random.randint(0,200)
    Seed.append(x)

#fp = open('vacc.csv', 'r')
#lines = fp.readlines()

#V = []
#for i in range(0,R):
#    l = lines[i].strip("\n").split(",")
#    V.append([]) 
    #x = random.randint(0,500)
#    for t in range(0,T):
#        V[i].append(l[t])

theta = []
for i in range(0,R):
    theta.append([])
    x = random.randint(5,9)/10.0
    y = (1-x)/2
    if(i+2 <= R-1):
      loc1 = i+2
    else:
      loc1 = i-2
    if(i-4 >= 0):
      loc2 = i-4
    else:
      loc2 = i+4
    #print(x,y)
    for j in range(0,R):
          if (i==j):
             theta[i].append(x)
          elif (j==loc1 or j == loc2):
             theta[i].append(y)
          else:
             theta[i].append(0)
fp = open('simulatedip.csv','w') 
fp.write(str(R)+","+str(T)+","+str(L)+","+str(b)+"\n")

for i in range(0,R):
    val = str(N[i])
    if(i == R-1):
       fp.write(val)
    else:
       fp.write(val+",")
fp.write("\n")
fp.write("0.3,0.3,0.5\n")

for i in range(0,R):
    val = str(Seed[i])
    if(i == R-1):
      fp.write(val)
    else:
      fp.write(val+",")
fp.write("\n")



#for i in range(0,R):
#    for t in range(0,T):
#        val = str(V[i][t])
#        if(i == R-1 and t == R-1):
#           fp.write(val)
#        else:
#           fp.write(val+",")
#fp.write("\n")


for i in range(0,R):
    for j in range(0,R):
        val = str(theta[i][j])
        if(i== R-1 and j == R-1):
           fp.write(val)
        else:
           fp.write(val+",")
fp.write("\n")

for i in range(0,T):
    if(i == T-1):
      fp.write("1000")
    else:
      fp.write("1000,")
fp.write("\n")
fp.write("5000000\n")


#fp.write("\n")       
#fp.write("10000\n")
#fp.write("8000\n")  
#fp.write("50000000\n")
#fp.write("0.3,0.3,0.5\n")
#fp.write("1500")

