import matplotlib.pyplot as plt
import matplotlib.patches as pat
import math as m
import sys
def yS(rho,theta):
    t = [0,556]
    re = []
    for i in t:
        ## this may ned to be jsut theta
        
        re.append(-(i*m.cos(theta)-rho)/m.sin(theta))
    return re
x,y = [],[]
"""
k = open("PrerhoPoints.txt").readlines()
xO,yO = [],[]
jf = 0
fig, ax = plt.subplots()
for i in k:
    z = i.split(",")
    jf+=1
    rho = float(z[0])
    theta = float(z[1])
    xO.append(rho)
    yO.append(theta)
print(jf)
f = open("rhoPoints.txt").readlines()
l = 0
for i in f:
    l+=1
    z = i.split(",")
     
    rho = float(z[0])
    theta = float(z[1])
    ax.add_artist(plt.Circle((rho,theta),radius=0.045,fill=False))
    
    x.append(rho)
    y.append(theta)
print(l)   

    rho = float(z[0])
    theta = float(z[1])
    if theta ==0:
        plt.plot([rho,rho],[0,1000])
    else:
        plt.plot([0,1000],yS(rho,theta))

    
j = open("tanPoints.txt").readlines()
for i in j:
    z = i.split(",")
    #plt.plot([float(z[0]),float(z[2]),float(z[4])],[float(z[1]),float(z[3]),float(z[5])])
    x.append(float(z[2]))
    y.append(float(z[3]))

sigX,sigY = [],[]
j = open("sigPoints.txt").readlines()
for i in j:
    z = i.split(" ")
    sigX.append(float(z[0]))
    sigY.append(float(z[1]))
"""
t = sys.argv[1]
f = open(t).readlines()
for i in f:
    t = i.split(" ")
    if t[0] == 'l':
        break
    else:
        x.append(float(t[1]))
        y.append(float(t[2]))
"""
lX, lY = [],[]
l = open("sigLine.txt").readlines()
for i in l:
    z = i.split("/")
    p1 = [int(i) for i in z[0].split(",")]
    p2 = [int(i) for i in z[1].split(",")]ls
    lX.append(p1[0])
    lX.append(p2[0])
    x.append(p1[0])
    x.append(p2[0])

    lY.append(p1[1])
    lY.append(p2[1])
    y.append(p1[1])
    y.append(p2[1])
"""
plt.scatter(x,y,color="b")
##plt.scatter(sigX,sigY,color="r")
mX, mY = [],[]
m = open("matching.txt").readlines()
b,f = False, False
for i in m:
    if i == "PairS\n":
        mX, mY = [],[]
    elif i[0]=="m":
        z = i.split(":")[1].split(",")
        ##plt.scatter(int(z[0]),int(z[1]),c='y')
        ##plt.scatter(int(z[2]),int(z[3]),c='g')
        ##plt.plot([int(z[0]),int(z[2])],[int(z[1]),int(z[3])])
    elif i == "PairE\n":
        ##print("set")
        if(b):
            plt.scatter(mX,mY,c='w')
        else:
            plt.scatter(mX,mY,c='m')
        mX, mY = [],[]
    elif i == "B\n":
        b = True
        f = False
    elif i == "F\n":
        b = False
        f = True
    else:
        z = i.split(",")
        mX.append(int(z[0]))
        mY.append(int(z[1]))




plt.ylim((0,max(y)+20))
plt.xlim((0,max(x)+20))

#plt.scatter(xO,yO,color="r")
##plt.plot(lX,lY)
plt.show()