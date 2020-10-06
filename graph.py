import matplotlib.pyplot as plt
import math as m
def yS(rho,theta):
    t = [0,556]
    re = []
    for i in t:
        ## this may ned to be jsut theta
        
        re.append(-(i*m.cos(theta)-rho)/m.sin(theta))
    return re
k = open("PrerhoPoints.txt").readlines()
xO,yO = [],[]
jf = 0
for i in k:
    z = i.split(",")
    jf+=1
    rho = float(z[0])
    theta = float(z[1])
    xO.append(rho)
    yO.append(theta)
print(jf)
f = open("rhoPoints.txt").readlines()
x,y = [],[]
l = 0
for i in f:
    l+=1
    z = i.split(",")
     
    rho = float(z[0])
    theta = float(z[1])
    x.append(rho)
    y.append(theta)
print(l)   
""" 
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

f = open("butterfly.obj").readlines()
for i in f:
    t = i.split(" ")
    if t[0] == 'l':
        break
    else:
        x.append(float(t[1]))
        y.append(float(t[2]))

plt.ylim((0,600))
plt.xlim((0,600))
"""
plt.scatter(xO,yO,color="r")
plt.scatter(x,y,color="b")
plt.show()