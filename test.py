import numpy as np
from itertools import permutations,combinations
x = np.random.randint(1,10,size=3);

b = np.random.randint(1,10,size=3);
permX = permutations(x)
permB = permutations(b)

for i in list(permX):
    for j in list(permB):
        print(f"x perm: {i}")
        print(f"b perm: {j}")
        temp = 1
        for t in range(3):
            temp *= i[t]**j[t]
            #print(i[t]**j[t])
        print(temp)