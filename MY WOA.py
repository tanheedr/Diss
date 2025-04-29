import time
import random
import math
from mpmath import iv

iv.dps = 60
iv.pretty = False

def compute_f(point):
    return math.exp(math.sin(50*point[0])) + math.sin(60*math.exp(point[1])) + math.sin(70*math.sin(point[0])) + math.sin(math.sin(80*point[1])) - math.sin(10*point[0]+10*point[1]) + 0.25*(point[0]**2 + point[1]**2)

n = 2

min_range = [-10, -10]
max_range = [10, 10]

start_time = time.time()


# WO (Whale Optimization)
#   - 'n' = dimension of optimization problem           int
#   - 'num_cyc' = number of cycles to iterate           int
#   - 'N' = number of whales                            int
#   - 'b' = spiral constant                             float or int

# In summary, before you input the bulk of your code, ensure that you:
# - import any (legal) modules you wish to use in the space provided below 
# - initialize your parameters in the space provided below
# - ensure that reserved variables have the correct type on termination.

num_cyc = 1000
N = 50
b = 0.8
print(num_cyc, N, b)

#Whale class with random initial positions
class whale:
    def __init__(self):
        self.x = [random.uniform(x[0], x[1]) for x in zip(min_range, max_range)]
        self.f = compute_f(self.x)
    
    def fit(self):
        self.f = compute_f(self.x)

#Size function
def size(a):
    return math.sqrt(sum([x ** 2 for x in a]))

#Ensure whales cannot move out of range of function
def feasible(x):
    for i in range(n):
        offset = max_range[i] - min_range[i]
        x.x[i] = (x.x[i] % offset) + min_range[i]

t = 1

whales = [whale() for _ in range(N)]
whaleBest = whale()

#Function to avoid python list weirdness
def updateBest(w):
    whaleBest.x = [i for i in w.x]
    whaleBest.fit()

updateBest(min(whales, key=lambda x: x.f))
whaleBest.fit()

while t < num_cyc:

    #a on sliding scale to 0 as t -> num_cyc
    a = 2 - (2 * (t / num_cyc))

    for w in whales:
        A = [random.uniform(-1 * a, a) for _ in range(n)]
        C = [random.uniform(0, 2) for _ in range(n)]
        #l also on sliding scale
        l = random.uniform(-1, 1 - ((2 * t) / num_cyc))
        p = random.random()

        if p < 0.5:
            if size(A) < 1:
                #Encircling
                D = [abs(y[0] - y[1]) for y in zip([x[0] * x[1] for x in zip(C, whaleBest.x)], w.x)]
                w.x = [y[0] - y[1] for y in zip(whaleBest.x, [x[0] * x[1] for x in zip(A, D)])]
                feasible(w)

            else:
                #Searching
                u = random.choice([x for x in whales if x != w])
                D = [abs(y[0] - y[1]) for y in zip([x[0] * x[1] for x in zip(C, u.x)], w.x)]
                w.x = [y[0] - y[1] for y in zip(u.x, [x[0] * x[1] for x in zip(A, D)])]
                feasible(w)
        
        else:
            #Bubble net attacking
            D = [abs(x[0] - x[1]) for x in zip(whaleBest.x, w.x)]
            cos = math.exp(b * l) * math.cos(2 * math.pi * l)
            w.x = [y[0] + y[1] for y in zip(whaleBest.x, [x * cos for x in D])]
            feasible(w)
        
        w.fit()
        
    minWhale = min(whales, key=lambda x: x.f)
    if minWhale.f < whaleBest.f:
        updateBest(minWhale)
    
    t += 1

min_fMain = whaleBest.f
minimumMain = whaleBest.x

now_time = time.time()
elapsed_timeMain = round(now_time - start_time, 1)

print("\nYou have found a minimum value of {0} and a minimum point of {1}.".format(min_fMain, minimumMain))
print("Your elapsed time was {0} seconds.\n".format(elapsed_timeMain))

#################################################################################################

X = iv.mpf(['-5','5'])
Y = iv.mpf(['-5','5'])
variables = [X,Y]

def compute_f(point):
    x,y = point[0], point[1]
    return 0.5 * (x**4 - 16*x**2 + 5*x + y**4 - 16*y**2 + 5*y)
    return iv.exp(iv.sin(50*point[0])) + iv.sin(60*iv.exp(point[1])) + iv.sin(70*iv.sin(point[0])) + iv.sin(iv.sin(80*point[1])) - iv.sin(10*point[0]+10*point[1]) + 0.25*(point[0]**2 + point[1]**2)

start_time = time.time()

print("Interval Case")
print(num_cyc, N, b)

#Whale class with random initial positions
class whale:
    def __init__(self):
        lst = []
        sideOne = [random.uniform(x.a,x.b) for x in variables]
        sideTwo = [random.uniform(x.a,x.b) for x in variables]
        if sideOne[0] <= sideTwo[0]:
            lst.append(iv.mpf([sideOne[0],sideTwo[0]]))
        else:
            lst.append(iv.mpf([sideTwo[0],sideOne[0]]))
        if sideOne[1] <= sideTwo[1]:
            lst.append(iv.mpf([sideOne[1],sideTwo[1]]))
        else:
            lst.append(iv.mpf([sideTwo[1],sideOne[1]]))
        self.x = [x for x in lst]
        self.f = compute_f(self.x)
    
    def fit(self):
        self.f = compute_f(self.x)

#Size function
def size(a):
    return math.sqrt(sum([x ** 2 for x in a]))

#Ensure whales cannot move out of range of function
def feasible(x):
    for i in range(n):
        a = max(variables[i].a, min(x.x[i].a, variables[i].b))
        b = min(variables[i].b, max(x.x[i].b, variables[i].a))
        x.x[i] = iv.mpf([a, b])

t = 1

whales = [whale() for _ in range(N)]
whaleBest = whale()

#Function to avoid python list weirdness
def updateBest(w):
    whaleBest.x = [i for i in w.x]
    whaleBest.fit()

updateBest(min(whales, key=lambda x: x.f.mid))
whaleBest.fit()

def perturb(whale):
    strength = 0.1*whale.f.delta
    for i in range(n):
        wiggle = random.uniform(-strength, strength)
        a = whale.x[i].a + wiggle
        b = whale.x[i].b + wiggle
        whale.x[i] = iv.mpf([min(a, b), max(a, b)])

def generalizedHukuharaDifference(A, B):
    #For A - B
    C_a = A.a - B.a
    C_b = A.b - B.b
    return iv.mpf([min([C_a,C_b]),max([C_a,C_b])])

print("NOW FOR REAL")

#Ensure whales cannot move out of range of function
def feasible(x, variables):
    for i in range(n):
        a = max(variables[i].a, min(x.x[i].a, variables[i].b))
        b = min(variables[i].b, max(x.x[i].b, variables[i].a))
        x.x[i] = iv.mpf([a, b])
        
def originalAlgo(domain):
    t = 1
    while t < num_cyc:
        #a on sliding scale to 0 as t -> num_cyc
        a = 2 - (2 * (t / num_cyc))

        for w in whales:
            A = [random.uniform(-1 * a, a) for _ in range(n)]
            C = [random.uniform(0, 2) for _ in range(n)]
            #l also on sliding scale
            l = random.uniform(-1, 1 - ((2 * t) / num_cyc))
            p = random.random()

            if p < 0.05:
                perturb(w)
                feasible(w,domain)
            if p < 0.5:
                if size(A) < 1:
                    #Encircling
                    D = [generalizedHukuharaDifference(y[0],y[1]) for y in zip([x[0] * x[1] for x in zip(C, whaleBest.x)], w.x)]
                    w.x = [generalizedHukuharaDifference(y[0],y[1]) for y in zip(whaleBest.x, [x[0] * x[1] for x in zip(A, D)])]
                    feasible(w,domain)

                else:
                    #Searching
                    u = random.choice([x for x in whales if x != w])
                    D = [generalizedHukuharaDifference(y[0],y[1]) for y in zip([x[0] * x[1] for x in zip(C, u.x)], w.x)]
                    w.x = [generalizedHukuharaDifference(y[0],y[1]) for y in zip(u.x, [x[0] * x[1] for x in zip(A, D)])]
                    feasible(w,domain)
            
            else:
                #Bubble net attacking
                D = [generalizedHukuharaDifference(x[0],x[1]) for x in zip(whaleBest.x, w.x)]
                cos = iv.exp(b * l) * iv.cos(2 * iv.pi * l)
                w.x = [y[0] + y[1] for y in zip(whaleBest.x, [x * cos for x in D])]
                feasible(w,domain)
                
            w.fit()
            
        minWhale = min(whales, key=lambda x: x.f.mid + 0.1*x.f.delta)
        if minWhale.f.mid < whaleBest.f.mid:
            updateBest(minWhale)
        
        t += 1
    min_f = whaleBest.f
    minimum = whaleBest.x
    return min_f, minimum

def repeat(domain,eps=1e-14):
    fstar = iv.mpf(['-10','10'])
    count = 0
    totalTime = 0
    while fstar.delta > eps:
        start = time.time()
        fstar,domain = originalAlgo(domain)
        end = time.time()
        timer = end - start
        print(domain, fstar, timer)
        count += 1
        totalTime += timer
    return eps, count-1, domain, fstar, totalTime

print(repeat(variables))
