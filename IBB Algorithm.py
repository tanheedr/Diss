import time
from mpmath import iv

iv.dps = 60
iv.pretty = False

def union(x,y):
    if x.a <= y.a:
        lower = x.a
    else:
        lower = y.a
    if x.b >= y.b:
        upper = x.b
    else:
        upper = y.b
    return iv.mpf([lower,upper])

def intersection(x,y):
    if x.a > y.b or x.b < y.a:
        return None
    if x.a <= y.a:
        lower = y.a
    else:
        lower = x.a
    if x.b >= y.b:
        upper = y.b
    else:
        upper = x.b
    return iv.mpf([lower,upper])

def extendedIntersection(A,B):
    if type(B) == list:
        arr = []
        for i in range(0,len(B)):
            arr.append(intersection(A,B[i]))
        return arr
    else:
        return intersection(A,B)

def generalizedHukuharaDifference(A, B):
    #For A - B
    C_a = A.a - B.a
    C_b = A.b - B.b
    return iv.mpf([min([C_a,C_b]),max([C_a,C_b])])

def extendedGeneralizedHukuharaDifference(A, B):
    #For A - B
    arr = []
    if type(B) != list and type(A) != list:
        return generalizedHukuharaDifference(A,B)
    elif type(B) == list and type(A) != list:
        for i in range(0, len(B)):
            C_a = A.a - B[i].a
            C_b = A.b - B[i].b
            arr.append(iv.mpf([min([C_a,C_b]),max([C_a,C_b])]))
        return arr
    elif type(A) == list and type(B) != list:
        for i in range(0, len(A)):
            C_a = A[i].a - B.a
            C_b = A[i].b - B.b
            arr.append(iv.mpf([min(C_a,C_b),max(C_a,C_b)]))
        return arr

def partialDeriv_x(f, *x):
    if len(x) == 2:
        z = x[0]
    else:
        z = x[2]
    y = x[1]
    x = x[0]
    #Centered form
    h = 0.1*x.delta if x.delta > 0 else 10*z.delta
    return generalizedHukuharaDifference(f(x+h,y),f(x-h,y)) / (2*h)

def partialDeriv_y(f, *x):
    if len(x) == 2:
        z = x[0]
    else:
        z = x[2]
    y = x[1]
    x = x[0]
    #Centered form
    h = 0.1*y.delta if y.delta > 0 else 10*z.delta
    return generalizedHukuharaDifference(f(x,y+h),f(x,y-h)) / (2*h)

def secondDerivative_3d(f, *x):
    if len(x) == 2:
        z = x[0]
    else:
        z = x[2]
    y = x[1]
    x = x[0]
    #Centered form
    h_x = 0.1*x.delta if x.delta > 0 else 10*z.delta
    h_y = 0.1*y.delta if y.delta > 0 else 10*z.delta
    f_xx = generalizedHukuharaDifference(partialDeriv_x(f,x+h_x/2,y),partialDeriv_x(f,x-h_x/2,y)) / (h_x)
    f_xy = generalizedHukuharaDifference(partialDeriv_y(f,x+h_x/2,y),partialDeriv_y(f,x-h_x/2,y)) / (h_x)
    f_yy = generalizedHukuharaDifference(partialDeriv_y(f,x,y+h_y/2),partialDeriv_y(f,x,y-h_y/2)) / (h_y)
    return f_xx, f_xy, f_yy

def f_objective(x_1, x_2):
    return iv.exp(iv.sin(50*x_1)) + iv.sin(60*iv.exp(x_2)) + iv.sin(70*iv.sin(x_1)) + iv.sin(iv.sin(80*x_2)) - iv.sin(10*(x_1+x_2)) + iv.mpf(['0.25','0.25'])*(x_1**2+x_2**2)

def mainPartial3d_x(f,x,y):
    h = 0.001*x.delta
    forward = f(x+h/2,y)
    backward = f(x-h/2,y)
    return (forward - backward) / h

def mainPartial3d_y(f,x,y):
    h = 0.001*y.delta
    forward = f(x,y+h/2)
    backward = f(x,y-h/2)
    return (forward - backward) / h

def mainDeriv3d_x(f,x,y):
    if x.delta > 0:
        return mainPartial3d_x(f,x,y)
    else:
        h = x/100000
        return (f(x+h,y) - f(x-h,y)) / (2*h)

def mainDeriv3d_y(f,x,y):
    if y.delta > 0:
        return mainPartial3d_y(f,x,y)
    else:
        h = y/100000
        return (f(x,y+h) - f(x,y-h)) / (2*h)

def mainSecondDerivative3d(f,x,y):
    h_x = 0.1*x.delta
    h_y = 0.1*x.delta
    f_xx = generalizedHukuharaDifference(mainDeriv3d_x(f,x+h_x/2,y),mainDeriv3d_x(f,x-h_x/2,y)) / (h_x)
    f_xy = generalizedHukuharaDifference(mainDeriv3d_y(f,x+h_x/2,y),mainDeriv3d_y(f,x-h_x/2,y)) / (h_x)
    f_yy = generalizedHukuharaDifference(mainDeriv3d_y(f,x,y+h_y/2),mainDeriv3d_y(f,x,y-h_y/2)) / (h_y)
    return f_xx, f_xy, f_yy

def hessianMatrix(f, *x):
    if len(x) == 2:
        z = x[0]
    else:
        z = x[2]
    y = x[1]
    x = x[0]
    secondDeriv = mainSecondDerivative3d(f, x, y)
    hessian = [[secondDeriv[0]],
               [secondDeriv[1]],
               [secondDeriv[1]],
               [secondDeriv[2]]]
    return hessian

def inverseMatrix(A):
    det = A[0][0] * A[3][0] - A[1][0] * A[2][0]
    A[0], A[3] = A[3], A[0]
    A[1][0] = -A[1][0]
    A[2][0] = -A[2][0]
    for i in range(0,4):
        A[i][0] = (1/det)*A[i][0]
    return A

def mainNewtonStep3d(f,x,y):
    hessian = hessianMatrix(f, x, y)
    hessianInv = inverseMatrix(hessian)
    c_x = x.mid
    c_y = y.mid
    hessMult_x = hessianInv[0][0] * mainDeriv3d_x(f, x, y) + \
                 hessianInv[1][0] * mainDeriv3d_y(f, x, y)
    hessMult_y = hessianInv[2][0] * mainDeriv3d_x(f, x, y) + \
                 hessianInv[3][0] * mainDeriv3d_y(f, x, y)
    N_x = extendedGeneralizedHukuharaDifference(c_x, hessMult_x)
    N_y = extendedGeneralizedHukuharaDifference(c_y, hessMult_y)
    L_x = extendedIntersection(x,N_x)
    L_y = extendedIntersection(y,N_y)
    
    temp = []
    if type(L_x) == list:
        for item in L_x:
            if item is not None:
                temp.append(item)
        L_x = temp
        if len(L_x) > 1:
            if L_x[0].a > L_x[1].a:
                temp = L_x[0]
                L_x[0] = L_x[1]
                L_x[1] = temp
    elif L_x is not None:
        L_x = [L_x]

    temp = []
    if type(L_y) == list:
        for item in L_y:
            if item is not None:
                temp.append(item)
        L_y = temp
        if len(L_y) > 1:
            if L_y[0].a > L_y[1].a:
                temp = L_y[0]
                L_y[0] = L_y[1]
                L_y[1] = temp
    elif L_y is not None:
        L_y = [L_y]

    return L_x, L_y

def mainMonotonicityTest3d(f,x,y):
    f_x = mainDeriv3d_x(f,x,y)
    f_y = mainDeriv3d_y(f,x,y)
    if 0 in f_x and 0 in f_y:
        return False
    else:
        return True

def mainConcavityTest3d(f,x,y):
    s = secondDerivative_3d(f,x,y)
    det = s[0]*s[2] - s[1]*s[1]
    if s[0].b < 0 and det.b < 0:
        return True
    else:
        return False

def mainMeanValueForm3d(f,x,y,c_x,c_y):
    return extendedIntersection(f(c_x,c_y) + mainDeriv3d_x(f,x,y)*(x-c_x) + mainDeriv3d_y(f,x,y)*(y-c_y), f(x,y))

def sort(L):
    return sorted(L, key = lambda x:x[2])

def mainMidptTest(L,ftilde):
    temp = []
    for item in L:
        if item[2].a <= ftilde.a:
            temp.append(item)
    return temp

def printL(L):
    for item in L:
        print(item)

#Refinements for better approx
def calculate(n,f,x,y):
    intervList = []
    intervReturned = []
    for i in range(1,n+1):
        xSubdiv = iv.mpf([x.a + (i-1)*(x.delta/n),x.a + i*(x.delta/n)])
        for j in range(1,n+1):
            ySubdiv = iv.mpf([y.a + (j-1)*(y.delta/n),y.a + j*(y.delta/n)])
            intervList.append([xSubdiv, ySubdiv])
    for item in intervList:
        intervReturned.append(f(item[0],item[1]))
    while len(intervReturned) != 1:
        z = union(intervReturned[0], intervReturned[1])
        intervReturned.pop(0)
        intervReturned.pop(0)
        intervReturned.insert(0,z)
    return intervReturned

def mainIBB3d(f,x,y,n=30,eps=1e-4):
    s = time.time()
    c_x = x.mid
    c_y = y.mid
    f_c = f(c_x,c_y)
    ftilde = f_c
    L = []
    R = []
    A, B, C, D = f(x.a,y.a), f(x.a,y.b), f(x.b,y.a), f(x.b,y.b)
    ftilde = min(ftilde, A, B, C, D)
    if ftilde.a >= A.a:
        L.append([x.a, y.a, A])
    if ftilde.a >= B.a:
        L.append([x.a, y.b, B])
    if ftilde.a >= C.a:
        L.append([x.b, y.a, C])
    if ftilde.a >= D.a:
        L.append([x.b, y.b, D])
    flag = False
    k = 0
    epsForN = 1e-2
    while flag==False:
        k+=1
        U = [[iv.mpf([x.a,c_x]), iv.mpf([y.a,c_y])],
             [iv.mpf([x.a,c_x]), iv.mpf([c_y,y.b])],
             [iv.mpf([c_x,x.b]), iv.mpf([y.a,c_y])],
             [iv.mpf([c_x,x.b]), iv.mpf([c_y,y.b])]]
        for i in range(0,4):
            if mainMonotonicityTest3d(f,U[i][0],U[i][1]) == True:
                continue
            f_u = calculate(n,f, U[i][0], U[i][1])[0]
            if ftilde < f_u.a:
                continue
            if mainConcavityTest3d(f,U[i][0],U[i][1]) == True:
                continue
            V = mainNewtonStep3d(f, U[i][0], U[i][1])
            if V is None:
                continue
            for j in range(0,len(V)//2):
                if mainMonotonicityTest3d(f, V[2*j][0], V[2*j+1][0]) == True:
                    continue
                c_vx = V[2*j][0].mid
                c_vy = V[2*j+1][0].mid
                f_v = calculate(n,f,V[2*j][0],V[2*j+1][0])[0]
                if ftilde >= f_v.a:
                    L.append([V[2*j][0],V[2*j+1][0], f_v.a])
                    L = sort(L)
        flag2 = False
        while len(L) != 0 and flag2 == False:
            x, y, f_x = L[0][0], L[0][1], L[0][2]
            L.remove(L[0])
            c_x = x.mid
            c_y = y.mid
            ftilde = min(ftilde, f(c_x,c_y))
            L = mainMidptTest(L,ftilde)
            fstar = iv.mpf([f_x,ftilde])
            if fstar.delta < eps and x.delta < eps and y.delta < eps:
                R.append([x, y, fstar])
                flag = True
                break
            else:
                flag2 = True
    return n, eps, k, time.time()-s, R

X = iv.mpf(['-10','10'])
Y = iv.mpf(['-10','10'])
start = time.time()
print("OBJECTIVE")
print(mainIBB3d(f_objective,X,Y))
end = time.time()
print(end-start)

listOfDifferentN = [mainIBB3d(f_objective,X,Y,5),
                    mainIBB3d(f_objective,X,Y,10),
                    mainIBB3d(f_objective,X,Y,20),
                    mainIBB3d(f_objective,X,Y,30),
                    mainIBB3d(f_objective,X,Y,40),
                    mainIBB3d(f_objective,X,Y,50),
                    mainIBB3d(f_objective,X,Y,75),
                    mainIBB3d(f_objective,X,Y,100)]
printL(listOfDifferentN)


##TEST FUNCS##

def fAckley(x,y):
    #Min: f(0,0) = 0
    return -20*iv.exp(-0.2*iv.sqrt(0.5*(x**2+y**2))) - iv.exp(0.5*(iv.cos(2*iv.pi*x)+iv.cos(2*iv.pi*x))) + iv.exp(1) + 20

X = iv.mpf(['-5','5'])
Y = iv.mpf(['-5','5'])
start = time.time()
print("ACKLEY")
print(mainIBB3d(fAckley,X,Y))
end = time.time()
print(end-start)

def fBeale(x,y):
    #Min: f(3,0.5) = 0
    return (1.5-x+x*y)**2 + (2.25-x+x*y**2)**2 + (2.625-x+x*y**3)**2

X = iv.mpf(['-4.5','4.5'])
Y = iv.mpf(['-4.5','4.5'])
start = time.time()
print("\nBEALE")
print(mainIBB3d(fBeale,X,Y))
end = time.time()
print(end-start)

def fMcCormick(x,y):
    #Min: f(-0.54719,-1.54719)=1.9133
    return iv.sin(x+y) + (x-y)**2 - 1.5*x + 2.5*y + 1

X = iv.mpf(['-1.5','4'])
Y = iv.mpf(['-3','4'])
start = time.time()
print("\nMCCORMICK")
print(mainIBB3d(fMcCormick,X,Y))
end = time.time()
print(end-start)

def fStyblinskiTang(x,y):
    #Min: -78.33234 < f(-2.903534,-2.903534) < -78.33232
    return 0.5 * (x**4 - 16*x**2 + 5*x + y**4 - 16*y**2 + 5*y)

X = iv.mpf(['-5','5'])
Y = iv.mpf(['-5','5'])
start = time.time()
print("\nSTYBLINSKI-TANG")
print(mainIBB3d(fStyblinskiTang,X,Y))
end = time.time()
print(end-start)
