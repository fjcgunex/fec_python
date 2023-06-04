import math
import numpy as np
import random
import bisect
import time
from operator import itemgetter, attrgetter
from os import system, remove
import pandas as pd
#-----
def propaga(dt):
    global vx, vy, x, y
    x = x + vx * dt
    y = y + vy * dt
#-----
def midedist2(i, j):
    global x, y
    dx = x[i] - x[j]
    dy = y[i] - y[j]
    dist2 = (np.power(dx, 2) + np.power(dy, 2)) - 4 * np.power(r, 2)
#-----
def tcol(i, j):
    global vx, vy, x, y
    dx = x[i] - x[j]
    dy = y[i] - y[j]
    dvx = vx[i] - vx[j]
    dvy = vy[i] - vy[j]
    drdv = dx * dvx + dy * dvy
    if(drdv > 0):
        vct = float('inf')
    else:
        dist2 = (np.power(dx, 2) + np.power(dy, 2)) - 4 * np.power(r, 2)
        raiz = np.power(drdv, 2) - dist2 * (np.power(dvx, 2) + np.power(dvy, 2))
        if(raiz < 0):
            vct = float('inf')
        else:
            vdt = dist2 / (np.sqrt(raiz) - drdv)
            xicol = x[i] + vx[i] * vdt
            yicol = y[i] + vy[i] * vdt
            xjcol = x[j] + vx[j] * vdt
            yjcol = y[j] + vy[j] * vdt
            if(np.abs(xicol) > lxr):
                vdt = float('inf')
            elif(np.abs(xjcol) > lxr):
                vdt = float('inf')
            elif(np.abs(yicol) > lyr):
                vdt = float('inf')
            elif(np.abs(yjcol) > lyr):
                vdt = float('inf')
            else:
                bisect.insort(listacol, [vdt, [i, j]])
#-----
def tpcol(i):
    global vx, vy, x, y
    if(vx[i] == 0):
        tx = float('inf')
    elif(vx[i] < 0):
        ltx = [-(lxr + x[i]) / vx[i], -1] # pared izq.
    elif(vx[i] > 0):
        ltx = [(lxr - x[i]) / vx[i], -3] # pared dcha.
    if(vy[i] == 0):
        ty = float('inf')
    elif(vy[i] < 0):
        lty = [-(lyr + y[i]) / vy[i], -2] # pared inf.
    elif(vy[i] > 0):
        lty = [(lyr - y[i]) / vy[i], -4] # pared sup.
    ltm = sorted([ltx, lty], key = itemgetter(0))
    vdt = ltm[0][0]
    im = ltm[0][1]
    bisect.insort(listacol, [vdt, [i, im]])
#-----
def pcolisiona(ii):
    global vx, vy, x, y
    if((ii[1] == -1) or (ii[1] == -3)):
        vx[ii[0]] = -vx[ii[0]]
    elif((ii[1] == -2) or (ii[1] == -4)):
        vy[ii[0]] = -vy[ii[0]]
    phi = np.sqrt(np.pi / 1000) * np.random.randn()
    vx[ii] = vx[ii] * np.cos(phi) - vy[ii] * np.sin(phi)
    vy[ii] = vx[ii] * np.sin(phi) + vy[ii] * np.cos(phi)
#-----
def colisiona(par):
    global vx, vy, x, y
    i = par[0]
    j = par[1]
    dx = x[i] - x[j]
    dy = y[i] - y[j]
    sigma_norma = np.sqrt(np.power(dx, 2) + np.power(dy, 2))
    sigmax = dx / sigma_norma
    sigmay = dy / sigma_norma
    gsigma = (vx[i] - vx[j]) * sigmax + (vy[i] - vy[j]) * sigmay
    vx[i] = vx[i] - 0.5 * (1 + alfa) * gsigma * sigmax
    vy[i] = vy[i] - 0.5 * (1 + alfa) * gsigma * sigmay
    vx[j] = vx[j] + 0.5 * (1 + alfa) * gsigma * sigmax
    vy[j] = vy[j] + 0.5 * (1 + alfa) * gsigma * sigmay
#-----
def write_micr_state(ja):
    global vx, vy, x, y
    print ("####### it: ########", it)
    print ("####### no. archivo: ########", ja)
    inum='{0:04d}'.format(ja)
    nombre='./uex/fisica_estadistica/ficheros_python/xy'+inum+'.txt'
    with open(nombre, 'w') as archivo:
        for i in range(npart):
            archivo.write('{0:10.2f} {1:10.2f}\n'.format(x[i], y[i]))
    archivo.closed
    nombre = './uex/fisica_estadistica/ficheros_python/vxvy'+inum+'.txt'
    with open(nombre, 'w') as archivo:
        for i in range(npart):
            archivo.write('{0:10.2f} {1:10.2f}\n'.format(vx[i], vy[i]))
    archivo.closed
#-----
def initialize_random():
    global vx, vy, x, y
    x[0] = random.uniform(-lxr, lxr)
    y[0] = random.uniform(-lyr, lyr)
    for i in range(1, npart):
        dr = False
        while (dr == False):
            x[i] = random.uniform(-lxr, lxr)
            y[i] = random.uniform(-lyr, lyr)
            for j in range(0, i):
                dr = (np.power((x[i] - x[j]), 2) + np.power((y[i] - y[j]), 2)) > (4 * np.power(r, 2))
                if(dr == False):
                    break
    for i in range(npart):
        vx[i] = np.random.randn() * np.power(npart, -0.5)
        vy[i] = np.random.randn() * np.power(npart, -0.5)
#-----
def calculate_averages(ja):
    global temp, a2, vx, vy, x, y
    temp[ja] = 0.
    a2[ja] = 0.
    for i in range(npart):
        vv = vx[i] * vx[i] + vy[i] * vy[i]
        temp[ja] = temp[ja] + vv
        a2[ja] = a2[ja] + np.power(vv, 2)
    temp[ja] = temp[ja] / npart
    a2[ja] = a2[ja] / (temp[ja] * temp[ja] * npart)
    a2[ja] = (a2[ja] - 2.0) * 0.5
#-----
def write_averages_evol():
    nombre='./uex/fisica_estadistica/ficheros_python/temp.txt'
    xy= pd.DataFrame(np.array([[temp[i], a2[i]] for i in range(len(temp))]))
    xy.to_csv(nombre, sep='\t', header = ['t', 'a2'] , index = False, float_format = '%8.5f')
#-----
# programa principal
LXMIN = 10
LXMAX = 100000
LYMIN = 10
LYMAX = 100000
NU = 0.5
NSIM = 100
#-----
while(True):
    lx = float(input("Longitud de la región en el eje OX: "))
    if((lx < LXMIN) or (lx > LXMAX)):
        print(f"Longitud fuera de rago, debe estar en [{LXMIN},{LXMAX}]")
    else:
        break
#-----
while(True):
    ly = float(input("Longitud de la región en el eje OY: "))
    if((ly < LYMIN) or (ly > LYMAX)):
        print(f"Longitud fuera de rago, debe estar en [{LYMIN},{LYMAX}]")
    else:
        break
#----
print("Para la lectura del radio,\n debe introducirse un valor tal que quepa al menos un disco\n en el menor cuadrado contenido en el rectángulo.")
while(True):
    r = float(input("Radio: "))
    if((r < 0) or (r > 0.5 * min(lx,ly))):
        print(f"Radio fuera de rango, debe estar en [{0},{0.5 * min(lx, ly)}]")
    else:
        break
#
lxr = lx * 0.5 - r
lyr = lx * 0.5 - r
#-----
while(True):
    npart = int(input("Número de partículas: "))
    areadiscos = npart * np.pi * np.power(r, 2)
    if((areadiscos >= NU * lx * ly) or (npart <= 0)):
        print(f"Intorduce otra cantidad, debe estar entre 0 y {np.floor(NU * lx * ly / (np.pi * np.power(r, 2)))}")
    else:
        break
#-----
nt = 100 * npart
#print(f"El número de pasos a dar será {nt}")
alfa = 1.0
tol = 1.0e-20
ncp = 1.0 * nt / npart
#-----
utermo = 1
#-----
vx = np.array([0. for i in range(npart)])
vy = np.array([0. for i in range(npart)])
#-----
x = np.array([0. for i in range(npart)])
y = np.array([0. for i in range(npart)])
#-----
temp = np.array([0. for i in range(nt + 1)])
a2 = np.array([0. for i in range(nt + 1)])
#-----
listacol = []
listacol_orden = []
ij = []
#-----
t = 0. # inicializa el tiempo
dt = 0.
it = 0
random.seed()
inicio = time.time()
#print("simulacion md, alfa = ", alfa)
#print("cols/part (total): ", ncp)
#print("its. entre snapshots: ", utermo)
#print("no. de archivos: ", nt / utermo)
initialize_random()
for i in range(npart - 1):
    for j in range(i + 1, npart):
        tcol(i, j)
#-----
for i in range(npart):
    tpcol(i)
#-----
Matvx = np.zeros((npart, nt +1))
Matvy = np.zeros((npart, nt +1))
Matvx[:,0] = vx
Matvy[:,0] = vy
Matvy = np.zeros((npart, nt +1))
for it in range(1, nt + 1):
    dt = listacol[0][0] * (1 - tol)
    ij = listacol[0][1]
    listacol = list(filter(lambda x : x[1][0] != ij[0], listacol))
    listacol = list(filter(lambda x : x[1][1] != ij[0], listacol))
    if(ij[1] > 0):
        listalcol = list(filter(lambda x : x[1][0] != ij[1], listacol))
        listacol = list(filter(lambda x : x[1][1] != ij[1], listacol))
    t = t + dt
    limit = range(len(listacol))
    c=[[listacol[i][0] - dt, listacol[i][1]] for i in limit]
    listacol = c
    propaga(dt)
    if(ij[1] < 0):
        pcolisiona(ij)
    else:
        colisiona(ij)
    i = ij[0]
    tpcol(i)
    for j in range(i):
        tcol(j, i)
    for j in range(i + 1, npart):
        tcol(i, j)
    if(ij[1] > 0):
        i = ij[1]
        tpcol(i)
        for j in range(i):
            tcol(j, i)
        for j in range(i + 1, npart):
            tcol(i, j)
    Matvx[:,it] = vx
    Matvy[:,it] = vy
while(True):
    j = int(input("Instante de tiempo para calcular la autocorrelación: "))
    if(lx > nt):
        print(f"Tiempo fuera de rago, debe estar en [{0},{nt}]")
    else:
        break
print(Matvx)
coefx = np.corrcoef(Matvx[:,0],Matvx[:,j])
coefx = coefx[0,1]
coefy = np.corrcoef(Matvx[:,0],Matvx[:,j])
coefy = coefy[0,1]
print(f"La correlación en el eje OX es {coefx}")
print(f"La correlación en el eje OY es {coefy}")