# -*- coding: utf-8 -*-
'''
Autor: Adrián Robles Arques
Física computacional

PROGRAMA PARA LA SIMULACIÓN DE REDES CRISTALINAS

FUNCIONES
-------------------------------------------------------------------------------
#1: RedCristalina(X,Y,Z,T,par, pred = 3.603, tred = 'cubico', dist = 'primitiva')
    
    Función encargada de generar una red de átomos de una dimensión determinada. La creación
    se realiza de forma iterativa a partir de una celda tipo básica que se repite en el espacio.
    Esto provoca que si se representa tal como aparece, falto un átomo en la arista opuesta a
    la inicial (0,0,0).
    
    Entradas:
        X,Y,Z: (int) Tamaño de la red, número de celdas a crear en cada una
                de las direcciones espaciales
                
        T: (float) Temperatura base de la red, para ajustar la distribución de velocidades.
        par: (tupla) Conjunto de constantes necesarias, masa de las partículas
                y constante de Boltzmann
                
        pred: (float, opt) parámetro de red, separación en amstrong entre los elementos de la red.
                predefinido para los datos del cobre, pred = 3.603
                
        tred: (string, opt) Tipo geométrico de la red
                            - 'cubrico': geometría cúbica
                            - Pendiente de implementar otras geometrías
                            
        dis: (string, opt) Distribución de los átomos dentro de la red:
                            - 'primitiva': Un átomo en cada vértice intersección
                            - 'BCC': Añadimos un átomo en el centro cada sección
                            - 'FCC': Añadimos un átomo en el centro de cada cara de la sección
                            - 'diamond': Distribución tipo de diamante, atomos apareados
    Salida:
        CoordX, CoordY, CoordZ, Velocidades: (Numpy array) posiciones de los átomos
            de la red, en amstrong, situando una arista de la red en la posición (0,0,0), así
            como un array con las velocidades en forma de vector de cada una de las partículas.
            
#2: Calcpot(Coord, tamaño, irc=1, sigma = 2.3151, epsilon = 0.167, pred = 3.603, mode = 'finito'):
    
    Función encargada del cálculo de las energías potenciales y de las fuerzas que actuan
    sobre cada una de las partículas. Calculada a partir de la divergencia del potencial.
    
    Entradas:
        Coord: (list(3,N)) Conjunto de arrays de las coordenadas X, Y y Z de la red, por separado.
            en Amstrongs.
            
        tamaño: (tupla(3,)) Conjunto de medida de la red, numero de celdas individuales en cada eje.
            p.ej: Si se crea una red con RedCristalina(2,2,2...) en este caso tamaño = (2,2,2)
            
        irc: (int) Tipo de aproximación para el cálculo de potencial, predefinida como 1 
            para cálculo de potencial de primeros vecinos, 2 para contar hasta segundo,
            3 para terceros, etc.
            
        sigma: (float, opt) Conductividad eléctrica característica del átomo en cuestión.
            inicializado con valores para el cobre. Usado para el cálculo de potencial.
            
        epsilon: (float, opt) Permitividad eléctrica de los átomos de la red. Inicializada
            para valores del cobre.
            
        pred: (float, opt) parámetro de red, separación en amstrong entre los elementos de la red.
                predefinido para los datos del cobre, pred = 3.603
                
        mode: (string, opt) Modo de cálculo de potenciales:
            
                - 'finito': Cálculo de potenciales para la red aislada, finita.
                            implica que el potencial para cada átomo cambiará en función
                            de su ubicación en la red y del tamaño de la misma.
                            
                - 'infinito': Emplea condiciones de contorno para, a pesar de tener
                            una red de tamaño finito, suponer que el cálculo de potencial
                            de cada átomo se realiza en una red igual pero que se extiende
                            indefinidamente. Todos los átomos deberían tener igual potencial
                            al estar rodeados por la misma cantidad de átomos en todas direcciones.
    Salidas: 
        Enerind: (list(N,)) Energía potencial individual de cada una de las partículas.
        
        Etotal: (float) Energía potencial total de la red.
        
        Fuerzas: (list(N,3)) Vector [Fx,Fy,Fz] de fuerzas actuando sobre cada partícula.
    
#3: Mov3D(Coord, Vel, Fvec, ini, time, h, pred = 3.603, modo = 'infinito'):
    
    Función encargada de resolver las ecuaciones diferenciales necesarias para representar
    el movimiento de las partículas de la red cristalina partiendo de las posiciones y
    velocidades iniciales de la misma, así como de las fuerzas originadas por la propia red
    debido a la interacción electromagnética entre átomos.
    
    Esta función requiere a su vez de la función Calcpot definida anteriormente, realiza
    llamadas a la misma para calculo de fuerzas a medida que evoluciona el tiempo.
    
    Entradas:
        
        Coord: (list(N)) Coordenadas iniciales de los átomos que conforman la red, en Amstrong.
        
        Vel: (list(3,N)) Velocidades iniciales de cada átomo separados por cada dirección.
        
        Fvec: (list(N,3)) Velocidades vectoriales iniciales para cada átomo. Coincide con el 
            formato de salida de  la función Calcpot para las Fuerzas (Calcpot[2]).
            
        ini: (tupla(3,)) Tupla que continene datos necesarios para operaciones internas y
            llamadas a la función Calcpot. Ini = (masa(eV), irc, tamaño) irc y tamaño como
            la entrada de Calcpot.
            
        time: (tupla(2,)) Tupla que contiene el tiempo inicial y final deseado. time= (t0,tf)
        
        h: (int) diferencia de tiempo entre cada evaluación. El número de instantes evaluados
            es igual a (tf - t0)/h
            
        pred (float, opt) Distancia entre átomos de la red en Amstrong. Inicializado para 
            valores del cobre.
            
        modo: (string, opt) Argumento opcional a pasar a la función de cálculo de potencial.
            inicializado con el modo continuo con condiciones de frontera definidas.
            
    Salidas:
        
        T: (numpy array) Array de tiempos usados para la simulación
        
        Coord: (tupla (3,T,N)) Conjunto de coordenadas divididas por cada dirección
            por cada instante de tiempo simulado para todas las partículas de la red.
            
        Vel: (tupla (3,T,N)) Conjunto de velocidades divididas por cada dirección
            por cada instante de tiempo simulado para todas las partículas de la red.
            
        Epot: (list (T,)) Lista de la energía potencial total de la red para cada instante
            de tiempo simulado.
            
        FUE: (list (T,N,3)) Lista de fuerzas en forma vectorial por cada instante de tiempo
            para cada una de las partículas de la red.
        
#4  Temperatua(velocidades,K,mode = 'vector'):
    
    Función encargada de calcular la temperatura y energía cinética de una red cristalina
    partiendo de la velocidad conocida de cada una de las partículas que la componen.
    
    Entradas:
        
        Velocidades: (lista (3,N) o lista(N,3)) Lista de las velocidades para cada
            partícula de la red. Admite dos modos de presentación de velocidades, según
            convenga. Se debe indicar el modo de presentación mediante la entreda mode.
        
        K: (float) Constante de Boltzmann en las unidades requeridas
        
        mode: (string, opt) Entrada que indica la forma de presentación de las velocidades.
            
                - 'vector': Indica que se pasa una lista(N,3) donde cada elemento de la lista es el vector
                        velocidad de cada partícula de la red. 
                        
                - 'lista': Indica que se pasa una lista (3,N) que debe descomponerse en tres listas
                        donde cada una representa las componentes de la velocidad en cada dirección
                        para cada una de las partículas de la red.
                        
    Salidas:
        
        T: (float) Temperatura de la red, en Kelvin
        
        Ecin: (float) Energía cinética de la red, unidades en función de K y velocidades.
-------------------------------------------------------------------------------    
'''
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as npran
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

#Vamos con la red, de tamaño X, Y, Z, pudiendo definir los parámetros y la distribución.
def RedCristalina(X,Y,Z,T,par,pred = 3.603 , tred = "cubico", dist = 'primitiva'):
  
    m,K = par
    
    if tred == 'cubico':
        Coordx = [] #Separamos las coordenadas
        Coordy = []
        Coordz = []
        
        velocidades = []

        if dist == 'primitiva':
          for i in range(X):
            for j in range(Y):
              for k in range(Z):
                Coordx.append(i)
                Coordy.append(j)
                Coordz.append(k)
                
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
    
        elif dist == 'BCC': #Añadimos un átomo central
          for i in range(X):
            for j in range(Y):
              for k in range(Z):
                Coordx.append(i)
                Coordy.append(j)
                Coordz.append(k)

                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
    
                Coordx.append(i+0.5)
                Coordy.append(j+0.5)
                Coordz.append(k+0.5) 
    
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
        elif dist == 'FCC': #Ahora cada cara tiene que tener un átomo en el centro
          for i in range(X):
            for j in range(Y):
              for k in range(Z):
                Coordx.append(i)
                Coordy.append(j)
                Coordz.append(k)
    
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
                Coordx.append(i+0.5)
                Coordy.append(j+0.5)
                Coordz.append(k)
                
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
                Coordx.append(i+0.5)
                Coordy.append(j)
                Coordz.append(k+0.5)
    
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
                Coordx.append(i)
                Coordy.append(j+0.5)
                Coordz.append(k+0.5)
                
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
    
        elif dist == 'diamond': #Este es como el anterior, pero con pares de átomos
          Coordx1 = []
          Coordy1 = []
          Coordz1 = []
        
          Coordx2 = []
          Coordy2 = []
          Coordz2 = []
          
          for i in range(X):
            for j in range(Y):
              for k in range(Z):
                Coordx1.append(i)
                Coordy1.append(j)
                Coordz1.append(k)
                
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
    
                Coordx1.append(i+0.5)
                Coordy1.append(j+0.5)
                Coordz1.append(k)
    
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
                Coordx1.append(i+0.5)
                Coordy1.append(j)
                Coordz1.append(k+0.5)
    
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
                Coordx1.append(i)
                Coordy1.append(j+0.5)
                Coordz1.append(k+0.5)
    
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
                Coordx2.append(i+0.25)
                Coordy2.append(j+0.25)
                Coordz2.append(k+0.25)
    
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
                Coordx2.append(i+0.75)
                Coordy2.append(j+0.75)
                Coordz2.append(k+0.25)
    
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
                Coordx2.append(i+0.75)
                Coordy2.append(j+0.25)
                Coordz2.append(k+0.75)
    
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
                Coordx2.append(i+0.25)
                Coordy2.append(j+0.75)
                Coordz2.append(k+0.75)
    
                vel = npran.random_sample((3,))-0.5
                velocidades.append(vel)
                
            Coordx = [Coordx1,Coordx2]
            Coordy = [Coordy1,Coordy2]
            Coordz = [Coordz1,Coordz2]
        
    ECINrand = 0
    N = len(velocidades)
    for i in velocidades:
        ECINrand += m/2 * (i[0]**2+i[1]**2+i[2]**2)
        
    Trand = 2/(3*K*N) * ECINrand
    facescala = np.sqrt(T/Trand)
    
    for i in range(N):
        velocidades[i] = facescala*velocidades[i]    
    
    return pred*np.array(Coordx), pred*np.array(Coordy), pred*np.array(Coordz), velocidades

#------------------------------------------------------------------------------
def Calcpot(Coord, tamaño, irc=1, sigma = 2.3151, epsilon = 0.167, pred = 3.603, mode = 'finito'): #Función para calcular el potencial
    N = len(Coord[0])
    X, Y, Z = Coord
    LX, LY, LZ = pred*np.array(tamaño)
    rc = irc*sigma

    enerind = [] #Energía individual de cada partícula
    Find = [] #Módulo de la fuerza por cada partícula
    FX, FY, FZ = ([],[],[]) #Fuerza como vector
    
    total = 0

    if mode == 'infinito':
      for i in range(N):
        e = 0
        Fx, Fy, Fz = (0,0,0)
        for j in range(N):
            if i != j: #Aplicamos las condiciones para evitar errores y las con. de contorno
              dx = X[i] - X[j]
              if dx > LX/2:
                  dx -= LX
              elif dx < - LX/2:
                  dx += LX
              dy = Y[i] - Y[j]
              if dy > LY/2:
                  dy -= LY
              elif dy < - LY/2:
                  dy += LY
              dz = Z[i] - Z[j]
              if dz > LZ/2:
                  dz -= LZ
              elif dz < - LZ/2:
                  dz += LZ
        
              dist = np.sqrt(dx**2 + dy**2 + dz**2)
        
              if dist < rc:

                V = 4*epsilon*((sigma/dist)**12-(sigma/dist)**6)
                total += V/2
                e += V
                
                dVdr = 4*epsilon*(-12*(sigma**12/dist**13) + 6*(sigma**6/dist**7)) #Derivada 
                Fx -= dVdr * (dx/dist) #Calculamos las componentes de la fuerza
                Fy -= dVdr * (dy/dist)
                Fz -= dVdr * (dz/dist)
                
        enerind.append(e/2)
        Find.append(0.5*np.sqrt((Fx**2+Fy**2+Fz**2)))
        FX.append(Fx)
        FY.append(Fy)
        FZ.append(Fz)
        
    elif mode == 'finito':
      for i in range(N-1):
        e = 0
        Fx, Fy, Fz = (0,0,0)
        for j in range(i+1,N):
          dx = X[i] - X[j]
          dy = Y[i] - Y[j]
          dz = Z[i] - Z[j]
    
          dist = np.sqrt(dx**2 + dy**2 + dz**2)
    
          if dist < rc:
            V = 4*epsilon*((sigma/dist)**12-(sigma/dist)**6) #Calculamos el potencial
            total += V
            e += V
            dVdr = 4*epsilon*(-12*(sigma**12*dist**-13) + 6*(sigma**6*dist**-7)) #Derivada 
            Fx -= dVdr * (dx/dist) #Calculamos las componentes de la fuerza
            Fy -= dVdr * (dy/dist)
            Fz -= dVdr * (dz/dist)
            
        enerind.append(e)
        Find.append((Fx**2+Fy**2+Fz**2)**0.5)
        FX.append(Fx)
        FY.append(Fy)
        FZ.append(Fz)
        
    return enerind, total, np.array([FX,FY,FZ])
#------------------------------------------------------------------------------
def Mov3D(Coord, Vel, Fvec, ini, time, h, pred = 3.603, modo = 'infinito'):
    t0,tf = time
    x, y, z = Coord
    vx, vy, vz = Vel
    m, irc,tamaño = ini
    N = int((tf-t0)/h)
    T = np.linspace(t0,tf,N+1)
    vectorx, vectory, vectorz, vectorvx, vectorvy, vectorvz = ([x],[y],[z],[vx],[vy],[vz])
    
    Fx,Fy,Fz = Fvec
    Ax,Ay,Az = (np.array(Fx)/m,np.array(Fy)/m,np.array(Fz)/m)
    Epot = [Calcpot([x,y,z],tamaño,irc,mode = modo)[1]]
    FUE = [Fvec]
    for i in range(N):
        xt1 = x + vx * h + Ax* h**2 /2
        yt1 = y + vy * h + Ay* h**2 /2
        zt1 = z + vz * h + Az* h**2 /2
        
        Coord = [xt1,yt1,zt1]
        potF = Calcpot(Coord,tamaño,irc,mode = modo)
        potenciales = potF[1]
        fuerzast1 = potF[2]
        Fxt1,Fyt1,Fzt1 = fuerzast1
        Axt1, Ayt1, Azt1 = np.array(Fxt1)/m , np.array(Fyt1)/m, np.array(Fzt1)/m
        vxt1 = vx + (Axt1 + Ax) *h/2
        vyt1 = vy + (Ayt1 + Ay) *h/2
        vzt1 = vz + (Azt1 + Az) *h/2
    
        vectorx.append(xt1)
        vectory.append(yt1)
        vectorz.append(zt1)
        vectorvx.append(vxt1)
        vectorvy.append(vyt1)
        vectorvz.append(vzt1)
        Epot.append(potenciales)
        FUE.append(fuerzast1)
        
        x,y,z,vx,vy,vz = (xt1,yt1,zt1,vxt1,vyt1,vzt1)
        Ax,Ay,Az = Axt1,Ayt1,Azt1
        
    return T, (vectorx,vectory,vectorz), (vectorvx,vectorvy,vectorvz), Epot, FUE

#------------------------------------------------------------------------------
def Temperatura(velocidades,K, mode = 'vector'):
    if mode == 'vector':
        N = len(velocidades)
        ECIN = 0
        for i in velocidades:
            ECIN += m/2 * (i[0]**2+i[1]**2+i[2]**2)
            
          
        T = 2/(3*K*N) * ECIN
    elif mode == 'listas':
        N = len(velocidades[0])
        ECIN = 0
        for i in range(N):
            ECIN += m/2 *(velocidades[0][i]**2 + velocidades[1][i]**2 + velocidades[2][i]**2)
            
        T = 2/(3*K*N) * ECIN
        
    return T, ECIN
#------------------------------------------------------------------------------
def ColisCob(Coord0,Vel0, Coord, Vel, Fvec, ini, time, h, pred = 3.603, modo = 'infinito'):
    coordnew = []
    velnew = []
    Fnew = []
    for i in range(len(Coord)):
        coordnew.append(np.append(Coord[i],Coord0[i]))
        velnew.append(np.append(Vel[i],Vel0[i]))
        Fnew.append(np.append(Fvec[i],0))
        
    return Mov3D(np.array(coordnew), np.array(velnew), np.array(Fnew), ini, time, h, pred, modo)
    
    
def AnimarRed(frame,X,Y,Z): #Esta función sirve para actualizar valores en cada frame
    Xrep,Yrep,Zrep = (X[frame],Y[frame],Z[frame])
    imagen._offsets3d = (Xrep,Yrep,Zrep)
    titulo.set_text('Evolución red cristalina, tiempo {0} s'.format(round(tiempos[frame],15)))
    
    return imagen, titulo

#------------------------------------------------------------------------------
m, K,sigma,epsilon,pred = (1.05e-25/16, 8.6e-5, 2.3151, 0.167, 3.603)

#Ahora vamos a ver si me hace el movimiento de las partículas
np.random.seed(1) #Vamos a fijar la semilla aleatoria
Prueba = RedCristalina(2,2,2,300,(m,K), dist='FCC') #Calculamos el estado inicial de la red
tamaño = (2,2,2)
temper = Temperatura(Prueba[3],K)
print(temper)
Potencial = Calcpot(Prueba[:3],tamaño,3, mode = 'infinito') #Calculamos la energía potencial inicial de la red
Fuerzas0 = Potencial[2] #Fuerzas iniciales de la red
irc = 3 #Aproximación a radio 3, terceros venicos
velocidades = Prueba[3] #Velocidades iniciales (aleatorias)
vx,vy,vz = ([],[],[])
for i in Prueba[3]: #Separamos los componentes de las velocidades
    vx.append(i[0])
    vy.append(i[1])
    vz.append(i[2])

#Definimos parámetros iniciales para Mov3D y los tiempos a usar  
parametros = (m , irc, tamaño )
tiempo=(0,1e-12)
Coord = (Prueba[0],Prueba[1],Prueba[2])
Vel = (np.array(vx),np.array(vy),np.array(vz))
h = 1e-15

#Evaluamos la función
TEST = Mov3D(Coord, Vel, Fuerzas0, parametros,tiempo,h)

#Desempaquetamos los datos
tiempos = TEST[0]
XT = TEST[1][0]
YT = TEST[1][1]
ZT = TEST[1][2]
VXT = TEST[2][0]
VYT = TEST[2][1]
VZT = TEST[2][2]

#Calculamos temperatura, energía cinética y potencial para las diferentes posiciones.
Temperaturas = []
Ecin = []
Epot = TEST[3]
for i in range(len(tiempos)):
    TeCIN = Temperatura([VXT[i],VYT[i],VZT[i]],K, mode = 'listas')
    Ecin.append(TeCIN[1])
    Temperaturas.append(TeCIN[0])
    
Fuerzas = TEST[4]   
    
plt.plot(tiempos,Temperaturas,label = 'temperatura')
plt.xlabel('tiempo (s)')
plt.ylabel('Temperatura (K)')
plt.title('Temperatura de la red')
plt.show()

EpotC = (Epot-Epot[0])
plt.plot(tiempos,Ecin,label = 'E. cin.')
plt.plot(tiempos,EpotC,label = 'E. pot.')
plt.plot(tiempos,np.array(Ecin)+EpotC, label = 'E. tot.')
plt.legend(loc='best')
plt.xlabel('tiempo (s)')
plt.ylabel('energía (eV)')
plt.title('Comparativa energías cinética y potencial')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.scatter(XT[0],YT[0],ZT[0])
ax.set_xlabel('Eje x')
ax.set_ylabel('Eje y')
ax.set_zlabel('Eje z')
plt.show()

#Vamos a representar el movimiento de una partícula
pasos = len(XT)
Xp1 = [XT[i][0] for i in range(pasos)]
Yp1 = [YT[i][0] for i in range(pasos)]
Zp1 = [ZT[i][0] for i in range(pasos)]

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot(Xp1,Yp1,Zp1)
ax.set_xlabel('Eje x (Ams)')
ax.set_ylabel('Eje y (Ams)')
ax.set_zlabel('Eje z (Ams)')
ax.set_title('movimiento de una partícula respecto a su posición inicial')
plt.show()

#¿Qué frames debemos coger?
#Lo definimos nosotros, cada cuantos frames calculados vamos a animar
fac = int(3)
index = []
#Cogemos los tiempos a representar y los índices de los mismos
for i in range(1,int(len(tiempos)/fac+1)):
    index.append(fac*i-fac)

#Vamos a probar a actualizar valores así:
Xrep, Yrep, Zrep = (XT[0],YT[0],ZT[0])
#Vamos a mostrar el cambio de posiciones en el tiempo
#Aquí tenemos que animar
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

#Introducimos los datos iniciales
imagen = ax.scatter(Xrep,Yrep,Zrep)

# Propiedades de la figura
titulo = ax.set_title('Evolución de red cristalina')
ax.set_xlabel('X (Am)')
ax.set_ylabel('Y (Am)')
ax.set_zlabel('Z (Am)')

#Llamamos a la función animación, la función de actualización
#Pasamos nuestra lista de índices para los frames de tiempo que queremos
#Las coordenadas son argumentos para actualizar valores
#Cambiamos el intervalo entre imágenes a 50 ms para que se vea bien que se mueve
animacion = animation.FuncAnimation(fig,AnimarRed,fargs = (XT,YT,ZT),frames = index
                                    ,repeat = True, interval = 50)

plt.show()