Física computacional

PROGRAMA PARA LA SIMULACIÓN DE REDES CRISTALINAS

FUNCIONES DEFINIDAS
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