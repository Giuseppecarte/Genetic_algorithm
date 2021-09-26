import numpy as np
from tensorflow.random import shuffle
import pandas as pd 
import matplotlib.pyplot as plt

#Definimos las variables globales 
prob_cruza = 0.9
tam_pob = 100
bits = 128
prob_mut = 1/bits
G_max = 100
gamma = [0 if i%2 == 0 else 1 for i in range(bits)]
#Población 0
Pinic = [list(np.random.randint(0,2,bits)) for i in range(tam_pob)]

##Estructura General --------------------------------------------
def Hamming(genotipo, gamma = gamma):
    assert len(genotipo)==len(gamma)
    peso = 0
    for j in range(len(genotipo)):
        if genotipo[j] == gamma[j]:
            peso += 1
        else: pass
    return peso

def padres(Pinic):
    #Hacemos las mezclas
    P1 = shuffle(Pinic).numpy().tolist()
    P2 = shuffle(P1).numpy().tolist() 
    #parte 2
    def pelea_padres(padres):
        padres_ganadores = []
        for j in list(range(0,len(padres),2)):
            if Hamming(padres[j]) <= Hamming(padres[j+1]):
                padres_ganadores.append(padres[j])
            else: padres_ganadores.append(padres[j+1])
        return padres_ganadores
    p1_prima = pelea_padres(P1)
    p2_prima = pelea_padres(P2)
    return p1_prima, p2_prima

p1_prima, p2_prima = padres(Pinic)

#Q
def cruza(padres1, padres2):
    assert len(padres1)==len(padres2)
    Q = []
    for i in range(len(padres1)):
        num_aleatorio = np.random.random(1)[0]
        #Aquí
        if num_aleatorio <= prob_cruza:
            aleatorio = np.random.randint(1, len(padres1))
            x = padres1[i]
            y = padres2[i]
            resto = abs(aleatorio - len(x))

            h1 = [*x[:aleatorio],*y[-resto:]]
            h2 = [*x[-resto:],*y[:aleatorio]]
            # original
            # h2 = [*y[:aleatorio],*x[-resto:]]
            Q.append(h1)
            Q.append(h2)
        else:
            Q.append(padres1[i])
            Q.append(padres2[i])
    return Q


def mutacion(Q):
    Q_prima = []
    for genotipo in Q:
        num_aleatorio = np.random.random(1)[0]
        #Aquí!!!!!!
        if num_aleatorio <= prob_mut:
            genotipo_prima = [0 if j==1 else 1 for j in genotipo ]
            Q_prima.append(genotipo_prima)
        else: Q_prima.append(genotipo)
    return Q_prima
#-------------------------------------------------------------

#Estrategias--------------------------------------------------
def mu_mas_l(Pinic=Pinic):
    conteo = 0
    result = []
    genotipos = []
    generaciones = []   
    while conteo < G_max:
        p1_prima, p2_prima = padres(Pinic)
        Q = cruza(p1_prima, p2_prima)
        concat = Pinic + Q
        eval_concat = [Hamming(line) for line in concat]
        datos = pd.DataFrame({'Genotipo':[line for line in concat], 'F(x)':eval_concat}).sort_values('F(x)')
        datos = datos.iloc[:int(len(concat)/2)]
        #print(f"Generación:{conteo} : {datos['Genotipo'].iloc[0]} : {datos['F(x)'].iloc[0]}")
        generaciones.append('Generación : '+str(conteo))
        genotipos.append(datos['Genotipo'].iloc[0])
        result.append(datos['F(x)'].iloc[0])
        Pinic = list(datos['Genotipo'])
        conteo += 1
    datos_finales = pd.DataFrame({'Generaciones':generaciones, 'Genotipos':genotipos, 'F(x)':result})
    return datos_finales

def m_l(Pinic = Pinic):
    conteo = 0
    result = []
    genotipos = []
    generaciones = []
    while conteo < G_max:
        p1_prima, p2_prima = padres(Pinic)
        Q_prima = mutacion(cruza(p1_prima, p2_prima))
        datos_1 = pd.DataFrame({'Genotipo':[genotipo for genotipo in Q_prima], 'F(x)':[Hamming(genotipo) for genotipo in Q_prima]}).sort_values('F(x)')
        #print(f"Generación:{conteo} : {datos_1['Genotipo'].iloc[0]} : {datos_1['F(x)'].iloc[0]}")
        generaciones.append('Generación : '+str(conteo))
        genotipos.append(datos_1['Genotipo'].iloc[0])
        result.append(datos_1['F(x)'].iloc[0])
        Pinic = list(datos_1['Genotipo'])
        conteo += 1
    datos_finales = pd.DataFrame({'Generaciones':generaciones, 'Genotipos':genotipos, 'F(x)':result})
    return datos_finales

def mu_l_elitismo(Pinic = Pinic):
    conteo = 0
    result = []
    genotipos = []
    generaciones = []
    while conteo < G_max:
        p1_prima, p2_prima = padres(Pinic)
        mejor_Pinic = [Hamming(genotipo) for genotipo in Pinic]
        datos_2 = pd.DataFrame({'Genotipo':[gen for gen in Pinic], 'F(x)':mejor_Pinic})
        #Sacamos el mejor de Pinic

        Q_prima = mutacion(cruza(p1_prima, p2_prima))
        datos_1= pd.DataFrame({'Genotipo':[genotipo for genotipo in Q_prima], 'F(x)':[Hamming(genotipo) for genotipo in Q_prima]}).sort_values('F(x)')
        datos_1 = datos_1.iloc[:-1]

        Q_bi_prima = datos_1.append(datos_2.iloc[0]).sort_values('F(x)')
        #print(f"Generación:{conteo} : {Q_bi_prima['Genotipo'].iloc[0]} : {Q_bi_prima['F(x)'].iloc[0]}")
        generaciones.append('Generación : '+str(conteo))
        genotipos.append(Q_bi_prima['Genotipo'].iloc[0])
        result.append(Q_bi_prima['F(x)'].iloc[0])


        Pinic = list(Q_bi_prima['Genotipo'])
        conteo +=1
    datos_finales = pd.DataFrame({'Generaciones':generaciones, 'Genotipos':genotipos, 'F(x)':result})
    return datos_finales
#-------------------------------------------------------------

#Gráficas ----------------------------------------------------
def graphs(datos_finales_1,datos_finales_2, datos_finales_3):
    plt.style.use('fivethirtyeight')
    plt.figure(figsize = [20,7])
    plt.plot(range(datos_finales_1.shape[0]),datos_finales_1['F(x)'],label = '$(\mu + \lambda)$',linewidth=1)
    plt.plot(range(datos_finales_2.shape[0]),datos_finales_2['F(x)'],label = '$(\mu \lambda)$',linewidth=1,color = 'green')
    plt.plot(range(datos_finales_2.shape[0]),datos_finales_3['F(x)'],label = '$(\mu \lambda) elitismo$',linewidth=1)
    plt.title('Resultados de las diferentes estrategias con '+str(bits)+' bits')
    plt.xlabel('Generaciones')
    plt.ylabel('$F(x)$')
    plt.legend()
    return plt.show()


def main():
    print('\n*********************************\n')
    print('Cargando las estrategias..........')
    print('\n*********************************\n')
    datos_finales_1 = mu_mas_l()
    datos_finales_2 = m_l()
    datos_finales_3 = mu_l_elitismo()
    print('Elige una de las siguientes estrategias: \n\n')
    print(' Estrategias                 Respuestas')
    print('* Mu + lambda                      (1)  ')
    print('* Mu Lambda                        (2)  ')
    print('* Mu Lambda con elitísmo           (3)   ')
    print('* Salir                            (0)')
    respuesta1 = input('\n Escribe tu respuesta -> ')

    if respuesta1 == '1':
        print(datos_finales_1)
        print(graphs(datos_finales_1,datos_finales_2, datos_finales_3))
        return main()
    if respuesta1 == '2':
        print(datos_finales_2)
        print(graphs(datos_finales_1,datos_finales_2, datos_finales_3))
        return main()
    if respuesta1 == '3':
        print(datos_finales_3)
        print(graphs(datos_finales_1,datos_finales_2, datos_finales_3))
        return main()
    if respuesta1 == '0':
        return None
    else: 
        print('No es una opción válida\n \n')
        return main()

print(main())