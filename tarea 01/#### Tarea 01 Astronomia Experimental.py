#### Tarea 01 Astronomia Experimental
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

def RMS(valores): # programación RMS para una lista de valores
    valor2 = 0
    for i in range(0, len(valores)):
        valor2 += (valores[i])**2
    return np.sqrt(valor2/len(valores))

# parte B
path = 'C:\\Users\\matia\\Desktop\\Astro Exp\\tarea-01-astro-exp\\tarea 01\\antdip_AE2023A.xlsx'
AD_data = pd.read_excel(path, sheet_name='Sheet1')
AD_data = AD_data.drop(['tau', 'valor_tau'], axis=1)
AD_data = AD_data.dropna()

el = AD_data['EL'].to_numpy()
ln_dP = AD_data['ln(dP)'].to_numpy()
sen = np.sin(np.radians(el))
plt.plot(1/sen, ln_dP, label='data entregada')
poly = np.polyfit(1/sen, ln_dP, 1) # con polyfit se aproxima un polinomio
def recta(x,m,n):
    return m*x + n
ajuste_lineal = recta(1/sen, poly[0], poly[1])
print('la pendiente corresponde a '+str(poly[0]))
plt.plot(1/sen, ajuste_lineal, label='ajuste lineal')
plt.title('Gráfica de 1/sen(EL) vs ln(dP)')
plt.xlabel('1/sen(EL)')
plt.ylabel('ln(dP)')
plt.legend()
plt.show()

#PARTE C
# exportación de la data
path = r'C:\\Users\\matia\\Desktop\\Astro Exp\\tarea-01-astro-exp\\tarea 01\\datos_espectros_2023_1' 
set_datos = [] # lista dónde se guardarán los datos de los espectros
for i in range(11,26):
    filename = "sdf_1"+str(i)+"_1"+str(i)
    v, T = np.genfromtxt(path + '\\' + filename, unpack=True, skip_header=108)
    set_datos.append([v,T])

T_promedios = [] # Aqui obtendremos los promedios y sus respectivos errores
RMS_total = []
for i in range(0,5): # Aquí recorremos los 5 espectros
    T_acumulado = []
    RMS_seq = []
    for N in range(0, 3): 
        v_1, T_1 = set_datos[i + N*5] # espectros de una misma coordenada
        T_acumulado.append(T_1) # acá iremos "acumulando" los espectros de una misma coordenada
        T_p = sum(T_acumulado)/len(T_acumulado)
        ruido = sigma_clip(T_p, masked=False) # obtenemos el ruido
        RMS_i = RMS(ruido)
        RMS_seq.append(RMS_i) # progresión del error c/r a N de una coordenada
    RMS_total.append(RMS_seq) # de la forma [coord1, coord2, coord3, coord4, coord5]
    T_promedios.append(T_p)

for i in range(0, len(RMS_total)):
    N = [1, 2, 3]
    plt.plot(N, RMS_total[i], label='coordenada '+str(i+1))
plt.title("Gráfica de RMS vs N para las 5 fuentes observadas")
plt.xlabel('cantidad de datos N')
plt.ylabel('Error RMS c/r al ruido')
plt.legend()
plt.show()

#gráfica de los 5 espectros promediados 
v, T_sin_usoxd = set_datos[0] 
for i in range(0, len(T_promedios)):
    plt.plot(v, T_promedios[i], label="Espectro coordenada "+str(i+1))
plt.title("Gráfica de los 5 espectros promediados")
plt.ylabel('Temperatura [K]', fontsize=12)
plt.xlabel(r'Velocidad [$\frac{km}{s}$]', fontsize=12)
plt.legend()
plt.show()

# Se define la funcion gaussiana para fitear
def f_gauss(x,T0,mean,stdv):
    return T0*np.exp(-((x-mean)**2)/(2*(stdv**2)))
for i in range(0, len(T_promedios)):
    T = T_promedios[i]
    T_max = np.amax(T) 
    for j in range(0, len(T)):
        if T[j] == T_max:
            v_medio = v[j]
    fg = [T_max, v_medio, 1] # valores para fitear, T_max, v_medio, desviacion estandar
    coefs,cov = curve_fit(f_gauss,v,T, p0=fg) # Se fitea
    t0,M,S = coefs[0],coefs[1],coefs[2]  # Se extraen los coeficientes fiteados
    plt.plot(v,f_gauss(v,t0, M,S), label="Espectro coordenada "+str(i+1)) # se grafica bajo la funcion gaussiana

plt.title('Espectro con gauss', fontsize=18)
plt.ylabel('Temperatura [K]', fontsize=12)
plt.xlabel(r'Velocidad [$\frac{km}{s}$]', fontsize=12)
plt.legend()
plt.show()

# Ahora la gráfica de las temperaturas máximas vs las coordenadas
T_max = []
for i in range(0, len(T_promedios)):
    T = T_promedios[i] 
    T_max_i = np.amax(T) # extraemos las temperaturas máximas
    T_max.append(T_max_i)

# la cruz es una representación simbólica de los datos representados con sus coordenadas
cruz = np.array([[0,T_max[0],0], [T_max[1],T_max[2],T_max[3]], [0,T_max[4],0]])
print(cruz)
coord_x_lii = [208.863495, 208.996002, 209.128510] # 12, 13, 14
coord_y_bii = [-19.260527, -19.385527, -19.510527] # 11, 13, 15

# se realiza Gauss con respecto a la coordenada lii
T_max_fit = [T_max[1],T_max[2],T_max[3]]
fg = [29.19566667, 208.996002, 1]
coefs,cov = curve_fit(f_gauss,coord_x_lii,T_max_fit, p0=fg) 
t0,M,S = coefs[0],coefs[1],coefs[2] 
X = np.linspace(208,210, 500)
plt.plot(X,f_gauss(X,t0, M,S), label='Fiteo de temperatura máxima bajo coordenada lii') 
T_maxx = np.amax(f_gauss(X,t0, M,S))
for i in range(0, len(X)):
    if f_gauss(X[i],t0, M,S) == T_maxx:
        X_max = X[i] # Se obtiene el punto máximo para lii
plt.title('Fiteo de temperatura máxima bajo coordenada lii')
plt.ylabel('Temperatura [K]', fontsize=12)
plt.xlabel('coordenada lii', fontsize=12)
plt.show()

# se realiza Gauss con respecto a la coordenada bii
T_max_fit = [T_max[0],T_max[2],T_max[4]]
fg = [29.19566667, -19.385527, 1]
coefs,cov = curve_fit(f_gauss,coord_y_bii,T_max_fit, p0=fg) 
t0,M,S = coefs[0],coefs[1],coefs[2] 
X = np.linspace(-20,-19, 500)
plt.plot(X,f_gauss(X,t0, M,S), label='Fiteo de temperatura máxima bajo coordenada bii') 
T_maxy = np.amax(f_gauss(X,t0, M,S))
for i in range(0, len(X)):
    if f_gauss(X[i],t0, M,S) == T_maxy:
        Y_max = X[i] # Se obtiene el punto máximo para bii
plt.title('Fiteo de temperatura máxima bajo coordenada bii')
plt.ylabel('Temperatura [K]', fontsize=12)
plt.xlabel('coordenada bii', fontsize=12)
plt.show()

# Se obtiene el punto (lii_1, bii_1) bajo el cual debería estar orión con T_maximo
print('la coordenada máxima es '+str(X_max)+str(Y_max))
plt.plot(X_max, Y_max, marker="o", color="red", label='ubicación real de orión')
plt.title('ubicación real de orión analizando la temperatura máxima')
plt.xlim(208,210)
plt.ylim(-19.45, -19.30)
plt.xlabel('coordenada lii')
plt.ylabel('coordenada bii')
plt.axhline(y = -19.385527) 
plt.axvline(x = 208.996002) 
plt.legend()
plt.show()


## Ahora la gráfica de las temperaturas integradas vs las coordenadas
T_int = []
v, T_sin_usoxd = set_datos[0] 
for i in range(0, len(T_promedios)):
    T = T_promedios[i]
    T_int_i = np.trapz(np.flip(T_promedios[i]), x=np.flip(v))
    T_int.append(T_int_i)
cruz = np.array([[0,T_int[0],0], [T_int[1],T_int[2],T_int[3]], [0,T_int[4],0]])
print(cruz)
coord_x_lii = [208.863495, 208.996002, 209.128510] # 12, 13, 14
coord_y_bii = [-19.260527, -19.385527, -19.510527] # 11, 13, 15

#para lii
T_int_fit = [T_int[1],T_int[2],T_int[3]]
fg = [140.07634617, 208.996002, 1]
coefs,cov = curve_fit(f_gauss,coord_x_lii,T_int_fit, p0=fg) 
t0,M,S = coefs[0],coefs[1],coefs[2] 
X = np.linspace(208,210, 500)
plt.plot(X,f_gauss(X,t0, M,S), label='Fiteo de temperatura máxima bajo coordenada lii') 
T_maxx = np.amax(f_gauss(X,t0, M,S))
for i in range(0, len(X)):
    if f_gauss(X[i],t0, M,S) == T_maxx:
        X_max = X[i]
plt.title('Fiteo de temperatura integrada bajo coordenada lii')
plt.ylabel('Temperatura [K]', fontsize=12)
plt.xlabel('coordenada lii', fontsize=12)
plt.show()

#para bii
T_int_fit = [T_int[0],T_int[2],T_int[4]]
fg = [140.07634617, -19.385527, 1]
coefs,cov = curve_fit(f_gauss,coord_y_bii,T_int_fit, p0=fg) 
t0,M,S = coefs[0],coefs[1],coefs[2] 
X = np.linspace(-20,-19, 500)
plt.plot(X,f_gauss(X,t0, M,S), label='Fiteo de temperatura máxima bajo coordenada bii') 
T_maxx = np.amax(f_gauss(X,t0, M,S))
for i in range(0, len(X)):
    if f_gauss(X[i],t0, M,S) == T_maxx:
        Y_max = X[i]
plt.title('Fiteo de temperatura integrada bajo coordenada bii')
plt.ylabel('Temperatura [K]', fontsize=12)
plt.xlabel('coordenada bii', fontsize=12)
plt.show()

## Se obtiene el punto (lii_1, bii_1) bajo el cual debería estar orión con T_integrado
print('la coordenada integrada es '+str(X_max)+str(Y_max))
plt.plot(X_max, Y_max, marker="o", color="red",  label='ubicación real de orión')
plt.title('ubicación real de orión analizando la temperatura integrada')
plt.xlim(208,210)
plt.ylim(-19.45, -19.30)
plt.xlabel('coordenada lii')
plt.ylabel('coordenada bii')
plt.axhline(y = -19.385527) 
plt.axvline(x = 208.996002) 
plt.legend()
plt.show()


