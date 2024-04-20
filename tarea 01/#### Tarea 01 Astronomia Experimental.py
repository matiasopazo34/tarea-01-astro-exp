#### Tarea 01 Astronomia Experimental
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

def RMS(valores):
    valor2 = 0
    for i in range(0, len(valores)):
        valor2 += (valores[i])**2
    return np.sqrt(valor2/len(valores))

# parte B
path = 'C:\\Users\\matia\\Desktop\\Astro Exp\\tarea-01-astro-exp\\tarea 01\\antdip_AE2023A.xlsx'
AD_data = pd.read_excel(path, sheet_name='Sheet1')
#print(AD_data['-SECZ'])
AD_data = AD_data.drop(['tau', 'valor_tau'], axis=1)
AD_data = AD_data.dropna()
# print(AD_data.head()) #base de datos limpia

sec = AD_data['-SECZ'].to_numpy()
ln_dP = AD_data['ln(dP)'].to_numpy()
plt.plot(sec, ln_dP)

poly = np.polyfit(sec, ln_dP, 1) # con polyfit se aproxima un polinomio
print(poly[0])

def recta(x,m,n):
    return m*x + n

ajuste_lineal = recta(sec, poly[0], poly[1])
plt.plot(sec, ajuste_lineal)
plt.show()

#PARTE C

path = r'C:\\Users\\matia\\Desktop\\Astro Exp\\tarea-01-astro-exp\\tarea 01\\datos_espectros_2023_1' # ojo con el r al inicio, es importante
set_datos = []
for i in range(11,26):
    filename = "sdf_1"+str(i)+"_1"+str(i)
    v, T = np.genfromtxt(path + '\\' + filename, unpack=True, skip_header=108)
    set_datos.append([v,T])

T_promedios = []
RMS_total = []
for i in range(0,5):
    T_acumulado = []
    RMS_seq = []
    for N in range(0, 3): #PARA UN DATO COINCIDENTE
        v_1, T_1 = set_datos[i + N*5] # datos coincidentes
        T_acumulado.append(T_1)
        T_p = sum(T_acumulado)/len(T_acumulado)
        ruido = sigma_clip(T_p, masked=False) # que devuelve sigma? xd
        RMS_i = RMS(ruido)
        RMS_seq.append(RMS_i)
    RMS_total.append(RMS_seq)
    T_promedios.append(T_p)

for i in range(0, len(RMS_total)):
    N = [1, 2, 3]
    plt.plot(N, RMS_total[i])
plt.title("Gráfica de RMS vs N para las 5 fuentes observadas")
plt.show()


#gráfica de los 5 promedios 
v, T_sin_usoxd = set_datos[0] 
for i in range(0, len(T_promedios)):
    plt.plot(v, T_promedios[i], label="T "+str(i))
plt.title("Gráfica de T_promediados")
plt.legend()
plt.show()

# ahora fitearemos gaussiana a los promedios
# Se define funcion gaussiana a fitear
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
    plt.plot(v,f_gauss(v,t0, M,S), label='Fiteo '+str(i+1)) # se grafica bajo el gauss
plt.title('Espectro con gauss', fontsize=18)
plt.ylabel('Temperatura [K]', fontsize=18)
plt.xlabel(r'Velocidad [$\frac{km}{s}$]', fontsize=18)
plt.legend()
plt.show()

#Hacer la cruz de T_max
T_max = []
for i in range(0, len(T_promedios)):
    T = T_promedios[i] 
    T_max_i = np.amax(T)
    T_max.append(T_max_i)
cruz = np.array([[0,T_max[0],0], [T_max[1],T_max[2],T_max[3]], [0,T_max[4],0]])
print(cruz)
coord_x_lii = [208.863495, 208.996002, 209.128510] # 12, 13, 14
coord_y_bii = [-19.260527, -19.385527, -19.510527] # 11, 13, 15

#Gauss para coordenada x T_max
T_max_fit = [T_max[1],T_max[2],T_max[3]]
fg = [29.19566667, 208.996002, 1]
coefs,cov = curve_fit(f_gauss,coord_x_lii,T_max_fit, p0=fg) 
t0,M,S = coefs[0],coefs[1],coefs[2] 
X = np.linspace(208,210, 500)
plt.plot(X,f_gauss(X,t0, M,S), label='Fiteo coordenada x T_max') 
T_maxx = np.amax(f_gauss(X,t0, M,S))
for i in range(0, len(X)):
    if f_gauss(X[i],t0, M,S) == T_maxx:
        X_max = X[i]
plt.title('Fiteo coordenada x T_max')
plt.show()


#Gauss para coordenada y T_max
T_max_fit = [T_max[0],T_max[2],T_max[4]]
fg = [29.19566667, -19.385527, 1]
coefs,cov = curve_fit(f_gauss,coord_y_bii,T_max_fit, p0=fg) 
t0,M,S = coefs[0],coefs[1],coefs[2] 
X = np.linspace(-20,-19, 500)
plt.plot(X,f_gauss(X,t0, M,S), label='Fiteo coordenada y T_max') 
T_maxy = np.amax(f_gauss(X,t0, M,S))
for i in range(0, len(X)):
    if f_gauss(X[i],t0, M,S) == T_maxy:
        Y_max = X[i]
plt.title('Fiteo coordenada y T_max')
plt.show()
#calcular/buscar puntos máximos 
print(X_max, Y_max)
plt.plot(X_max, Y_max, marker="o", color="red", label='ubicación real de orión bajo T_max')
plt.title('ubicación real de orión bajo T_max')
plt.axhline(y = -19.385527) 
plt.axvline(x = 208.996002) 
plt.legend()
plt.show()


#repetir lo mismo con T_int
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
plt.plot(X,f_gauss(X,t0, M,S), label='Fiteo coordenada x T_int') 
T_maxx = np.amax(f_gauss(X,t0, M,S))
for i in range(0, len(X)):
    if f_gauss(X[i],t0, M,S) == T_maxx:
        X_max = X[i]
plt.title('Fiteo coordenada x T_int')
plt.show()

#para bii
T_int_fit = [T_int[0],T_int[2],T_int[4]]
fg = [140.07634617, -19.385527, 1]
coefs,cov = curve_fit(f_gauss,coord_y_bii,T_int_fit, p0=fg) 
t0,M,S = coefs[0],coefs[1],coefs[2] 
X = np.linspace(-20,-19, 500)
plt.plot(X,f_gauss(X,t0, M,S), label='Fiteo coordenada y T_int') 
T_maxx = np.amax(f_gauss(X,t0, M,S))
for i in range(0, len(X)):
    if f_gauss(X[i],t0, M,S) == T_maxx:
        Y_max = X[i]
plt.title('Fiteo coordenada y T_int')
plt.show()
print(X_max, Y_max)
plt.plot(X_max, Y_max, marker="o", color="red",  label='ubicación real de orión bajo T_int')
plt.title('ubicación real de orión bajo T_int')
plt.axhline(y = -19.385527) 
plt.axvline(x = 208.996002) 
plt.legend()
plt.show()
