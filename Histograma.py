import matplotlib.pyplot as plot
import numpy as np

with open("hist.dat", "r") as file:
    data = file.read().splitlines()
    velocidades = [float(dato) for dato in data]

intervalos = np.histogram_bin_edges(velocidades, bins='auto')
intervalos_redondeados = np.round(intervalos)
bins=40
plot.hist(x=velocidades, bins=bins, color='#86CEF8', rwidth=0.85, density=True)
#plot.xlabel('Distancia r')
#plot.ylabel('Densidad de particulas g(r)')
plot.xlabel('Velocidad')
plot.ylabel('Frecuencia')

x = np.linspace(min(velocidades), max(velocidades), 100)  # Generar puntos en el rango de valores de edades
# Definir las funciones por partes
#funcion_parte1 =  2.0 * np.pi * x / 100  # Función para x en [0, 5]
#funcion_parte2 = lambda x: 2.0 * np.pi * x *(1-4/np.pi*np.arccos(10/(2*x)))/ 100  # Función para x en [5, 10]

# Calcular los valores de la función por partes
#y_parte1 = np.piecewise(x, [x >= 0, x < 5], [funcion_parte1, 0])
#y_parte2 = np.piecewise(x, [x >= 5, x <= 10], [funcion_parte2, 0])
#funcion = (1/(0.71*2.0*3.141592)**0.5)* np.exp(-x**2/(2.0*0.71))
funcion = (1/(0.71))*x* np.exp(-x**2/(2.0*0.71))

# Graficar la función
#plot.plot(x, y_parte1, color='red', label='x=[0,5]')
#plot.plot(x, y_parte2, color='red', label='x=[5,10]')
plot.plot(x, funcion, color='red')

plot.show() #dibujamos el histograma