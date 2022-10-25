# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 12:34:56 2022

@author: Simple
"""

import os
import numpy as np
from tqdm import trange
from scipy.fft import fft

temps=[200,300,400]
intervalos=[20,10,5]

def espectrografia_LAMMPS(direccion, temperatura, ts):
    
    #Correr la simulacion con los parametros especificados
    
    print('Paso 0: Corre simulacion')
    os.system('cd ' + direccion + ' && lmp -in in_2nm_1ratio50_{}K_{}TS.airebo'.format(temperatura,ts))
    print('Simulacion terminada')
    
    vel=open(direccion+r'\vel_2nm_1ratio50_{}K_{}TS.dat'.format(temperatura,ts),'r') #Archivo de velocidades de LAMMPS
    info=vel.readlines() #Lineas del archivo
    n_atoms=int(info[3]) #Numero de atomos provisto por el archivo
    time=10000 #Este valor lo ponemos nosotros en la tarjeta de entrada para LAMMPS
    timestep=ts #Cada timestep se toma una medida de las velocidades de las particulas
    timesteps=int(time/timestep) #Numero de instantes de tiempo en el que medimos las velocidades

    vel_coordenadas=np.ndarray(shape=(timesteps,n_atoms,3)) #El array donde vienen todas las velocidades
    
    print('Paso 1: Guarda las velocidades')

    times=np.zeros(timesteps)
    for t in trange(timesteps):
        times[t]=t*timestep
        for n in range(n_atoms):
            v_atom=info[9+n+(9+n_atoms)*t].split()
            for c in range(3):
                vel_coordenadas[t][n][c]=v_atom[c]
                 
    
    print('Paso 2: Calcula la funcion de autocorrelacion para todos los tiempos')
    
    #Nota: Siguiendo el codigo de fortran se calcula la funcion de autocorrelacion
    
    Z=np.zeros(timesteps)   #Funcion de autocorrelacion dependiente del tiempo
    
    print('Paso 2.1: Define el denominador de la funcion de autocorrelacion')

    denominador=0

    for t0 in trange(timesteps):
        for n in range(n_atoms):
            denominador+=np.sqrt(vel_coordenadas[t0][n][0]**2+vel_coordenadas[t0][n][0]**2+vel_coordenadas[t0][n][0]**2)
            
    print('Paso 2.2: Define el numerador de la funcion de autocorrelacion para todos los tiempos')

    numerador=np.zeros(timesteps)

    for t in trange(timesteps):
        num=0
        for t0 in range(timesteps):
            for n in range(n_atoms):
                if t+t0>499:
                    num+=0
                else:
                    num+=np.sqrt(abs(vel_coordenadas[t0+t][n][0]*vel_coordenadas[t0][n][0])+abs(vel_coordenadas[t0+t][n][1]*vel_coordenadas[t0][n][1])+abs(vel_coordenadas[t0+t][n][2]*vel_coordenadas[t0][n][2]))
        numerador[t]=num
        
    Z=numerador/denominador
    
    fft_funcion_autocorrelacion=fft(Z)
    
    onda_vs_intensidad=open(direccion+r'\i_vs_k_2nm_1ratio50_{}K_{}TS.dat'.format(temperatura,ts),"w")
    
    norm = [float(i)/max(fft_funcion_autocorrelacion) for i in fft_funcion_autocorrelacion]

    for i in range(len(times)):
        onda_vs_intensidad.write(str(times[i]) + " " + str(norm[i].real) + "\n")
        
    onda_vs_intensidad.close()

espectrografia_LAMMPS(r'C:\Users\Simple\Documents\U\LAMMPS\Simul+Datos\Simul full',350,20)
espectrografia_LAMMPS(r'C:\Users\Simple\Documents\U\LAMMPS\Simul+Datos\Simul full',350,10)
espectrografia_LAMMPS(r'C:\Users\Simple\Documents\U\LAMMPS\Simul+Datos\Simul full',350,5)
espectrografia_LAMMPS(r'C:\Users\Simple\Documents\U\LAMMPS\Simul+Datos\Simul full',300,20)
espectrografia_LAMMPS(r'C:\Users\Simple\Documents\U\LAMMPS\Simul+Datos\Simul full',300,10)
espectrografia_LAMMPS(r'C:\Users\Simple\Documents\U\LAMMPS\Simul+Datos\Simul full',300,5)