# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import scipy.integrate as sci
import re
#url = '/Users/Madaim3b/Documents/INAOE/proyecto/Muestra_MgII/'
url='/Users/madai/OneDrive/Documentos/INAOE/proyecto/Muestra_MgII/'
lista = glob.glob(url+'muestra*')

maximos=[]

asimetrias=[]

v=[]

AIs=[]
L25=[]
CIs=[]
L80=[]
FWHMs=[]

vrestas=[]
fang=[]
fanch=[]

sumaflujos=[]
residuos=[]
IDs=[]


for i in range(len(lista)):
    #Para extraer el ID de los nombres de los archivos
    ID=int(re.findall("[0-9]+", lista[i])[0])
    IDs.append(ID)
    archivo = pd.read_csv(lista[i], delim_whitespace=True, skiprows=1,\
    names=['lambda','angosta','ancha','perfil'])
    sumaflujo=archivo['angosta']+archivo['ancha']
    archivo['sumaflujo']=sumaflujo #agrega columna a archivo original
    sumaflujos.append(sumaflujo)
    
    residuo=archivo['perfil']-archivo['sumaflujo']
    archivo['residuo']=residuo
    residuos.append(residuo)
    
    maximo = max(archivo['sumaflujo'])
    maximos.append(maximo) #?
    
    Fmed=maximo/2
    F1_4=maximo/4
    F8_10=maximo*8/10
    F3_4=maximo*3/4
    

    indice_max = np.argmax(archivo['sumaflujo'])
    Lambda_max = archivo['lambda'][indice_max]
    
    part1 = archivo[archivo['lambda']<= Lambda_max]
    
    part2 = archivo[archivo['lambda']> Lambda_max]
    
    part2_inv=part2.sort_values(by='lambda',ascending=False).\
    reset_index(drop=True) 
    
    Lambdamed1=np.interp(Fmed,part1['sumaflujo'],part1['lambda'])
    Lambdamed2=np.interp(Fmed,part2_inv['sumaflujo'],part2_inv['lambda'])
    
    FWHM=Lambdamed2-Lambdamed1
    
    Lambda1_4=np.interp(F1_4,part1['sumaflujo'],part1['lambda'])
    Lambda1_4_2=np.interp(F1_4,part2_inv['sumaflujo'],part2_inv['lambda'])
    Parametro2=Lambda1_4_2-Lambda1_4
    Lambda_entre2=Parametro2/2
    Lambda25=Lambda1_4+Lambda_entre2
    
    Lambda8_10=np.interp(F8_10,part1['sumaflujo'],part1['lambda'])
    Lambd8_10_2=np.interp(F8_10,part2_inv['sumaflujo'],part2_inv['lambda'])
    Parametro3=Lambd8_10_2-Lambda8_10
    Lambda8_entre2=Parametro3/2
    Lambda80=Lambda8_entre2+Lambda8_10 
    
    
    A=(Lambda25-Lambda80)/FWHM
    
    asimetrias.append(A) #?
    L25.append(Lambda25)
    L80.append(Lambda80)
    FWHMs.append(FWHM)


    c=3E5 
    L0=2798

    v=(c*(archivo['lambda']-L0))/L0 
    archivo['v']=v #Intento agregar columna de velocidades 
    v_peak = archivo['v'][indice_max]

    
    part1v = archivo[archivo['v']<= v_peak]
    
    part2v = archivo[archivo['v']> v_peak]
    
    part2_inv_v=part2v.sort_values(by='lambda',ascending=False).\
    reset_index(drop=True)
    
    vB1_4=np.interp(F1_4,part1['sumaflujo'],part1v['v'])
    vR1_4=np.interp(F1_4,part2_inv['sumaflujo'],part2_inv_v['v'])
    
    AI=(vR1_4+vB1_4-(2*v_peak))/(vR1_4-vB1_4)
    C1_4=(vR1_4+vB1_4)/2
    
    AIs.append(AI)
    
    vB3_4=np.interp(F3_4,part1['sumaflujo'],part1v['v'])
    vR3_4=np.interp(F3_4,part2_inv['sumaflujo'],part2_inv_v['v'])
    
    CI=(vR3_4-vB3_4)/(vR1_4-vB1_4)
    CIs.append(CI)
    
  
    indice_max_ang = np.argmax(archivo['angosta'])
    v_peak_ang = archivo['v'][indice_max_ang]
  
    indice_max_anch = np.argmax(archivo['ancha'])
    v_peak_anch = archivo['v'][indice_max_anch]
    
    vresta=v_peak_anch-v_peak_ang
    
    vrestas.append(vresta)
    
    
    f_ang=sci.trapz(archivo.angosta,x=archivo['lambda'])#integra la angosta  (x dice de donde a donde integrar)
    f_anch=sci.trapz(archivo.ancha,x=archivo['lambda'])

    fang.append(f_ang)
    fanch.append(f_anch)


#dl25=pd.DataFrame({'L25':L25})
#dl25.to_csv('datosMontecarlo.csv',index=False)


df=pd.DataFrame({'name':lista,'A':asimetrias,'AI':AIs,'CI':CIs,'dV':vrestas,'Fang':fang,'Fanch':fanch})
df.to_csv('parametros_asimetria.txt',sep=' ',index=False) #Creo un archivo csv para guaradar mis datos
bines=np.arange(-0.35,0.45,0.05)
plt.hist(asimetrias, color="gray",bins=bines)
plt.xlabel("A")
plt.ylabel("N")      
plt.title("Parámetro de Asimetrías")
plt.show()   
    

    
bines1=np.arange(-0.4,0.4,0.05)
plt.hist(AIs,color="gray",bins=bines1)
plt.xlabel("AI")   
plt.ylabel("N")   
plt.title("Índice de Asimetrías")
plt.show()  


plt.hist(vrestas,color="gray",bins=15,range=(-4000,3500))
plt.title("Diferencia de velocidades")
plt.show()  
    

binesc=np.arange(0,0.5,0.045)
plt.hist(CIs, color='pink',bins=binesc)
plt.xlabel("CIs")   
plt.ylabel("N")   
plt.title("Curtosis")
plt.savefig('curtosis.pdf', bbox_inches='tight') 
plt.show()

binesv=np.arange(-4000,3500,500)
plt.hist(vrestas, color='gray',bins=binesv)
plt.xlabel("dv [km/s]")   
plt.ylabel("N")   
plt.title("Diferencia de velocidad central linea ancha y angosta")
plt.show()
    
    
plt.plot(archivo['v'],archivo['perfil'], label ='perfil')
plt.plot(archivo['v'],archivo['ancha'], label ='Linea ancha')
plt.plot(archivo['v'],archivo['angosta'], label ='linea angosta')
plt.legend()
plt.xlabel("Velocidad [km/s]")   # Establece el título del eje x
plt.ylabel("Flujo")   # Establece el título del eje y
plt.show()

#url_1 = '/Users/Madaim3b/Documents/INAOE/proyecto/Muestra_MgII/muestra_VIAI_231'
url_1='/Users/madai/OneDrive/Documentos/INAOE/proyecto/Muestra_MgII/muestra_VIAI_974'

archivo = pd.read_csv(url_1,delim_whitespace=True,skiprows=1,\
                      names=['lambda','angosta','ancha','perfil'])  
plt.plot(archivo['lambda'],archivo['perfil'], label ='perfil')
plt.plot(archivo['lambda'],archivo['ancha'], label ='Linea ancha')
plt.plot(archivo['lambda'],archivo['angosta'], label ='linea angosta')
plt.title("Componentes línea de emisión Mg II (ID 974)")
plt.legend()
plt.xlabel("Longitud de onda (Å)")   
plt.ylabel("Flujo normalizado (u.a.)")   
#plt.savefig('objeto974.pdf', bbox_inches='tight') 
plt.show()




#url_1 = '/Users/Madaim3b/Documents/INAOE/proyecto/Muestra_MgII/muestra_VIAI_2033'
url_1='/Users/madai/OneDrive/Documentos/INAOE/proyecto/Muestra_MgII/muestra_VIAI_974'
archivo = pd.read_csv(url_1,delim_whitespace=True,skiprows=1,\
                      names=['lambda','angosta','ancha','perfil'])  
    
#df=pd.DataFrame({'name':lista,'sumaflujo':sumaflujos,'residuo':residuos})
#df.to_csv('residuo.txt',sep=' ',index=False) #? 

#df = pd.DataFrame(archivo,columns=[archivo['sumaflujo']]) 
#sumaflujo = [[archivo['sumaflujo']]

sumaflujo=archivo['angosta']+archivo['ancha']
archivo['sumaflujo']=sumaflujo

residuo=archivo['perfil']-archivo['sumaflujo']
archivo['residuo']=residuo


plt.plot(archivo['lambda'],archivo['perfil'], label ='perfil')
plt.plot(archivo['lambda'],archivo['sumaflujo'], label ='Ajuste')
plt.plot(archivo['lambda'],archivo['residuo'], label ='residuo')
plt.title("Ajuste línea de emisión Mg II (ID 974)")
plt.legend()
plt.xlabel("Longitud de onda (Å)")   
plt.ylabel("Flujo normalizado (u.a.)")   
#plt.savefig('objeto974_ajustecool.pdf', bbox_inches='tight') 
plt.show()



dfN=pd.DataFrame({'L25':L25,'L80':L80,'FWHM':FWHM })
dfN.to_csv('variables_parametros_asimetria.txt',sep=' ',index=False)





