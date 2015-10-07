# El programa simula una copolimerización lineal de dos tipos de monómero. 
# Esta build nunca creó materia en todas las pruebas que he hecho.
# Fue bastante divertido y frustrante.
# Disculpe el cochinero, no tengo funciones definidas porque no funcionaban. No tengo idea porqué. 
# Y sí, fue tan difícil debuggearlo como parece.
# Nunca pude hacer que contara el número máximo de monómeros consecutivos.

import random # Para poder usar numeros al azar.
import math # Para hacer calculos estadísticos.
import re # Para los datos curiosos.
import itertools

# Los componentes deben de ser strings.
A="A" # Monómero A.
B="B" # Monómero B.
I="I" # Iniciador.


i=int(input("[I] = ")) # Concentración del componente X, [X].
a=int(input("[A] = "))
b=int(input("[B] = "))
MM_I=float(input("Masa molecular I = "))
MM_A=float(input("Masa molecular A = "))
MM_B=float(input("Masa molecular B = "))
Lambda=int(input("Lambda = ")) # Longitud de la cadena cinética.

# No quiero que el programa marque error cuando se equivoca el usuario, sino que deje al usuario intentar de nuevo.
print("\nAegúrese que las probabilidades sean iguales a 1. Por cuestiones del programa, ciertas combinaciones de doubles\nno pueden ser comparadas con un double de 1.\n")
while 1==1:
    prob_dism=float(input("P(Dismutación) = "))
    prob_trans=float(input("P(Transferencia) = "))
    prob_recomb=float(input("P(Recombinación) = "))
    if (prob_dism<0 or prob_trans<0 or prob_recomb<0):
        print("No puede haber probabilidades negativas. Esto no es cuántica.")
    elif prob_dism>=0 and prob_trans>=0 and prob_recomb>=0:
        break

print("\nCuando hay muchas cadenas y una diferencia muy grande entre constantes, el programa tarda en terminar de correr.\nCuando pasa esto es recomendable incrementar el número componentes (manteniendo las proporciones deseadas)\no reducir la brecha entre las constantes dividiendo ambas entre un mismo número arbitrario.\n")
while 1==1:
    kia=float(input("K_IA = ")) # Los coeficientes de velocidad K_ij. 
    kib=float(input("K_IB = "))
    if kia==0 and kib==0:
        print("\nUna de las dos constantes de iniciación debe ser mayor que cero,\nde lo contrario no se lleva a cabo la reacción.\n")
    else:
        break
while 1==1:
    kaa=float(input("K_AA = "))
    kab=float(input("K_AB = "))
    kbb=float(input("K_BB = "))
    kba=float(input("K_BA = "))
    if kaa==0 and kab==0 and kbb==0 and kba==0:
        print("Alguna de las cosntantes de propagación debe ser mayor que cero,\n de lo contrario no hay propagación.")
    else:
        break
datoscuriosos=(input("Quiere ver datos curiosos? y/n? "))
if datoscuriosos == "y":
    rep_a=int(input("Cuantas repeticiones de A quiere que cuente el programa? "))
    rep_b=int(input("Cuantas repeticiones de B quiere que cuente el programa? "))
    rep_ia=int()
    rep_ib=int()
    rep_ai=int()
    rep_bi=int()
    rep_ab=int()
    rep_ba=int()
    maxrepa=int()
    maxreob=int()
else:
    datoscuriosos=="n"

# Coeficientes y probabilidades para hacer los cálculos. Probabilidades Bayesianas.
ra=kaa/kab
rb=kbb/kba
p_aa=ra*a/(ra*a+b) # Probabilidad que reaccione A con A.
p_ab=1-p_aa
p_bb=rb*b/(rb*b+a)
p_ba=1-p_bb
p_ia=kia*a/(kia*a+kib*b) # Probabilidad que el iniciador reaccione con A.
p_ib=1-p_ia

# Variables estrella que sacan de apuros y ciclos. "C- C- C- COMBO BREAKER!" - http://www.youtube.com/watch?v=ezdptxvyAxo
CcccomboBreaker=5
CcccomboBreaker2=5

# Arreglos unidimensionales para guardar las cadenas para despus iterar sobre ellas y hacer estadísticas.
cadenas_dism=[]
cadenas_trans=[]
cadenas_recomb=[]
cadena=I

# Comienza lo bueno.
while i>0 and (a>0 or b>0) or (len(cadena)>1 and(a>0 or b>0)): # El ciclo no empieza si no hay iniciador o si no hay una cadena que este siendo formada por transferencia.
    CcccomboBreaker=5 # Reseteo de variable para que no interfiera.
    if CcccomboBreaker2==3: # Sirve para que no se hagan las transferencias no se hagan ciclos infinitos. Así termina forzadamente.
        break     
        # Rutina de terminación. Se repite cada vez que termina una cadena.
    if 1<=-1 or (a<=0 and b<=0): # En este caso es para que el programa cuente la cadena cuando ya se acabaron los reactivos.
        if CcccomboBreaker==1:
            break
        while len(cadenas_trans)==0 and len(cadenas_dism)==0: # dismutación
            cadenas_dism.append(cadena)
            CcccomboBreaker=2
            break
        x_t=random.uniform(0,1)
        while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
            cadenas_dism.append(cadena)
            CcccomboBreaker=2
            break
        while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0):#transferencia
            CcccomboBreaker=2
            CcccomboBreaker2=8
            x_r_int=random.randint(0,1)                        
            if a<=0 and b<=0:
                cadenas_trans.append(cadena)
                CcccomboBreaker2=3
                break
            elif len(cadenas_dism)>0 and x_r_int==0:
                if len(cadenas_dism)==1:
                    x_p_c=0
                elif len(cadenas_dism)>1:
                    x_p_c=random.randint(0,len(cadenas_dism)-1)
                cadenas_trans.append(cadena)
                cadena2=cadenas_dism[x_p_c]                        
                cadena=cadena2+"QSY" 
                cadenas_trans.append(cadena)
                cadenas_dism.pop(x_p_c)
                break
            elif len(cadenas_trans)>0 and x_r_int==1:
                if len(cadenas_trans)==1:
                    x_p_c=0
                elif len(cadenas_trans)>1:
                    x_p_c=random.randint(0,len(cadenas_trans)-1) 
                    cadenas_trans.append(cadena)
                    cadena2=cadenas_trans[x_p_c]                        
                    cadena=cadena2+"QSY"
                    cadenas_trans.append(cadena)
                    cadenas_trans.pop(x_p_c)
                    break
                else:
                    continue
        while x_t>(1-prob_recomb): #recombinacion
            CcccomboBreaker=2
            CcccomboBreaker2=3 #
            x_r_int=random.randint(0,1) 
            if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                cadenas_dism.append(cadena) 
                break                     
            elif len(cadenas_dism)>0 and x_r_int==0:
                if len(cadenas_dism)==1:
                    x_p_c=0
                elif len(cadenas_dism)>1:
                    x_p_c=random.randint(0,len(cadenas_dism)-1)
                cadena2=cadenas_dism[x_p_c]                        
                cadena=cadena+cadena2[-1::-1]
                cadenas_recomb.append(cadena)
                cadenas_dism.pop(x_p_c)
                break
            elif len(cadenas_trans)>0 and x_r_int==1:
                if len(cadenas_trans)==1:
                    x_p_c=0
                elif len(cadenas_trans)>1:
                    x_p_c=random.randint(0,len(cadenas_trans)-1) 
                    cadena2=cadenas_trans[x_p_c]                        
                    cadena=cadena+cadena2[-1::-1]
                    cadenas_recomb.append(cadena)
                    cadenas_trans.pop(x_p_c)
                    break 
            else:
                continue   
    elif CcccomboBreaker2==6: # Como recombinación agarra una cadena que ya esta terminada,
        cadenas_dism.append(cadena) # necesita cadenas en el array de dismutacion o transferencia. 
        cadena="I" # Cuando no hay cadenas en ninguno de las dos arrays, la cadena termina por dismutacion 
        i=i-1 # y se pone en el array de dism.
        CcccomboBreaker2=5
        CcccomboBreaker=5
        x_l=2
    elif CcccomboBreaker2==2: # Sirve para que las transferencias entren al ciclo de nuevo.
        CcccomboBreaker2=5
        if a>0 or b>0:        
            x_l=2 # Reseteo de variable para que no interfiera.
        elif a<=0 and b<=0:
            break
    else: # Es para las cadenas normales.
        i=i-1        
        cadena=I
        x_l=2 # Para que no haga interferencia y se resetee.
        if i<=-1 or (a<=0 and b<=0):
            break    
    while p_ia<=1: # Inicia la rutina Iniciador + Monómero.
        if a<=0 and b<=0 and cadena=="I":
            CcccomboBreaker2==3
            break
        if CcccomboBreaker==1: # Para salir de todos los ciclos cuando hay terminación.
            break
        else:
            x_ia=random.uniform(0,1) # Genera un float al azar entre 0 y 1.
            if (x_ia<=p_ia and a>0) or (a>0 and cadena[len(cadena)-1::]==A): # Solo si hay productos. Si el numero al va desde 0 hasta el valor de la probabilidad que I reaccione con A. O si la cadena por transferencia termina en A.
                if cadena=="I":
                    cadena=cadena + A # A la cadena 1 se le suma un A.
                    a=a-1
                if i<=-1 or (b>0 and a>0): # Para que no me de infinito.
                    p_aa=ra*a/(ra*a+b) # Probabilidad que reaccione A con A.
                    p_ab=1-p_aa
                    p_bb=rb*b/(rb*b+a)
                    p_ba=1-p_bb
                    p_ia=kia*a/(kia*a+kib*b) # Probabilidad que el iniciador reaccione con A
                    p_ib=1-p_ia
                x_l=random.uniform(0,1)
                while x_l<=1/Lambda and len(cadenas_trans)==0 and len(cadenas_dism)==0:
                        cadenas_dism.append(cadena)
                        CcccomboBreaker=1
                        break
                while x_l<=1/Lambda and len(cadena)==2: # Para que solo las cadenas nuevas estén sujetas a terminar, no las que acaban de reentrar.
                    if CcccomboBreaker==1:
                        break
                    while len(cadenas_trans)==0 and len(cadenas_dism)==0: # Como las otras terminaciones requieren de cadenas, cuando no hay cadenas para que agarren, la cadena termina por dismutación.
                        cadenas_dism.append(cadena)
                        CcccomboBreaker=1
                        break
                    x_t=random.uniform(0,1)
                    while x_t<=prob_dism: # Dismutación. 
                        cadenas_dism.append(cadena)
                        CcccomboBreaker=1
                        break
                    while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0): # Transferencia.
                        CcccomboBreaker=1  # Para salir de los ciclos y entrar al ciclo deseado para resetear sus valores.
                        CcccomboBreaker2=2 
                        x_r_int=random.randint(0,1)  # Volado para ver si se saca una cadena de las que terminaron por dismutación o transferencia.
                        if a<=0 and b<=0:
                            break
                        elif len(cadenas_dism)>0 and x_r_int==0:
                            if len(cadenas_dism)==1:
                                x_p_c=0
                            elif len(cadenas_dism)>1:
                                x_p_c=random.randint(0,len(cadenas_dism)-1)                            
                            cadenas_trans.append(cadena) # Poner la cadena que termina en el array de cadenas por transferencia.
                            cadena2=cadenas_dism[x_p_c] # Sacar la cadena correspondiente del arreglo correspondiente.
                            cadena=cadena2 # La cadena es ahora la cadena a la que se le transfirio el sitio activo.
                            cadenas_dism.pop(x_p_c) # Borrar la cadena de dicho arreglo.
                            break
                        elif len(cadenas_trans)>0 and x_r_int==1:
                            if len(cadenas_trans)==1:
                                x_p_c=0
                            elif len(cadenas_trans)>1:
                                x_p_c=random.randint(0,len(cadenas_trans)-1)
                            cadenas_trans.append(cadena)
                            cadena2=cadenas_trans[x_p_c]
                            cadena=cadena2
                            cadenas_trans.pop(x_p_c)
                            break
                        else:
                            continue
                    while x_t>(1-prob_recomb): #recombinacion
                        CcccomboBreaker=1
                        x_r_int=random.randint(0,1)
                        if a<=0 and b<=0:
                            break
                        if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                            CcccomboBreaker2=6
                            break
                        elif len(cadenas_dism)>0 and x_r_int==0:
                            if len(cadenas_dism)==1:
                                x_p_c=0
                            elif len(cadenas_dism)>1:
                                x_p_c=random.randint(0,len(cadenas_dism)-1) 
                            cadena2=cadenas_dism[x_p_c]
                            cadena=cadena+cadena2[-1::-1]
                            cadenas_recomb.append(cadena)
                            cadenas_dism.pop(x_p_c)
                            break
                        elif len(cadenas_trans)>0 and x_r_int==1:
                            if len(cadenas_trans)==1:
                                x_p_c=0
                            elif len(cadenas_dism)>1:
                                x_p_c=random.randint(0,len(cadenas_trans)-1) 
                            cadena2=cadenas_trans[x_p_c]                            
                            cadena=cadena+cadena2[-1::-1]
                            cadenas_recomb.append(cadena)
                            cadenas_trans.pop(x_p_c)
                            break 
                        else:
                            continue
                if x_l>1/Lambda or len(cadena)>2:                    
                    while p_aa<=1: # Inicia la subrutina para ver que reacciona con la cadena creciente.
                        if CcccomboBreaker==1:
                            break
                        else:                        
                            x_aa=random.uniform(0,1)
                            if x_aa<=p_aa and a>0: # Aqui comienza el programa cuando la cadena viene por transferencia.
                                cadena=cadena + A # A con A.
                                a=a-1
                                if b>0 and a>0:
                                    p_aa=ra*a/(ra*a+b)
                                    p_ab=1-p_aa
                                    p_bb=rb*b/(rb*b+a)
                                    p_ba=1-p_bb
                                    p_ia=kia*a/(kia*a+kib*b)
                                    p_ib=1-p_ia
                                x_l=random.uniform(0,1)
                                while x_l<=1/Lambda and len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                    cadenas_dism.append(cadena)
                                    CcccomboBreaker=1
                                    break
                                while x_l<=1/Lambda:
                                    if CcccomboBreaker==1:
                                        break
                                    x_t=random.uniform(0,1)
                                    while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0): 
                                        cadenas_dism.append(cadena)
                                        CcccomboBreaker=1
                                        break
                                    while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0): # Transferencia.
                                        CcccomboBreaker=1  
                                        CcccomboBreaker2=2 
                                        x_r_int=random.randint(0,1)    
                                        if a<=0 and b<=0:
                                            break
                                        elif len(cadenas_dism)>0 and x_r_int==0:
                                            if len(cadenas_dism)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_dism)-1)                            
                                            cadenas_trans.append(cadena)
                                            cadena2=cadenas_dism[x_p_c]                 
                                            cadena=cadena2
                                            cadenas_dism.pop(x_p_c)
                                            break
                                        elif len(cadenas_trans)>0 and x_r_int==1:
                                            if len(cadenas_trans)==1:
                                                x_p_c=0
                                            elif len(cadenas_trans)>1:
                                                x_p_c=random.randint(0,len(cadenas_trans)-1)
                                            cadenas_trans.append(cadena)
                                            cadena2=cadenas_trans[x_p_c]                                            
                                            cadena=cadena2
                                            cadenas_trans.pop(x_p_c)
                                            break
                                        else:
                                            continue
                                    while x_t>(1-prob_recomb): #recombinacion
                                        CcccomboBreaker=1
                                        x_r_int=random.randint(0,1)
                                        if a<=0 and b<=0:
                                            break
                                        if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                            CcccomboBreaker2=6
                                            break
                                        elif len(cadenas_dism)>0 and x_r_int==0:
                                            if len(cadenas_dism)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_dism)-1) 
                                            cadena2=cadenas_dism[x_p_c]                                            
                                            cadena=cadena+cadena2[-1::-1]
                                            cadenas_recomb.append(cadena)
                                            cadenas_dism.pop(x_p_c)
                                            break
                                        elif len(cadenas_trans)>0 and x_r_int==1:
                                            if len(cadenas_trans)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                            cadena2=cadenas_trans[x_p_c]                                            
                                            cadena=cadena+cadena2[-1::-1]
                                            cadenas_recomb.append(cadena)
                                            cadenas_trans.pop(x_p_c)
                                            break 
                                        else:
                                            continue
                                if x_l>1/Lambda:
                                    continue
                            elif x_aa>p_aa and b>0: # A con B.
                                cadena=cadena + B
                                b=b-1
                                if b>0 and a>0:
                                    p_aa=ra*a/(ra*a+b)
                                    p_ab=1-p_aa
                                    p_bb=rb*b/(rb*b+a)
                                    p_ba=1-p_bb
                                    p_ia=kia*a/(kia*a+kib*b)
                                    p_ib=1-p_ia
                                x_l=random.uniform(0,1)
                                while x_l<=1/Lambda and len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                    cadenas_dism.append(cadena)
                                    CcccomboBreaker=1
                                    break
                                while x_l<=1/Lambda:
                                    if CcccomboBreaker==1:
                                        break
                                    while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                        cadenas_dism.append(cadena)
                                        CcccomboBreaker=1
                                        break
                                    x_t=random.uniform(0,1)
                                    while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                                        cadenas_dism.append(cadena)
                                        CcccomboBreaker=1
                                        break
                                    while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0): # Transferencia.
                                        CcccomboBreaker=1
                                        CcccomboBreaker2=2
                                        x_r_int=random.randint(0,1)
                                        if a<=0 and b<=0:
                                            break
                                        elif len(cadenas_dism)>0 and x_r_int==0:
                                            if len(cadenas_dism)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_dism)-1)                            
                                            cadenas_trans.append(cadena) 
                                            cadena2=cadenas_dism[x_p_c] 
                                            cadenas_dism.pop(x_p_c) 
                                            cadena=cadena2 
                                            break
                                        elif len(cadenas_trans)>0 and x_r_int==1:
                                            if len(cadenas_trans)==1:
                                                x_p_c=0
                                            elif len(cadenas_trans)>1:
                                                x_p_c=random.randint(0,len(cadenas_trans)-1)
                                            cadenas_trans.append(cadena)
                                            cadena2=cadenas_trans[x_p_c]                                            
                                            cadena=cadena2
                                            cadenas_trans.pop(x_p_c)
                                            break
                                        else:
                                            continue
                                    while x_t>(1-prob_recomb): #recombinacion
                                        CcccomboBreaker=1
                                        x_r_int=random.randint(0,1) 
                                        if a<=0 and b<=0:
                                            break
                                        if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                            CcccomboBreaker2=6
                                            break
                                        elif len(cadenas_dism)>0 and x_r_int==0:
                                            if len(cadenas_dism)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_dism)-1) 
                                            cadena2=cadenas_dism[x_p_c]                                            
                                            cadena=cadena+cadena2[-1::-1]
                                            cadenas_recomb.append(cadena)
                                            cadenas_dism.pop(x_p_c)
                                            break
                                        elif len(cadenas_trans)>0 and x_r_int==1:
                                            if len(cadenas_trans)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                            cadena2=cadenas_trans[x_p_c]                                            
                                            cadena=cadena+cadena2[-1::-1]
                                            cadenas_recomb.append(cadena)
                                            cadenas_trans.pop(x_p_c)
                                            break
                                        else:
                                            continue
                                if x_l>1/Lambda:
                                    while p_bb<=1: # Inicia subrutina para ver que se pega con B.
                                        if CcccomboBreaker==1:
                                            break
                                        x_bb=random.uniform(0,1)
                                        if x_bb<=p_bb and b>0: # B con B.
                                            cadena=cadena + B
                                            b=b-1
                                            if b>0 and a>0:
                                                p_aa=ra*a/(ra*a+b)
                                                p_ab=1-p_aa
                                                p_bb=rb*b/(rb*b+a)
                                                p_ba=1-p_bb
                                                p_ia=kia*a/(kia*a+kib*b)
                                                p_ib=1-p_ia
                                            x_l=random.uniform(0,1)
                                            while x_l<=1/Lambda and len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                                cadenas_dism.append(cadena)
                                                CcccomboBreaker=1
                                                break
                                            while x_l<=1/Lambda:
                                                 if CcccomboBreaker==1:
                                                     break
                                                 while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                                     cadenas_dism.append(cadena)
                                                     CcccomboBreaker=1
                                                     break
                                                 x_t=random.uniform(0,1)
                                                 while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                                                     cadenas_dism.append(cadena)
                                                     CcccomboBreaker=1
                                                     break
                                                 while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0): # Transferencia.
                                                     CcccomboBreaker=1
                                                     CcccomboBreaker2=2 
                                                     x_r_int=random.randint(0,1)
                                                     if a<=0 and b<=0:
                                                         break
                                                     elif len(cadenas_dism)>0 and x_r_int==0:
                                                         if len(cadenas_dism)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_dism)-1)                            
                                                         cadenas_trans.append(cadena)
                                                         cadena2=cadenas_dism[x_p_c]                                                       
                                                         cadena=cadena2
                                                         cadenas_dism.pop(x_p_c)
                                                         break
                                                     elif len(cadenas_trans)>0 and x_r_int==1:
                                                         if len(cadenas_trans)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_trans)>1:
                                                             x_p_c=random.randint(0,len(cadenas_trans)-1)
                                                         cadenas_trans.append(cadena)
                                                         cadena2=cadenas_trans[x_p_c]                                                         
                                                         cadena=cadena2
                                                         cadenas_trans.pop(x_p_c)
                                                         break
                                                     else:
                                                         continue
                                                 while x_t>(1-prob_recomb): #recombinacion
                                                     CcccomboBreaker=1
                                                     x_r_int=random.randint(0,1)
                                                     if a<=0 and b<=0:
                                                         break
                                                     if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                                         CcccomboBreaker2=6
                                                         break
                                                     elif len(cadenas_dism)>0 and x_r_int==0:
                                                         if len(cadenas_dism)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_dism)-1) 
                                                         cadena2=cadenas_dism[x_p_c]                                                         
                                                         cadena=cadena+cadena2[-1::-1]
                                                         cadenas_recomb.append(cadena)
                                                         cadenas_dism.pop(x_p_c)
                                                         break
                                                     elif len(cadenas_trans)>0 and x_r_int==1:
                                                         if len(cadenas_trans)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                                         cadena2=cadenas_trans[x_p_c]                                                         
                                                         cadena=cadena+cadena2[-1::-1]
                                                         cadenas_recomb.append(cadena)
                                                         cadenas_trans.pop(x_p_c)
                                                         break
                                                     else:
                                                         continue
                                            if x_l>1/Lambda:
                                                continue
                                        elif x_bb>p_bb and a>0: # B con A.
                                            cadena=cadena + A
                                            a=a-1
                                            if b>0 and a>0:
                                                p_aa=ra*a/(ra*a+b)
                                                p_ab=1-p_aa
                                                p_bb=rb*b/(rb*b+a)
                                                p_ba=1-p_bb
                                                p_ia=kia*a/(kia*a+kib*b)
                                                p_ib=1-p_ia
                                            else:
                                                continue
                                            x_l=random.uniform(0,1)
                                            while x_l<=1/Lambda and len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                                cadenas_dism.append(cadena)
                                                CcccomboBreaker=1
                                                break
                                            while x_l<=1/Lambda:
                                                 if CcccomboBreaker==1:
                                                     break
                                                 while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                                     cadenas_dism.append(cadena)
                                                     CcccomboBreaker=1
                                                     break
                                                 x_t=random.uniform(0,1)
                                                 while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                                                     cadenas_dism.append(cadena)
                                                     CcccomboBreaker=1
                                                     break
                                                 while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0): # Transferencia.
                                                     CcccomboBreaker=1
                                                     CcccomboBreaker2=2 
                                                     x_r_int=random.randint(0,1)
                                                     if a<=0 and b<=0:
                                                         break
                                                     elif len(cadenas_dism)>0 and x_r_int==0:
                                                         if len(cadenas_dism)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_dism)-1)                            
                                                         cadenas_trans.append(cadena)
                                                         cadena2=cadenas_dism[x_p_c]                                                       
                                                         cadena=cadena2
                                                         cadenas_dism.pop(x_p_c)
                                                         break
                                                     elif len(cadenas_trans)>0 and x_r_int==1:
                                                         if len(cadenas_trans)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_trans)>1:
                                                             x_p_c=random.randint(0,len(cadenas_trans)-1)
                                                         cadenas_trans.append(cadena)
                                                         cadena2=cadenas_trans[x_p_c]                                                         
                                                         cadena=cadena2
                                                         cadenas_trans.pop(x_p_c)
                                                         break
                                                     else:
                                                         continue
                                                 while x_t>(1-prob_recomb): #recombinacion
                                                     CcccomboBreaker=1
                                                     x_r_int=random.randint(0,1)
                                                     if a<=0 and b<=0:
                                                         break
                                                     if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                                         CcccomboBreaker2=6
                                                         break
                                                     elif len(cadenas_dism)>0 and x_r_int==0:
                                                         if len(cadenas_dism)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_dism)-1) 
                                                         cadena2=cadenas_dism[x_p_c]                                                         
                                                         cadena=cadena+cadena2[-1::-1]
                                                         cadenas_recomb.append(cadena)
                                                         cadenas_dism.pop(x_p_c)
                                                         break
                                                     elif len(cadenas_trans)>0 and x_r_int==1:
                                                         if len(cadenas_trans)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                                         cadena2=cadenas_trans[x_p_c]                                                         
                                                         cadena=cadena+cadena2[-1::-1]
                                                         cadenas_recomb.append(cadena)
                                                         cadenas_trans.pop(x_p_c)
                                                         break
                                                     else:
                                                         continue
                                            if x_l>1/Lambda:
                                                break
                                        elif a<=0 and b<=0:
                                            if CcccomboBreaker==1:
                                                break
                                            while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                                cadenas_dism.append(cadena)
                                                CcccomboBreaker=1
                                                break
                                            x_t=random.uniform(0,1)
                                            while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                                                cadenas_dism.append(cadena)
                                                CcccomboBreaker=1
                                                break
                                            while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0):#transferencia
                                                CcccomboBreaker=1
                                                CcccomboBreaker2=3
                                                x_r_int=random.randint(0,1)                        
                                                if a<=0 and b<=0:
                                                    break
                                                elif len(cadenas_dism)>0 and x_r_int==0:
                                                     if len(cadenas_dism)==1:
                                                        x_p_c=0
                                                     elif len(cadenas_dism)>1:
                                                        x_p_c=random.randint(0,len(cadenas_dism)-1)                                                     
                                                     cadena2=cadenas_dism[x_p_c]                                                     
                                                     cadena=cadena2
                                                     cadenas_dism.append(cadena)
                                                     cadenas_dism.pop(x_p_c)
                                                     break
                                                elif len(cadenas_trans)>0 and x_r_int==1:
                                                     if len(cadenas_trans)==0:
                                                         break
                                                     elif len(cadenas_trans)==1:
                                                        x_p_c=0
                                                     elif len(cadenas_trans)>1:
                                                        x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                                     cadenas_trans.append(cadena)
                                                     cadena2=cadenas_trans[x_p_c]                                                     
                                                     cadena=cadena2
                                                     cadenas_trans.pop(x_p_c)
                                                     break
                                                else:
                                                    continue
                                            while x_t>(1-prob_recomb): #recombinacion
                                                CcccomboBreaker=1
                                                x_r_int=random.randint(0,1)
                                                if a<=0 and b<=0:
                                                    break
                                                if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                                    cadenas_dism.append(cadena)
                                                    CcccomboBreaker2=3
                                                    break
                                                elif len(cadenas_dism)>0 and x_r_int==0:
                                                    if len(cadenas_dism)==1:
                                                        x_p_c=0
                                                    elif len(cadenas_dism)>1:
                                                        x_p_c=random.randint(0,len(cadenas_dism)-1)
                                                    cadena2=cadenas_dism[x_p_c]                                                    
                                                    cadena=cadena+cadena2[-1::-1]
                                                    cadenas_recomb.append(cadena)
                                                    cadenas_dism.pop(x_p_c)
                                                    break
                                                elif len(cadenas_trans)>0 and x_r_int==1:
                                                    if len(cadenas_trans)==0:
                                                        break
                                                    elif len(cadenas_trans)==1:
                                                        x_p_c=0
                                                    elif len(cadenas_trans)>1:
                                                        x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                                    cadena2=cadenas_trans[x_p_c]                                                    
                                                    cadena=cadena+cadena2[-1::-1]
                                                    cadenas_recomb.append(cadena)
                                                    cadenas_trans.pop(x_p_c)
                                                    break
                                                else:
                                                    continue
                            elif a<=0 and b<=0:
                                if CcccomboBreaker==1:
                                    break
                                while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                    cadenas_dism.append(cadena)
                                    CcccomboBreaker=1
                                    break
                                x_t=random.uniform(0,1)
                                while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                                    cadenas_dism.append(cadena)
                                    CcccomboBreaker=1
                                    break
                                while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0):#transferencia
                                    CcccomboBreaker=1
                                    x_r_int=random.randint(0,1)                        
                                    if a<=0 and b<=0:
                                        break
                                    elif len(cadenas_dism)>0 and x_r_int==0:
                                        if len(cadenas_dism)==1:
                                            x_p_c=0
                                        elif len(cadenas_dism)>1:
                                            x_p_c=random.randint(0,len(cadenas_dism)-1)
                                        cadenas_trans.append(cadena)
                                        cadena2=cadenas_dism[x_p_c]                                        
                                        cadena=cadena2
                                        cadenas_dism.pop(x_p_c)
                                        break
                                    elif len(cadenas_trans)>0 and x_r_int==1:
                                        if len(cadenas_trans)==1:
                                            x_p_c=0
                                        elif len(cadenas_trans)>1:
                                            x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                        cadenas_trans.append(cadena)
                                        cadena2=cadenas_trans[x_p_c]                                        
                                        cadena=cadena2
                                        cadenas_trans.pop(x_p_c)
                                        break
                                    else:
                                        continue
                                while x_t>(1-prob_recomb): #recombinacion
                                    CcccomboBreaker=1
                                    x_r_int=random.randint(0,1) 
                                    if a<=0 and b<=0:
                                        break
                                    if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                        cadenas_dism.append(cadena)
                                        CcccomboBreaker2=3
                                        break
                                    elif len(cadenas_dism)>0 and x_r_int==0:
                                        if len(cadenas_dism)==1:
                                            x_p_c=0
                                        elif len(cadenas_dism)>1:
                                            x_p_c=random.randint(0,len(cadenas_dism)-1)
                                        cadena2=cadenas_dism[x_p_c]                                        
                                        cadena=cadena+cadena2[-1::-1]
                                        cadenas_recomb.append(cadena)
                                        cadenas_dism.pop(x_p_c)
                                        break
                                    elif len(cadenas_trans)>0 and x_r_int==1:
                                        if len(cadenas_trans)==1:
                                            x_p_c=0
                                        elif len(cadenas_trans)>1:
                                            x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                        cadena2=cadenas_trans[x_p_c]                                        
                                        cadena=cadena+cadena2[-1::-1]
                                        cadenas_recomb.append(cadena)
                                        cadenas_trans.pop(x_p_c)
                                        break
                                    else:
                                        continue
            elif (x_ia>p_ia and b>0) or (b>0 and cadena[len(cadena)-1::]==B): # Las mismas anotaciones aplican para este lado del ciclo, solo que se comienza añadiendo B.
                if cadena=="I":
                    cadena=cadena + B
                    b=b-1               
                if b>0 and a>0:
                    p_aa=ra*a/(ra*a+b)
                    p_ab=1-p_aa
                    p_bb=rb*b/(rb*b+a)
                    p_ba=1-p_bb
                    p_ia=kia*a/(kia*a+kib*b)
                    p_ib=1-p_ia
                x_l=random.uniform(0,1)
                while x_l<=1/Lambda and len(cadenas_trans)==0 and len(cadenas_dism)==0:
                    cadenas_dism.append(cadena)
                    CcccomboBreaker=1
                    break
                while x_l<=1/Lambda and len(cadena)==2:
                    if CcccomboBreaker==1:
                        break
                    while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                        cadenas_dism.append(cadena)
                        CcccomboBreaker=1
                        break
                    x_t=random.uniform(0,1)
                    while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0): # Dismutación.
                        cadenas_dism.append(cadena)
                        CcccomboBreaker=1
                        break
                    while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0): # Transferencia.
                        CcccomboBreaker=1
                        CcccomboBreaker2=2 
                        x_r_int=random.randint(0,1)
                        if a<=0 and b<=0:
                            break
                        elif len(cadenas_dism)>0 and x_r_int==0:
                            if len(cadenas_dism)==1:
                                x_p_c=0
                            elif len(cadenas_dism)>1:
                                x_p_c=random.randint(0,len(cadenas_dism)-1)                            
                            cadenas_trans.append(cadena)
                            cadena2=cadenas_dism[x_p_c]                          
                            cadena=cadena2
                            cadenas_dism.pop(x_p_c)
                            break
                        elif len(cadenas_trans)>0 and x_r_int==1:
                            if len(cadenas_trans)==1:
                                x_p_c=0
                            elif len(cadenas_trans)>1:
                                x_p_c=random.randint(0,len(cadenas_trans)-1)
                            cadenas_trans.append(cadena)
                            cadena2=cadenas_trans[x_p_c]
                            cadenas_trans.pop(x_p_c)
                            cadena=cadena2                            
                            break
                        else:
                            continue
                    while x_t>(1-prob_recomb): #recombinacion
                        CcccomboBreaker=1
                        x_r_int=random.randint(0,1)
                        if a<=0 and b<=0:
                            break
                        if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                            CcccomboBreaker2=6
                            break                       
                        elif len(cadenas_dism)>0 and x_r_int==0:
                            if len(cadenas_dism)==1:
                                x_p_c=0
                            elif len(cadenas_dism)>1:
                                x_p_c=random.randint(0,len(cadenas_dism)-1) 
                            cadena2=cadenas_dism[x_p_c]                            
                            cadena=cadena+cadena2[-1::-1]
                            cadenas_recomb.append(cadena)
                            cadenas_dism.pop(x_p_c)
                            break
                        elif len(cadenas_trans)>0 and x_r_int==1:
                            if len(cadenas_trans)==1:
                                x_p_c=0
                            elif len(cadenas_dism)>1:
                                x_p_c=random.randint(0,len(cadenas_trans)-1) 
                            cadena2=cadenas_trans[x_p_c]                            
                            cadena=cadena+cadena2[-1::-1]
                            cadenas_recomb.append(cadena)
                            cadenas_trans.pop(x_p_c)
                            break
                        else:
                            continue
                if x_l>1/Lambda or len(cadena)>2:
                    while p_bb<=1:
                        if CcccomboBreaker==1:
                            break
                        else:                        
                            x_bb=random.uniform(0,1)
                            if x_bb<=p_bb and b>0:
                                cadena=cadena + B
                                b=b-1
                                if b>0 and a>0:
                                    p_aa=ra*a/(ra*a+b)
                                    p_ab=1-p_aa
                                    p_bb=rb*b/(rb*b+a)
                                    p_ba=1-p_bb
                                    p_ia=kia*a/(kia*a+kib*b)
                                    p_ib=1-p_ia
                                x_l=random.uniform(0,1)
                                while x_l<=1/Lambda and len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                    cadenas_dism.append(cadena)
                                    CcccomboBreaker=1
                                    break
                                while x_l<=1/Lambda:
                                    if CcccomboBreaker==1:
                                        break
                                    while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                        cadenas_dism.append(cadena)
                                        CcccomboBreaker=1
                                        break
                                    x_t=random.uniform(0,1)
                                    while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                                        cadenas_dism.append(cadena)
                                        CcccomboBreaker=1
                                        break
                                    while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0): # Transferencia.
                                        CcccomboBreaker=1
                                        CcccomboBreaker2=2 
                                        x_r_int=random.randint(0,1)
                                        if a<=0 and b<=0:
                                            break
                                        elif len(cadenas_dism)>0 and x_r_int==0:
                                            if len(cadenas_dism)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_dism)-1)                            
                                            cadenas_trans.append(cadena)
                                            cadena2=cadenas_dism[x_p_c]                                          
                                            cadena=cadena2
                                            cadenas_dism.pop(x_p_c)
                                            break
                                        elif len(cadenas_trans)>0 and x_r_int==1:
                                            if len(cadenas_trans)==1:
                                                x_p_c=0
                                            elif len(cadenas_trans)>1:
                                                x_p_c=random.randint(0,len(cadenas_trans)-1)
                                            cadenas_trans.append(cadena)
                                            cadena2=cadenas_trans[x_p_c]                                            
                                            cadena=cadena2
                                            cadenas_trans.pop(x_p_c)
                                            break
                                        else:
                                            continue
                                    while x_t>(1-prob_recomb): #recombinacion
                                        CcccomboBreaker=1
                                        x_r_int=random.randint(0,1)
                                        if a<=0 and b<=0:
                                            break
                                        if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                            CcccomboBreaker2=6
                                            break
                                        elif len(cadenas_dism)>0 and x_r_int==0:
                                            if len(cadenas_dism)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_dism)-1) 
                                            cadena2=cadenas_dism[x_p_c]                                            
                                            cadena=cadena+cadena2[-1::-1]
                                            cadenas_recomb.append(cadena)
                                            cadenas_dism.pop(x_p_c)
                                            break
                                        elif len(cadenas_trans)>0 and x_r_int==1:
                                            if len(cadenas_trans)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                            cadena2=cadenas_trans[x_p_c]                                            
                                            cadena=cadena+cadena2[-1::-1]
                                            cadenas_recomb.append(cadena)
                                            cadenas_trans.pop(x_p_c)
                                            break
                                        else:
                                            continue
                                if x_l>1/Lambda:
                                    continue
                            elif x_bb>p_bb and a>0:
                                cadena=cadena + A
                                a=a-1
                                if b>0 and a>0:
                                    p_aa=ra*a/(ra*a+b)
                                    p_ab=1-p_aa
                                    p_bb=rb*b/(rb*b+a)
                                    p_ba=1-p_bb
                                    p_ia=kia*a/(kia*a+kib*b)
                                    p_ib=1-p_ia
                                x_l=random.uniform(0,1)
                                while x_l<=1/Lambda and len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                    cadenas_dism.append(cadena)
                                    CcccomboBreaker=1
                                    break
                                while x_l<=1/Lambda:
                                    if CcccomboBreaker==1:
                                        break
                                    while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                        cadenas_dism.append(cadena)
                                        CcccomboBreaker=1
                                        break
                                    x_t=random.uniform(0,1)
                                    while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                                        cadenas_dism.append(cadena)
                                        CcccomboBreaker=1
                                        break
                                    while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0): # Transferencia.
                                        CcccomboBreaker=1
                                        CcccomboBreaker2=2 
                                        x_r_int=random.randint(0,1)
                                        if a<=0 and b<=0:
                                            break
                                        elif len(cadenas_dism)>0 and x_r_int==0:
                                            if len(cadenas_dism)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_dism)-1)                            
                                            cadenas_trans.append(cadena)
                                            cadena2=cadenas_dism[x_p_c]                                          
                                            cadena=cadena2
                                            cadenas_dism.pop(x_p_c)
                                            break
                                        elif len(cadenas_trans)>0 and x_r_int==1:
                                            if len(cadenas_trans)==1:
                                                x_p_c=0
                                            elif len(cadenas_trans)>1:
                                                x_p_c=random.randint(0,len(cadenas_trans)-1)
                                            cadenas_trans.append(cadena)
                                            cadena2=cadenas_trans[x_p_c]                                            
                                            cadena=cadena2
                                            cadenas_trans.pop(x_p_c)
                                            break
                                        else:
                                            continue
                                    while x_t>(1-prob_recomb): #recombinacion
                                        CcccomboBreaker=1
                                        x_r_int=random.randint(0,1) 
                                        if a<=0 and b<=0:
                                            break
                                        if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                            CcccomboBreaker2=6
                                            break
                                        elif len(cadenas_dism)>0 and x_r_int==0:
                                            if len(cadenas_dism)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_dism)-1) 
                                            cadena2=cadenas_dism[x_p_c]                                            
                                            cadena=cadena+cadena2[-1::-1]
                                            cadenas_recomb.append(cadena)
                                            cadenas_dism.pop(x_p_c)
                                            break
                                        elif len(cadenas_trans)>0 and x_r_int==1:
                                            if len(cadenas_trans)==1:
                                                x_p_c=0
                                            elif len(cadenas_dism)>1:
                                                x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                            cadena2=cadenas_trans[x_p_c]                                            
                                            cadena=cadena+cadena2[-1::-1]
                                            cadenas_recomb.append(cadena)
                                            cadenas_trans.pop(x_p_c)
                                            break
                                        else:
                                            continue
                                if x_l>1/Lambda:
                                    while p_aa<=1:
                                        if CcccomboBreaker==1:
                                            break
                                        x_aa=random.uniform(0,1)
                                        if x_aa<=p_aa and a>0:
                                            cadena=cadena + A
                                            a=a-1
                                            if b>0 and a>0:
                                                p_aa=ra*a/(ra*a+b)
                                                p_ab=1-p_aa
                                                p_bb=rb*b/(rb*b+a)
                                                p_ba=1-p_bb
                                                p_ia=kia*a/(kia*a+kib*b)
                                                p_ib=1-p_ia
                                            x_l=random.uniform(0,1)
                                            while x_l<=1/Lambda and len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                                cadenas_dism.append(cadena)
                                                CcccomboBreaker=1
                                                break
                                            while x_l<=1/Lambda:
                                                 if CcccomboBreaker==1:
                                                     break
                                                 while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                                     cadenas_dism.append(cadena)
                                                     CcccomboBreaker=1
                                                     break
                                                 x_t=random.uniform(0,1)
                                                 while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                                                     cadenas_dism.append(cadena)
                                                     CcccomboBreaker=1
                                                     break
                                                 while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0):
                                                     CcccomboBreaker=1
                                                     CcccomboBreaker2=2
                                                     x_r_int=random.randint(0,1)
                                                     if a<=0 and b<=0:
                                                         break
                                                     elif len(cadenas_dism)>0 and x_r_int==0:
                                                         if len(cadenas_dism)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_dism)-1)                            
                                                         cadenas_trans.append(cadena)
                                                         cadena2=cadenas_dism[x_p_c]                                                        
                                                         cadena=cadena2
                                                         cadenas_dism.pop(x_p_c)
                                                         break
                                                     elif len(cadenas_trans)>0 and x_r_int==1:
                                                         if len(cadenas_trans)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_trans)>1:
                                                             x_p_c=random.randint(0,len(cadenas_trans)-1)
                                                         cadenas_trans.append(cadena)
                                                         cadena2=cadenas_trans[x_p_c]                                                         
                                                         cadena=cadena2
                                                         cadenas_trans.pop(x_p_c)
                                                         break
                                                     else:
                                                         continue
                                                 while x_t>(1-prob_recomb): #recombinacion
                                                     CcccomboBreaker=1
                                                     x_r_int=random.randint(0,1) 
                                                     if a<=0 and b<=0:
                                                         break
                                                     if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                                         CcccomboBreaker2=6
                                                         break
                                                     elif len(cadenas_dism)>0 and x_r_int==0:
                                                         if len(cadenas_dism)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_dism)-1) 
                                                         cadena2=cadenas_dism[x_p_c]                                                         
                                                         cadena=cadena+cadena2[-1::-1]
                                                         cadenas_recomb.append(cadena)
                                                         cadenas_dism.pop(x_p_c)
                                                         break
                                                     elif len(cadenas_trans)>0 and x_r_int==1:
                                                         if len(cadenas_trans)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                                         cadena2=cadenas_trans[x_p_c]                                                         
                                                         cadena=cadena+cadena2[-1::-1]
                                                         cadenas_recomb.append(cadena)
                                                         cadenas_trans.pop(x_p_c)
                                                         break
                                                     else:
                                                         continue
                                            if x_l>1/Lambda:
                                                continue
                                        elif x_aa>p_aa and b>0:
                                            cadena=cadena + B
                                            b=b-1
                                            if b>0 and a>0:
                                                p_aa=ra*a/(ra*a+b)
                                                p_ab=1-p_aa
                                                p_bb=rb*b/(rb*b+a)
                                                p_ba=1-p_bb
                                                p_ia=kia*a/(kia*a+kib*b)
                                                p_ib=1-p_ia
                                            x_l=random.uniform(0,1)
                                            while x_l<=1/Lambda and len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                                cadenas_dism.append(cadena)
                                                CcccomboBreaker=1
                                                break
                                            while x_l<=1/Lambda:
                                                 if CcccomboBreaker==1:
                                                     break
                                                 while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                                     cadenas_dism.append(cadena)
                                                     CcccomboBreaker=1
                                                     break
                                                 x_t=random.uniform(0,1)
                                                 while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0): # Dismutación
                                                     cadenas_dism.append(cadena)
                                                     CcccomboBreaker=1
                                                     break
                                                 while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0): # Transferencia.
                                                     CcccomboBreaker=1
                                                     CcccomboBreaker2=2 
                                                     x_r_int=random.randint(0,1)
                                                     if a<=0 and b<=0:
                                                         break
                                                     elif len(cadenas_dism)>0 and x_r_int==0:
                                                         if len(cadenas_dism)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_dism)-1)                            
                                                         cadenas_trans.append(cadena)
                                                         cadena2=cadenas_dism[x_p_c]                                                        
                                                         cadena=cadena2
                                                         cadenas_dism.pop(x_p_c)
                                                         break
                                                     elif len(cadenas_trans)>0 and x_r_int==1:
                                                         if len(cadenas_trans)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_trans)>1:
                                                             x_p_c=random.randint(0,len(cadenas_trans)-1)
                                                         cadenas_trans.append(cadena)
                                                         cadena2=cadenas_trans[x_p_c]                                                         
                                                         cadena=cadena2
                                                         cadenas_trans.pop(x_p_c)
                                                         break
                                                     else:
                                                         continue
                                                 while x_t>(1-prob_recomb): #recombinacion
                                                     CcccomboBreaker=1
                                                     x_r_int=random.randint(0,1)
                                                     if a<=0 and b<=0:
                                                         break
                                                     if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                                         CcccomboBreaker2=6
                                                         break
                                                     elif len(cadenas_dism)>0 and x_r_int==0:
                                                         if len(cadenas_dism)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_dism)-1) 
                                                         cadena2=cadenas_dism[x_p_c]                                                         
                                                         cadena=cadena+cadena2[-1::-1]
                                                         cadenas_recomb.append(cadena)
                                                         cadenas_dism.pop(x_p_c)
                                                         break
                                                     elif len(cadenas_trans)>0 and x_r_int==1:
                                                         if len(cadenas_trans)==1:
                                                             x_p_c=0
                                                         elif len(cadenas_dism)>1:
                                                             x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                                         cadena2=cadenas_trans[x_p_c]                                                         
                                                         cadena=cadena+cadena2[-1::-1]
                                                         cadenas_recomb.append(cadena)
                                                         cadenas_trans.pop(x_p_c)
                                                         break
                                                     else:
                                                         continue
                                            if x_l>1/Lambda:
                                                break
                                        elif a<=0 and b<=0:
                                            if CcccomboBreaker==1:
                                                break
                                            while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                                cadenas_dism.append(cadena)
                                                CcccomboBreaker=1
                                                break
                                            x_t=random.uniform(0,1)
                                            while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                                                cadenas_dism.append(cadena)
                                                CcccomboBreaker=1
                                                break
                                            while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0):#transferencia
                                                CcccomboBreaker=1
                                                x_r_int=random.randint(0,1)                        
                                                if a<=0 and b<=0:
                                                    break
                                                elif len(cadenas_dism)>0 and x_r_int==0:
                                                     if len(cadenas_dism)==1:
                                                        x_p_c=0
                                                     elif len(cadenas_dism)>1:
                                                        x_p_c=random.randint(0,len(cadenas_dism)-1)                                                     
                                                     cadena2=cadenas_dism[x_p_c]                                                     
                                                     cadena=cadena2
                                                     cadenas_dism.append(cadena)
                                                     cadenas_dism.pop(x_p_c)
                                                     break
                                                elif len(cadenas_trans)>0 and x_r_int==1:
                                                     if len(cadenas_trans)==0:
                                                         break
                                                     elif len(cadenas_trans)==1:
                                                        x_p_c=0
                                                     elif len(cadenas_trans)>1:
                                                        x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                                     cadenas_trans.append(cadena)
                                                     cadena2=cadenas_trans[x_p_c]                                                     
                                                     cadena=cadena2
                                                     cadenas_trans.pop(x_p_c)
                                                     break
                                                else:
                                                    continue
                                            while x_t>(1-prob_recomb): #recombinacion
                                                CcccomboBreaker=1
                                                x_r_int=random.randint(0,1) 
                                                if a<=0 and b<=0:
                                                    break
                                                if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                                    cadenas_dism.append(cadena)
                                                    break
                                                elif len(cadenas_dism)>0 and x_r_int==0:
                                                    if len(cadenas_dism)==1:
                                                        x_p_c=0
                                                    elif len(cadenas_dism)>1:
                                                        x_p_c=random.randint(0,len(cadenas_dism)-1)
                                                    cadena2=cadenas_dism[x_p_c]                                                    
                                                    cadena=cadena+cadena2[-1::-1]
                                                    cadenas_recomb.append(cadena)
                                                    cadenas_dism.pop(x_p_c)
                                                    break
                                                elif len(cadenas_trans)>0 and x_r_int==1:
                                                    if len(cadenas_trans)==0:
                                                        break
                                                    elif len(cadenas_trans)==1:
                                                        x_p_c=0
                                                    elif len(cadenas_trans)>1:
                                                        x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                                    cadena2=cadenas_trans[x_p_c]                                                    
                                                    cadena=cadena+cadena2[-1::-1]
                                                    cadenas_recomb.append(cadena)
                                                    cadenas_trans.pop(x_p_c)
                                                    break
                                                else:
                                                    continue
                            elif a<=0 and b<=0:
                                if CcccomboBreaker==1:
                                    break
                                while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                                    cadenas_dism.append(cadena)
                                    CcccomboBreaker=1
                                    break
                                x_t=random.uniform(0,1)
                                while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                                    cadenas_dism.append(cadena)
                                    CcccomboBreaker=1
                                    break
                                while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0):#transferencia
                                    CcccomboBreaker=1
                                    x_r_int=random.randint(0,1)                        
                                    if a<=0 and b<=0:
                                        break
                                    elif len(cadenas_dism)>0 and x_r_int==0:
                                        if len(cadenas_dism)==1:
                                            x_p_c=0
                                        elif len(cadenas_dism)>1:
                                            x_p_c=random.randint(0,len(cadenas_dism)-1)
                                        cadenas_trans.append(cadena)
                                        cadena2=cadenas_dism[x_p_c]                                        
                                        cadena=cadena2
                                        cadenas_dism.pop(x_p_c)
                                        break
                                    elif len(cadenas_trans)>0 and x_r_int==1:
                                        if len(cadenas_trans)==1:
                                            x_p_c=0
                                        elif len(cadenas_trans)>1:
                                            x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                        cadenas_trans.append(cadena)
                                        cadena2=cadenas_trans[x_p_c]                                        
                                        cadena=cadena2
                                        cadenas_trans.pop(x_p_c)
                                        break
                                    else:
                                        continue
                                while x_t>(1-prob_recomb): #recombinacion
                                    CcccomboBreaker=1
                                    x_r_int=random.randint(0,1)  
                                    if a<=0 and b<=0:
                                        break
                                    if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                                        cadenas_dism.append(cadena)
                                        break
                                    elif len(cadenas_dism)>0 and x_r_int==0:
                                        if len(cadenas_dism)==1:
                                            x_p_c=0
                                        elif len(cadenas_dism)>1:
                                            x_p_c=random.randint(0,len(cadenas_dism)-1)
                                        cadena2=cadenas_dism[x_p_c]                                        
                                        cadena=cadena+cadena2[-1::-1]
                                        cadenas_recomb.append(cadena)
                                        cadenas_dism.pop(x_p_c)
                                        break
                                    elif len(cadenas_trans)>0 and x_r_int==1:
                                        if len(cadenas_trans)==1:
                                            x_p_c=0
                                        elif len(cadenas_trans)>1:
                                            x_p_c=random.randint(0,len(cadenas_trans)-1) 
                                        cadena2=cadenas_trans[x_p_c]                                        
                                        cadena=cadena+cadena2[-1::-1]
                                        cadenas_recomb.append(cadena)
                                        cadenas_trans.pop(x_p_c)
                                        break
                                    else:
                                        continue
            elif a<=0 and b<=0:
                if CcccomboBreaker==1:
                    break
                while len(cadenas_trans)==0 and len(cadenas_dism)==0:
                    cadenas_dism.append(cadena)
                    CcccomboBreaker=1
                    break
                x_t=random.uniform(0,1)
                while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
                    cadenas_dism.append(cadena)
                    CcccomboBreaker=1
                    break
                while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0):#transferencia
                    CcccomboBreaker=1
                    x_r_int=random.randint(0,1)                        
                    if a<=0 and b<=0:
                        break
                    elif len(cadenas_dism)>0 and x_r_int==0:
                        if len(cadenas_dism)==1:
                            x_p_c=0
                        elif len(cadenas_dism)>1:
                            x_p_c=random.randint(0,len(cadenas_dism)-1)
                        cadenas_trans.append(cadena)
                        cadena2=cadenas_dism[x_p_c]                        
                        cadena=cadena2+"QSY" 
                        cadenas_trans.append(cadena)
                        cadenas_dism.pop(x_p_c)
                        break
                    elif len(cadenas_trans)>0 and x_r_int==1:
                        if len(cadenas_trans)==1:
                            x_p_c=0
                        elif len(cadenas_trans)>1:
                            x_p_c=random.randint(0,len(cadenas_trans)-1) 
                        cadenas_trans.append(cadena)
                        cadena2=cadenas_trans[x_p_c]                        
                        cadena=cadena2+"QSY"
                        cadenas_trans.append(cadena)
                        cadenas_trans.pop(x_p_c)
                        break
                    else:
                        continue
                while x_t>(1-prob_recomb): #recombinacion
                    CcccomboBreaker=1
                    x_r_int=random.randint(0,1) 
                    if a<=0 and b<=0:
                        break
                    if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                        cadenas_dism.append(cadena)
                        break                       
                    elif len(cadenas_dism)>0 and x_r_int==0:
                        if len(cadenas_dism)==1:
                            x_p_c=0
                        elif len(cadenas_dism)>1:
                            x_p_c=random.randint(0,len(cadenas_dism)-1)
                        cadena2=cadenas_dism[x_p_c]                        
                        cadena=cadena+cadena2[-1::-1]
                        cadenas_recomb.append(cadena)
                        cadenas_dism.pop(x_p_c)
                        break
                    elif len(cadenas_trans)>0 and x_r_int==1:
                        if len(cadenas_trans)==1:
                            x_p_c=0
                        elif len(cadenas_trans)>1:
                            x_p_c=random.randint(0,len(cadenas_trans)-1) 
                        cadena2=cadenas_trans[x_p_c]                        
                        cadena=cadena+cadena2[-1::-1]
                        cadenas_recomb.append(cadena)
                        cadenas_trans.pop(x_p_c)
                        break 
                    else:
                        continue
                    
# Función "Atrapa mierdas". Esta función toma la última cadena que a veces se queda sin anexar, 
# y si no se encuentra dentro de los arrays del polímero. La pone en el polímero de acuerdo 
# a las probabilidades de terminación.
for i in range(0,1):
    if cadena in cadenas_dism or cadena in cadenas_trans or cadena in cadenas_recomb:
        break
    elif cadena==I:
        break
    else:
        while len(cadenas_trans)==0 and len(cadenas_dism)==0: # dismutación
            cadenas_dism.append(cadena)
            CcccomboBreaker=2
            break
        x_t=random.uniform(0,1)
        while x_t<=prob_dism or (len(cadenas_trans)==0 and len(cadenas_dism)==0):
            cadenas_dism.append(cadena)
            CcccomboBreaker=2
            break
        while x_t>prob_dism and x_t<=(1-prob_recomb) and (len(cadenas_dism)>0 or len(cadenas_trans)>0):#transferencia
            CcccomboBreaker=2
            CcccomboBreaker2=8
            x_r_int=random.randint(0,1)                        
            if a<=0 and b<=0:
                cadenas_trans.append(cadena)
                CcccomboBreaker2=3
                break
            elif len(cadenas_dism)>0 and x_r_int==0:
                if len(cadenas_dism)==1:
                    x_p_c=0
                elif len(cadenas_dism)>1:
                    x_p_c=random.randint(0,len(cadenas_dism)-1)
                cadenas_trans.append(cadena)
                cadena2=cadenas_dism[x_p_c]                        
                cadena=cadena2+"QSY" 
                cadenas_trans.append(cadena)
                cadenas_dism.pop(x_p_c)
                break
            elif len(cadenas_trans)>0 and x_r_int==1:
                if len(cadenas_trans)==1:
                    x_p_c=0
                elif len(cadenas_trans)>1:
                    x_p_c=random.randint(0,len(cadenas_trans)-1) 
                    cadenas_trans.append(cadena)
                    cadena2=cadenas_trans[x_p_c]                        
                    cadena=cadena2+"QSY"
                    cadenas_trans.append(cadena)
                    cadenas_trans.pop(x_p_c)
                    break
                else:
                    continue
        while x_t>(1-prob_recomb): #recombinacion
            CcccomboBreaker=2
            CcccomboBreaker2=3 
            x_r_int=random.randint(0,1) 
            if len(cadenas_dism)==0 and len(cadenas_recomb)==0:
                cadenas_dism.append(cadena) 
                break                     
            elif len(cadenas_dism)>0 and x_r_int==0:
                if len(cadenas_dism)==1:
                    x_p_c=0
                elif len(cadenas_dism)>1:
                    x_p_c=random.randint(0,len(cadenas_dism)-1)
                cadena2=cadenas_dism[x_p_c]                        
                cadena=cadena+cadena2[-1::-1]
                cadenas_recomb.append(cadena)
                cadenas_dism.pop(x_p_c)
                break
            elif len(cadenas_trans)>0 and x_r_int==1:
                if len(cadenas_trans)==1:
                    x_p_c=0
                elif len(cadenas_trans)>1:
                    x_p_c=random.randint(0,len(cadenas_trans)-1) 
                    cadena2=cadenas_trans[x_p_c]                        
                    cadena=cadena+cadena2[-1::-1]
                    cadenas_recomb.append(cadena)
                    cadenas_trans.pop(x_p_c)
                    break 
            else:
                continue


# Estadísticas.     
print("__________________________\n\n")

# Dismutación.     
print("Estadísticas.\n")
print("__________________________\n")
print("Dismutación")
if len(cadenas_dism)>0: # Para que solo si hay cadenas por dismutación. Hay uno para cada tipo de terminación.
    peso_d_I=[] # Se van a ir guardando los pesos.
    peso_d_A=[]
    peso_d_B=[]
    peso_cadenas_d=[]
    peso_total_d=float()
    peso_total_d_Xp=float()
    num_d_I=int()
    num_d_A=int()
    num_d_B=int()
    for i in range(len(cadenas_dism)): # Va a iterar a lo largo del arreglo.
        peso_d_i=cadenas_dism[i].count(I)*MM_I # Variable temporal para guardar el peso hasta el proximo ciclo.
        peso_d_a=cadenas_dism[i].count(A)*MM_A
        peso_d_b=cadenas_dism[i].count(B)*MM_B
        peso_cadenas_d.append(peso_d_a+peso_d_b+peso_d_i) # Suma los pesos de los componentes y los pone en un nuev arreglo.
        peso_total_d=peso_total_d+peso_cadenas_d[i] # Para calcular Xn. Va sumando los pesos totales con cada iteración. 
        peso_total_d_Xp=peso_total_d_Xp+peso_cadenas_d[i]**2 # Para calcular Xp. Va sumando los pesos totales cuadrados con cada iteración.
        peso_d_I.append(cadenas_dism[i].count(I)*MM_I) # Guarda el peso de cada componente en un array. Por si se quiere usar posteriormente.
        peso_d_A.append(cadenas_dism[i].count(A)*MM_A)
        peso_d_B.append(cadenas_dism[i].count(B)*MM_B)
        num_d_I=num_d_I+cadenas_dism[i].count(I)
        num_d_A=num_d_A+cadenas_dism[i].count(A)
        num_d_B=num_d_B+cadenas_dism[i].count(B)
    Xn_d=peso_total_d/len(peso_cadenas_d) # Cáclulo de Xn.
    Xp_d=peso_total_d_Xp/(peso_total_d)
    polimol_d=Xp_d/Xn_d
    Lambda_d=(num_d_A+num_d_B)/len(cadenas_dism)
    print("pesos A = "+str(peso_d_A))
    print("pesos B = "+str(peso_d_B))
    print("pesos I = "+str(peso_d_I))
    print("pesos = "+str(peso_cadenas_d))
    print("peso total = "+str(peso_total_d))
    print("peso total Xp= "+str(peso_total_d_Xp))
    print("Xn Dismutación = "+str(Xn_d))
    print("Xp Dismutación = "+str(Xp_d))
    print("Polimolecularidad Dismutación = "+str(polimol_d))
    print("Lambda Dismutación = "+str(Lambda_d))
    print("\n")
else:
    print("No se toma en cuenta porque no hay.\n") # En caso de que no haya. Así no arroja error.
    peso_d_I=0
    peso_d_A=0
    peso_d_B=0
    peso_total_d=0
    peso_total_d_Xp=0
    num_d_I=0
    num_d_A=0
    num_d_B=0
    Xn_d=0
    Xp_d=0
    polimol_d=0
    Lambda_d=0
print("__________________________\n")

# Transferencia.
print("Transferencia")
if len(cadenas_trans)>0:
    peso_t_B=[]
    peso_t_A=[]
    peso_t_I=[]
    peso_cadenas_t=[]
    peso_total_t=float()
    peso_total_t_Xp=float()
    num_t_I=int()
    num_t_A=int()
    num_t_B=int()
    for i in range(len(cadenas_trans)):
        peso_t_i=cadenas_trans[i].count(I)*MM_I
        peso_t_a=cadenas_trans[i].count(A)*MM_A
        peso_t_b=cadenas_trans[i].count(B)*MM_B
        peso_cadenas_t.append(peso_t_a+peso_t_b+peso_t_i)
        peso_total_t=peso_total_t+peso_cadenas_t[i]
        peso_total_t_Xp=peso_total_t_Xp+peso_cadenas_t[i]**2
        peso_t_I.append(cadenas_trans[i].count(I)*MM_I)
        peso_t_A.append(cadenas_trans[i].count(A)*MM_A)
        peso_t_B.append(cadenas_trans[i].count(B)*MM_B)
        num_t_I=num_t_I+cadenas_trans[i].count(I)
        num_t_A=num_t_A+cadenas_trans[i].count(A)
        num_t_B=num_t_B+cadenas_trans[i].count(B)
    Xn_t=peso_total_t/len(peso_cadenas_t)
    Xp_t=peso_total_t_Xp/(peso_total_t)
    polimol_t=Xp_t/Xn_t
    Lambda_t=(num_t_A+num_t_B)/len(cadenas_trans)
    print("pesos A = "+str(peso_t_A))
    print("pesos B = "+str(peso_t_B))
    print("pesos I = "+str(peso_t_I))
    print("pesos = "+str(peso_cadenas_t))
    print("peso total = "+str(peso_total_t))
    print("peso total Xp= "+str(peso_total_t_Xp))
    print("Xn Transferencia = "+str(Xn_t))
    print("Xp Transferencia = "+str(Xp_t))
    print("Polimolecularidad Transferencia = "+str(polimol_t))
    print("Lambda de transferencia = "+str(Lambda_t))
    print("\n")
else:
    print("No se toma en cuenta porque no hay.")
    peso_t_I=0
    peso_t_A=0
    peso_t_B=0
    peso_total_t=0
    peso_total_t_Xp=0
    num_t_I=0
    num_t_A=0
    num_t_B=0
    Xn_t=0
    Xp_t=0
    polimol_t=0
    Lambda_t=0
print("__________________________\n")

# Recombinación
print("Recombinación")
if len(cadenas_recomb)>0:
    peso_r_B=[]
    peso_r_A=[]
    peso_r_I=[]
    peso_cadenas_r=[]
    peso_total_r=float()
    peso_total_r_Xp=float()
    num_r_I=int()
    num_r_A=int()
    num_r_B=int()
    for i in range(len(cadenas_recomb)):
        peso_r_i=cadenas_recomb[i].count(I)*MM_I
        peso_r_a=cadenas_recomb[i].count(A)*MM_A
        peso_r_b=cadenas_recomb[i].count(B)*MM_B
        peso_cadenas_r.append(peso_r_a+peso_r_b+peso_r_i)
        peso_total_r=peso_total_r+peso_cadenas_r[i]
        peso_total_r_Xp=peso_total_r_Xp+peso_cadenas_r[i]**2
        peso_r_I.append(cadenas_recomb[i].count(I)*MM_I)
        peso_r_A.append(cadenas_recomb[i].count(A)*MM_A)
        peso_r_B.append(cadenas_recomb[i].count(B)*MM_B)
        num_r_I=num_r_I+cadenas_recomb[i].count(I)
        num_r_A=num_r_A+cadenas_recomb[i].count(A)
        num_r_B=num_r_B+cadenas_recomb[i].count(B)
    Xn_r=peso_total_r/len(peso_cadenas_r)
    Xp_r=peso_total_r_Xp/(peso_total_r)
    polimol_r=Xp_r/Xn_r
    Lambda_r=(num_r_A+num_r_B)/len(cadenas_recomb)
    print("pesos A = "+str(peso_r_A))
    print("pesos B = "+str(peso_r_B))
    print("pesos I = "+str(peso_r_I))
    print("pesos = "+str(peso_cadenas_r))
    print("peso total = "+str(peso_total_r))
    print("peso total Xp= "+str(peso_total_r_Xp))
    print("Xn recombinación = "+str(Xn_r))
    print("Xp Recombinación = "+str(Xp_r))
    print("Polimolecularidad Recombinación = "+str(polimol_r))
    print("Lambda de recombinación = "+str(Lambda_r))
    print("\n")
else:
    print("No se toma en cuenta porque no hay.")
    peso_r_I=0
    peso_r_A=0
    peso_r_B=0
    peso_total_r=0
    peso_total_r_Xp=0
    num_r_I=0
    num_r_A=0
    num_r_B=0
    Xn_r=0
    Xp_r=0
    polimol_r=0
    Lambda_r=0
print("__________________________\n")

# Totales. Se toman en cuenta todas las terminaciones.
print("Estadísticas del polímero.")
Xn_Total=(peso_total_d+peso_total_t+peso_total_r)/(len(cadenas_dism)+len(cadenas_trans)+len(cadenas_recomb))
Xp_Total=(peso_total_d_Xp+peso_total_t_Xp+peso_total_r_Xp)/(peso_total_d+peso_total_t+peso_total_r)
Polimol_Total=Xp_Total/Xn_Total
Lambda_Total=(num_d_A+num_d_B+num_t_A+num_t_B+num_r_A+num_r_B)/(len(cadenas_dism)+len(cadenas_trans)+len(cadenas_recomb))
print("El número de iniciadores que hay en el polímero es: "+str(num_d_I+num_t_I+num_r_I))
print("El número de Monómeros A que hay en el polímero es: "+str(num_d_A+num_t_A+num_r_A))
print("El número de Monómeros B que hay en el polímero es: "+str(num_d_B+num_t_B+num_r_B))
print("Xn del polímero = "+str(Xn_Total))
print("Xp del polímero = "+str(Xp_Total))
print("Polimolecularidad del polímero = "+str(Polimol_Total))
print("Lambda del polímero = "+str(Lambda_Total))
print("Número total de cadenas = "+str(len(cadenas_dism)+len(cadenas_trans)+len(cadenas_recomb)))
print("_________________________________________________\n")

print("Datos curiosos.\n")
repa=rep_a*"A"
repb=rep_b*"B"
# Todas las combinaciones posibles para que nunca marque error.
if len(cadenas_dism)>0 and len(cadenas_trans)>0 and len(cadenas_recomb)>0: # Así nunca marca error.
    cadenas_Totales=cadenas_dism+cadenas_trans+cadenas_recomb
elif len(cadenas_trans)>0 and len(cadenas_recomb)>0 and len(cadenas_dism)==0:
    cadenas_Totales=cadenas_trans+cadenas_recomb
elif len(cadenas_dism)>0 and len(cadenas_trans)>0 and len(cadenas_recomb)==0:
    cadenas_Totales=cadenas_dism+cadenas_trans
elif len(cadenas_dism)>0 and len(cadenas_recomb)>0 and len(cadenas_trans)==0:
    cadenas_Totales=cadenas_dism+cadenas_recomb
elif len(cadenas_dism)>0 and len(cadenas_trans)==0 and len(cadenas_recomb)==0:
    cadenas_Totales=cadenas_dism
elif len(cadenas_dism)==0 and len(cadenas_trans)>0 and len(cadenas_recomb)==0:
    cadenas_Totales=cadenas_trans
elif len(cadenas_dism)==0 and len(cadenas_trans)==0 and len(cadenas_recomb)>0:
    cadenas_Totales=cadenas_recomb
    
# No pude hacer que contara las repeticiones máximas de cada monómero.
if len(cadenas_Totales)>0 and datoscuriosos=="y":
    rep__a=int()
    rep__b=int()
    rep__ab=int()
    rep__ba=int()
    rep__ia=int()
    rep__ai=int()
    rep__ib=int()
    rep__bi=int()
    rep__a_max=[]
    reo__b_max=[]
    rep__at=int()
    rep__bt=int()
    rep__abt=int()
    rep__bat=int()
    rep__iat=int()
    rep__ait=int()
    rep__ibt=int()
    rep__bit=int()
   #rep__max_at=[]
   #rep__max_bt=[]
   #As=[]
   #Bs=[]
    rep_max=[]
    lastchar = ''
    charcount = 0
    tmpcount = 1
    for i in range(len(cadenas_Totales)):
        rep__a=cadenas_Totales[i].count(repa)
        rep__b=cadenas_Totales[i].count(repb)
        rep__ab=cadenas_Totales[i].count("AB")
        rep__ba=cadenas_Totales[i].count("BA")
        rep__ia=cadenas_Totales[i].count("IA")
        rep__ai=cadenas_Totales[i].count("AI")
        rep__ib=cadenas_Totales[i].count("IB")
        rep__bi=cadenas_Totales[i].count("BI")
        #As=len(re.findall(r'A+A',cadenas_Totales[i]))   Codigo fallido :(
        #rep__max_a=len(max(As[i]))
        #Bs=len(re.findall(r'B+B',cadenas_Totales[i]))
        #rep__max_b=len(max(Bs[i]))
        #As.append(len(re.findall(r'A',cadenas_Totales[i])))
        #Bs.append(len(re.findall(r'B',cadenas_Totales[i])))
        rep__at=rep__at+rep__a
        rep__bt=rep__bt+rep__b
        rep__abt=rep__abt+rep__ab
        rep__bat=rep__bat+rep__ba
        rep__iat=rep__iat+rep__ia
        rep__ait=rep__ait+rep__ai
        rep__ibt=rep__ibt+rep__ib
        rep__bit=rep__bit+rep__bi
        #rep__max_at.append(As) # Arreglo con la incidencia máxima de A repetida
        #rep__max_bt.append(Bs)
    #rep___max_at=max(As)
    #rep___max_bt=max(Bs)
    print("Número veces que hay "+str(rep_a)+" A consecutivos en el polímero = "+str(rep__at))
    print("Número veces que hay "+str(rep_b)+" B consecutivos en el polímero "+str(rep__bt))
    print("Número de repeticiones de AB = "+str(rep__abt))
    print("Número de repeticiones de BA = "+str(rep__bat))
    print("Número de repeticiones de IA = "+str(rep__iat))
    print("Número de repeticiones de AI = "+str(rep__ait))
    print("Número de repeticiones de IB = "+str(rep__ibt))
    print("Número de repeticiones de BI = "+str(rep__bit))
    #print("Repeticiones máximas de monómero por cadena = "+str(rep_max))
    #print("Número máximo de veces que B se repite en cada cadena = "+str(Bs))
    #print("Número máximo de repeticiones de A en el polímero = "+str(rep___max_at))
    #print("Número máximo de repeticiones de B en el polímero = "+str(rep___max_bt))
    #print("Número de veces que se repite el número máximo de repeticiones de A en el polímero = "+str(rep__max_at.count(rep___max_at)))
    #print("Número de veces que se repite el número máximo de repeticiones de B en el polímero = "+str(rep__max_bt.count(rep___max_bt)))

#print(max(len(list(y)) for (c,y) in itertools.groupby(cadenas_Totales) if c=='B'))
#print(max(len(list(y)) for (c,y) in itertools.groupby(cadenas_Totales) if c=='A'))

print("Cadenas\n")
print("Dismutación\n")
print(cadenas_dism)
print("Número de cadenas por dismutación = "+str(len(cadenas_dism))+"\n__________________________________")
print("Transferencia\n")
print(cadenas_trans)
print("Número de cadenas por transferencia = "+str(len(cadenas_trans))+"\n__________________________________")
print("Recombinación\n")
print(cadenas_recomb)
print("Número de cadenas por recombinación = "+str(len(cadenas_recomb))+"\n__________________________________")
print("Probabilidades Finales")
print("Probabilidad de I con A = "+str(p_ia))
print("Probabilidad de I con B = "+str(p_ib))
print("Probabilidad de A con A = "+str(p_aa))
print("Probabilidad de A con B = "+str(p_ab))
print("Probabilidad de B con B = "+str(p_bb))
print("Probabilidad de B con A = "+str(p_ba))
print("\n")
print("Constantes Utilizadas")
print("k_IA = "+str(kia))
print("k_IB = "+str(kib))
print("k_AA = "+str(kaa))
print("k_AB = "+str(kab))
print("k_BB = "+str(kbb))
print("k_BA = "+str(kba))
print("r_A = "+str(ra))
print("r_B = "+str(rb))
