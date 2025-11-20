import numpy as np

numero_iteraciones = 2
Y = np.array(
    [
        [20-50j, -10+20j, -10+30j],
        [-10+20j, 26-52j, -16+32j],
        [-10+30j, -16+32j, 26-62j],
    ])

V_slack=np.array([1.05+0j])
V_pq=np.array([1+0j])
V_pv=np.array([1.04+0j])
V = np.concatenate((V_slack, V_pq, V_pv)) #Programa válido para redes con un solo nodo PV y lo colocamos el final
V_PU_sin_corregir=[]

S_pq=np.array([-4-2.5j])
P_pv=np.array([2])

for n in range(0, numero_iteraciones+1):
    c=0
    sumQ=0

    while c<Y.shape[0]:
        sumQ = sumQ + Y[Y.shape[0]-1][c]*V[c]
        c+=1
    S_pv = np.array([P_pv[0]-np.imag(sumQ*np.conjugate(V[Y.shape[0]-1]))*1j])
    S_i = np.concatenate((S_pq, S_pv))
    I_slack = 0 + 0j
    for a in range(Y.shape[0]):
        I_slack += Y[0, a] * V[a]
    S_slack = V_slack * np.conjugate(I_slack)
    print(f"Iteración {n}:")
    for t in range(1, Y.shape[0]+1):
         if t==1:
            print(f"Nodo slack: V_{t} = {V[t-1]:.4f}", "------------------ ",f"S_{t}={S_slack[0]:.4f}")
         elif 1<t<Y.shape[0]:
            print(f"Nodo PQ: V_{t} = {V[t-1]:.4f}", "---------------------- ",f"S_{t}={S_i[t-2]:.4f}")
         elif t==Y.shape[0]:
            print(f"Nodo PV: V_{t} correcto = {V[t-1]:.4f}", "-----------",f"S_{t}={S_i[len(S_i)-1]:.4f}")
        
    if n>0:
        sum3=0    
        for z in range(0, Y.shape[0]):
            sum3 = sum3 + abs(V[z]-Vanterior[z])
            #print(V[z], "-", Vanterior[z], V[z]-Vanterior[z],abs(V[z]-Vanterior[z]))
        print(f"Error = {sum3:.5f}")
    Vanterior = np.array(V)       
    for i in range(2, Y.shape[0]+1):
        sum1, sum2 =0,0
        for j in range(1, i):
            sum1=sum1+Y[i-1][j-1]*V[j-1]

        for j in range(i+1, Y.shape[0]+1):
            sum2=sum2+Y[i-1][j-1]*V[j-1]

        V_i=(np.conjugate(S_i[i-2])/np.conjugate(V[i-1])-sum1-sum2)/Y[i-1][i-1]
        V[i-1]=V_i
    
    V_PU_sin_corregir.append(V[Y.shape[0]-1])    
    if n>0:
        print(f"Nodo PV: V_{Y.shape[0]} sin corregir = {V_PU_sin_corregir[n-1]:.4f}")

    V[Y.shape[0]-1]=V_pv[0]*(V[Y.shape[0]-1]/np.abs(V[Y.shape[0]-1]))
    print("")
    V = np.round(V, 4)
    
