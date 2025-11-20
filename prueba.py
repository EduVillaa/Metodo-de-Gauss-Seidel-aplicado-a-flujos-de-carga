import numpy as np

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

S_pq=np.array([-4-2.5j])
P_pv=np.array([2])



c=0
sumQ=0
while c<Y.shape[0]:
    sumQ = sumQ + Y[Y.shape[0]-1][c]*V[c] 
    c+=1
S_pv = np.array([P_pv[0]-np.imag(sumQ*np.conjugate(V_pv[0]))*1j])

S_i = np.concatenate((S_pq, S_pv))    

for i in range(2, Y.shape[0]+1):
    sum1, sum2 =0,0
    for j in range(1, i):
        sum1=sum1+Y[i-1][j-1]*V[j-1]

    for j in range(i+1, Y.shape[0]+1):
        sum2=sum2+Y[i-1][j-1]*V[j-1]

    V_i=(np.conjugate(S_i[i-2])/V[i-1]-sum1-sum2)/Y[i-1][i-1]

    V[i-1]=V_i

print(np.round(V, 4))

