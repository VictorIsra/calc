import sys
import numpy as np
import scipy as sp
from math import floor
from scipy.sparse import linalg
import matplotlib.pyplot as plt
from numpy import linalg as la

comprimento_grid = float(input("entre com o comprimento do grid (ex: 32): "))# apenas o comprimento do grid, nao o total de pontos contidos nele
temperatura_l = [ float(t) for t in input("digite as 6 temperaturas pros lados da malha L (ex: 100 90 72.5 65 -14 0): ").split()  ]

if(comprimento_grid <3):
	print("***COMPRIMENTO DO GRID DEVE SER NO MÍNIMO 3!***\n")
	sys.exit()
elif len(temperatura_l) != 6:
		print("***É NECESSÁRIO PASSAR AS 6 TEMPERATURAS!***\n")
		sys.exit()

for index, t in enumerate(temperatura_l):
	if t == 1.:
		temperatura_l[index] = t+1

grid = np.ones((comprimento_grid, comprimento_grid))#criando grid
cordenadas_contorno = []	#quardará todas as codenadas do grid que estiverem ma borda da malha em L ou dos pontos fora dela

#vertices da malha em L dentro do grid ( cordenadas no estilo matriz, nao cartesiano)
v = [ (0,0), (0,floor(2/3*(comprimento_grid-1))), ((floor(comprimento_grid-1)/2), floor(2/3*(comprimento_grid-1))), ((floor(comprimento_grid-1)/2), floor((comprimento_grid-1)/2)),(comprimento_grid-1,comprimento_grid-1),(comprimento_grid-1,0) ]

#preenche o grid com as temperaturas de contorno e com 1's dentro de L 
grid[ 0, :] = temperatura_l[0]
grid[ 1:v[2][0], v[1][1] ] = temperatura_l[1]
grid[ v[2][0]-1, v[1][1]-1: ] = temperatura_l[2]
grid[ v[2][0]-1:,-1] = temperatura_l[3]
grid[ -1, : ] = temperatura_l[4]
grid[1: , 0 ] = temperatura_l[5]
grid[ :v[2][0]-1, v[1][1]: ] = temperatura_l[0]

qtdade_pontos_l = 0 #guardará a quantidade de pontos no interior da malha em L
for index, ponto in np.ndenumerate(grid): #conta quantos pontos existem dentro da malha em L (desconsidera bordas)
	if(ponto == 1):
		qtdade_pontos_l = qtdade_pontos_l + 1
	else:
		cordenadas_contorno = cordenadas_contorno + [ index ]
			
A = np.zeros((qtdade_pontos_l,qtdade_pontos_l))#matriz A terá as dimensoes da quantidade de pontos contidos na malha L
B = np.arange(qtdade_pontos_l, dtype=float)   #vetor B que conterá as temperaturas de contorno dos pontos da malha L
temp_fronteira = 0 #recebera a temperatura resultante nas fronteiras de um ponto de L

u = 0
for index, x in np.ndenumerate(grid):
	if(index not in cordenadas_contorno): #checa se ponto pertence a malha L
		grid[index] = u #atribui um rotúlo(u) para os pontos que pertençam a malha L
		u = u + 1
linha_A = 0
for index, u in np.ndenumerate(grid):	 #calcula as diferencas finitas pros pontos no interior da malha em L	
		if(index not in cordenadas_contorno ):
			#para cada ponto de L é guardado os valores ( indices pares ) e coordenadas ( indices impares ) dos
			#quatro pontos ao seu redor
			u =  grid[(index[0], index[1])], (index[0], index[1]), grid[( index[0] , index[1]+1 )],( index[0] , index[1]+1 ), grid[( index[0]-1 , index[1] )],  (index[0]-1 , index[1] ), grid[( index[0] , index[1]-1 )],( index[0] , index[1]-1 ), grid[( index[0]+1 , index[1] )], (index[0]+1, index[1]) 
			for indice,l_u in enumerate(u):
				if(indice % 2 != 0):
					if(l_u in cordenadas_contorno):#checa se ponto é de fronteira e caso seja guarda sua temperatura
						temp_fronteira = temp_fronteira + u[indice-1]	
						
				else:
					if(u[indice+1] not in cordenadas_contorno):
						A[linha_A, l_u ] = 1	
						
			B[linha_A] = temp_fronteira
			temp_fronteira = 0				
			linha_A = linha_A + 1

diagonal_principal = -5*np.eye(qtdade_pontos_l)	
A = A + diagonal_principal
B = -1*B# para as temperaturas ficarem negativas ao passar pro outro lado da equacao
#print("dimensoes de A:",np.shape(A),"MATRIZ A:", A, "\ndimensoes de B:",np.shape(B),"\nVETOR B:", B,"\n")
print("CALCULANDO...")
x =sp.sparse.linalg.spsolve(A,B)#resolve o sistema Ax=b
#print("B", B)
#print("dimensoes de x:",np.shape(x),"VETOR x:", x,"\n")
i=0
for index, pontos in np.ndenumerate(grid):#completa os pontos dentra da malha com as temperaturas calculadas pelo sistema Ax=b
	if(index not in cordenadas_contorno):
		grid[index] = x[i]
		i = i + 1		
#print("dimensoes do grid:",np.shape(grid) ,"temperaturas no interior de L e do grid:", grid)
print("\nTOTAL DE PONTOS NO INTERIOR DA MALHA EM L:" , qtdade_pontos_l )

#plota o grid com a malha em L preenchida com suas temperaturas internas'''
axis0 = axis1 = np.linspace(0, comprimento_grid, comprimento_grid)
X, Y = np.meshgrid(axis0,axis1)

plt.pcolor(X, Y, grid)
plt.axis("scaled")
plt.gca().invert_yaxis()
plt.colorbar()
plt.xlabel(r"$X$ (m)")
plt.ylabel(r"$Y$ (m)")
plt.title("$TEMPERATURA$ $T(x,y)$ $NO$ $INTERIOR$ $DA$ $MALHA$ $L$ $EM$ $UM$ $GRID$ $DE$ %dx%d: $(m)$\n$TOTAL$ $DE$ $PONTOS$ $N$O $INTERIOR$ $DA$ $MALHA:$ %d $PONTOS$\n" % 	(comprimento_grid,comprimento_grid,qtdade_pontos_l))
plt.show()
