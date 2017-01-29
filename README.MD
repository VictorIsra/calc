**requisitos:**<br />python3+<br />Numpy<br />Scipy<br />matplotlib
#Calcula as temperaturas nos pontos do interior de uma malha em forma de L qualquer
	
	**A ideia básica foi gerar um grid quadrado de comprimento determinado pelo usuario e desenhar uma
		figura em forma de L nesse grid. O programa sempre fatia o grid em um L de acordo com as proporcoes contidas no vetor 'v'. (linha 66)
		 A figura em L é delimitada ao se marcar no grid as temperaturas de contorno determinadas pelo usuario e definindo que o L é tudo dentro
		desta reigao. Depois o grid é só percorrido nos pontos que pertencam a L. estes pontos sao marcados com um
		rótulo u, além de cada ponto em L guardar o valor dos quatro pontos ao seu redor.
		
			A matriz A, do sistema Ax=b terá ORDEM igual as da quantidade de pontos que pertencam a L, de modo que as temperatuas
		em cada ponto de L em um mesmo nivel serão calculadas linha por linha (equacoes),e cada rotulo de L contido nessas equacoes
		indicará o pertencimento ou nao equacao. Sempre que pertencer, será marcado um '1' na cordenada correspondente em A, de modo que asssim
		o programa irá gerar a matriz A completando-a linha a linha. 
		
		O vetor B associa uma temperatura de contorno (ou sua ausencia se for homogêneo) a cada um dos pontos que pertencam a L
		
			tendo em maos 'A' e 'B', o programa calcula a temperatura de cada ponto em L, de modo que resulta em um vetor 'x' com as temperatuas
		contidas no interior de L na ordem dos rótulos, resolvendo o problema.
		
		sistema gerado/resolvido: Ax=B
				
			o programa plota o grid final com a malha em L preenchida com suas temperaturas internas calculadas, bem como
		o número de pontos dentro da malha. vale a pena observar que a entrada "comprimento do grid" nada tem a ver com a
		quantidade de pontos no interior de L, se nao o fato de serem diretamente proporcionais.
		
			o programa mostra a corretude das soluções mostrando que os vetores A,x e B gerados possuem as dimensoes corretas, além de 
		mostrar cada uma das matriz/vetores (A,x,B).É aconselhado entrar com um comprimento pequeno para poder
		observar as matrizes geradas como um todo
			por motivos de praticidade, as matriz/vetores 'A', 'x' e 'B' calculadas nao foram printadas, mas caso o leitor  queira checar
		a corretude da saída, basta descomentá-las e printa-las:
		
			LINHA PARA PRINTAR MATRIZ/VETORES 'A','x' e 'B:  linhas: 113 e 117 (descomente essas linhas para ver a matrize os e vetores gerados )
		**
