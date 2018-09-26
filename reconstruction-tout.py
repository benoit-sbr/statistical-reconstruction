# Chargement des bibliothèques
from io import BytesIO
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

# Chargement des fonctions
from fonctions import * # Packages and modules should have short, all-lowercase names 
# http://www.python.org/dev/peps/pep-0008/

# 1. vecteur des expressions g_{c, r}(t)
# mes nombres issus de excel ont ',' comme séparateur avant la partie décimale
gbinaire        = open('donnees/donnees_tous_genes_gmkA.txt').read().replace(',','.').encode()
g_tout          = np.genfromtxt(BytesIO(gbinaire), missing_values = 'NA', filling_values = np.nan, skip_header = 1)
#g_reforme_tout  = g_tout[ : , 3 : ].reshape(2, 4, 11, 17)	# vraiment tout
#g_reforme_tout  = g_tout[ : , 4 : 16 ].reshape(2, 4, 11, 12)	# grand opéron
g_reforme_tout  = g_tout[ : , 4 : 8 ].reshape(2, 4, 11, 4)	# gènes 2613 à 2616
#g_reforme_tout  = g_tout[ : , 8 ].reshape(2, 4, 11, 1)		# gène 2617
# 2. vecteur des temps
T		= np.array([0., 0.5, 1., 3.5, 6.5, 9.25, 12.5, 18.5, 22.5, 26.5, 30.5])
# 3. petit amusement
C, R, N, G      = g_reforme_tout.shape
# 4. condition
condition = 0 # alginate
# condition = 1 # maltose
# 5. expressions initiales
# on prend g_{condition, r}(0) dans le vecteur g_reforme_tout et on fait la moyenne sur les réplicats
mu0             = np.nanmean(g_reforme_tout[condition, :, 0, :], axis = 0)

###################################
# initial guess pour l'optimisation
x0 = np.ones(G*5 + T.size - 2)
bnds = tuple((0., None) for x in x0)

# résultats de l'optimisation
res1 = optimize.minimize(ml, x0, bounds = bnds, args = (C, R, N, G, condition, g_reforme_tout, mu0, T))
eta		= np.zeros(T.size - 1)
eta[0]		= 1.
eta[1 : ]	= res1.x[-10 : -1]
theta		= np.zeros((4, G))
Sigma		= np.zeros(G)
for gene_index in range(G): # gene_index varie de 0 à G-1
  theta[ : , gene_index]	= res1.x[5*gene_index : 5*gene_index + 4]
  Sigma[gene_index]		= res1.x[5*gene_index + 4]
  alpha = theta[0, gene_index]
  beta  = theta[1, gene_index]
  gamma = theta[2, gene_index]
  delta = theta[3, gene_index]
  sigma = Sigma[gene_index]

# théorie
tempsvisu = np.arange(T[0], T[-1], 0.01)				# TODO : faire un choix pertinent du pas de temps
resultat  = theorie(tempsvisu, T, eta, mu0[3], theta[ : , 3])

# mesures
observation = mesure(g_reforme_tout, 3)

# 1. Create a figure of size 8x6 inches, 80 dots per inch
figure1 = plt.figure(figsize = (8, 6), dpi = 80)
figure1.canvas.set_window_title('My title')

#plt.subplot(1,2,1)
# 2. Plot (T, resultat)
plt.plot(tempsvisu, resultat, '-', label = 'mu alginate')
plt.plot(T, observation, 'o', label = 'g alginate')
plt.legend()
plt.grid(True)
"""
plt.subplot(1,2,2)
# 2. Plot (T, resultatM)
plt.plot(tempsvisu, resultatM, '-', label = 'mu maltose')
plt.plot(T, observationM, 'o', label = 'g maltose')
plt.legend()
plt.grid(True)
"""
# 3. Show result on screen
plt.show()
"""
# prepare data to write in a file
data = np.column_stack([tempsvisu, resultatA])

with open('fichier-test.txt', 'wb') as f: # here you open the ascii file
  np.savetxt(f, data, fmt = '%.2f')  # here the ascii file is populated
"""
