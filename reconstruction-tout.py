# Chargement des bibliothèques
from io import BytesIO
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit, minimize

# Chargement des fonctions
from fonctions import * # Packages and modules should have short, all-lowercase names 
# http://www.python.org/dev/peps/pep-0008/

# 1. vecteur des expressions g_{c, r}(t)
# mes nombres issus de excel ont ',' comme séparateur avant la partie décimale
gbinaire        = open('donnees/donnees_tous_genes_gmkA.txt').read().replace(',','.').encode()
g_tout          = np.genfromtxt(BytesIO(gbinaire), missing_values = 'NA', filling_values = np.nan, skip_header = 1)
#g_temp_tout  = g_tout[ : , 3: ].reshape(2, 4, 11, 17)	# vraiment tout
g_temp_tout  = g_tout[ : , 4:16 ].reshape(2, 4, 11, 12)	# grand opéron
#g_temp_tout  = g_tout[ : , 4:8 ].reshape(2, 4, 11, 4)	# gènes 2613 à 2616
g_temp_TF    = g_tout[ : , 8 ].reshape(2, 4, 11)	# gène 2617
# 1.1 amélioration ? Les derniers temps ne sont pas bons...
nb_pas_temps_enleves = 0
#g_reforme_tout = g_temp_tout[ : , : , :-nb_pas_temps_enleves , : ]
#g_reforme_TF   = g_temp_TF[ : , : , :-nb_pas_temps_enleves ]
g_reforme_tout = g_temp_tout[ : , : , : , : ]
g_reforme_TF   = g_temp_TF[ : , : , : ]

# 2. vecteur des temps
T		= np.array([0., 0.5, 1., 3.5, 6.5, 9.25, 12.5, 18.5, 22.5, 26.5, 30.5])
# 3. petit amusement
C, R, N, G      = g_reforme_tout.shape
# 4. condition
condition = 0 # alginate
# condition = 1 # maltose
# 5. expressions initiales
# on prend g_{condition, r}(0)
# et on fait la moyenne sur les réplicats
mu0             = np.nanmean(g_reforme_tout[condition, :, 0, :], axis = 0)

# initial guess pour l'optimisation
x0 = np.zeros(G*5 + T.size - 2 - nb_pas_temps_enleves)
# initial guess pour les alphas
x0[0 : G*5 : 5] = 0.
# initial guess pour les betas
# on prend g_{condition, r}(1)
# et on fait la moyenne sur les réplicats
x0[1 : G*5 + 1 : 5] = np.nanmean(g_reforme_tout[condition, :, 1, :], axis = 0)

# initial guess pour les deltas
# on prend g_{condition, r}(3:7)
# on approche cela par une exponentielle décroissante
# et on fait la moyenne sur les réplicats
def expo(t, delta):
    global g3
    return g3 * np.exp(-delta * t)
for gene_index in range(G): # gene_index varie de 0 à G-1
    xdata = T[3:8] - T[3]
    ydata = g_reforme_tout[condition, : , 3:8, gene_index]
    for replicat in (0, 1, 3): # et non pas in range(R) car le réplicat #2 a un NaN que je ne sais pas encore comment traiter
        g3 = g_reforme_tout[condition, replicat, 3, gene_index]
        popt, pcov = curve_fit(expo, xdata, ydata[replicat])
    x0[3 + gene_index*5] = np.nanmean(popt)
# initial guess pour les sigmas
x0[4 : G*5 + 4 : 5] = 1.
# initial guess pour les etas
g_TF_renormalise = np.zeros((R, N-2))
for replicat in range(R):
    g_TF_renormalise[replicat, :] = (g_reforme_TF[condition, replicat, 1:-1] + g_reforme_TF[condition, replicat, 2:]) / (2 * g_reforme_TF[condition, replicat, 0])
x0[-T.size+2+nb_pas_temps_enleves :] = np.nanmean(g_TF_renormalise, axis = 0)
# initial guess pour les gammas
x0[2 : G*5 + 2 : 5] = x0[-T.size+2 :].max() / 2

bnds = tuple((0., None) for x in x0)

# résultats de l'optimisation
res1 = minimize(ml, x0, bounds = bnds, args = (C, R, N, G, condition, g_reforme_tout, mu0, T))
eta		= np.zeros(T.size - 1)
eta[0]		= 1.
eta[1 : ]	= res1.x[-eta.size : -1]
theta		= np.zeros((R, G))
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
#plt.plot(T[0:-nb_pas_temps_enleves], observation, 'o', label = 'g alginate')
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
