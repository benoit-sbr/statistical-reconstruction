# Chargement des bibliothèques
from io import BytesIO
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

# Chargement des fonctions
from gradient import * # Packages and modules should have short, all-lowercase names 
# http://www.python.org/dev/peps/pep-0008/

# 1. vecteur des expressions g_{c, r}(t)
gene		= '2616'
# mes nombres issus de excel ont ',' comme séparateur avant la partie décimale
gbinaire        = open('donnees/gene_'+gene+'_gmkA.txt').read().replace(',','.').encode()
g_gene		= np.genfromtxt(BytesIO(gbinaire), missing_values = 'NA', filling_values = np.nan, skip_header = 1)
g_reforme_gene	= g_gene.reshape(2, 4, 11)
# 2. vecteur des temps
T = np.array([0, 0.5, 1, 3.5, 6.5, 9.25, 12.5, 18.5, 22.5, 26.5, 30.5])
# 3. petit amusement
C, R, N, G = g_reforme_gene.shape + (1,) # append to a tuple
# 4. condition
condition = 0 # alginate
# condition = 1 # maltose
# 5. expression initiale
# on prend g_{condition, r}(0) dans le vecteur g_reforme_gene et on fait la moyenne sur les réplicats
mu0             = np.nanmean(g_reforme_gene[condition, :, 0])
# 6. vecteur des eta valeurs des TFX utilisées comme proxy pour les TFA
gbinaire        = open('donnees/gene_2617_gmkA.txt').read().replace(',','.').encode()
g_TF		= np.genfromtxt(BytesIO(gbinaire), missing_values = 'NA', filling_values = np.nan, skip_header = 1)
g_reforme_TF	= np.nanmean(g_TF.reshape(2, 4, 11), axis = 1)
eta		= g_reforme_TF[condition, :]

###################################
# fonctions qui pourraient être mises dans un autre fichier...
def theorie(tempsvisu):
    Nvisu = tempsvisu.size
    coef = np.interp(tempsvisu, T, eta)
    resultat = np.zeros(Nvisu)
    for nvisu in range(Nvisu): # entrée dans la boucle sur les temps
        tvisu = tempsvisu[nvisu]
        # calcul de l'integrale
        tauvisu = tempsvisu[:nvisu+1]
        etatau  = coef[:nvisu+1]
        integrale = np.trapz(np.exp(delta*tauvisu) * etatau, tauvisu)
        # calcul de mu_c(t)
        resultat[nvisu] = (mu0 - alpha/delta) * np.exp(-delta*tvisu) + alpha/delta + beta*np.exp(-delta*tvisu)/delta*integrale
    return resultat

def mesure(temps):
    # on prend g_{0, r}(t) dans le vecteur g_reforme_gene et on fait la moyenne sur les réplicats
    resultatA = np.nanmean(g_reforme_gene[0, :, :], axis = 0)
    # on prend g_{1, r}(t) dans le vecteur g_reforme_gene et on fait la moyenne sur les réplicats
    resultatM = np.nanmean(g_reforme_gene[1, :, :], axis = 0)
    return resultatA, resultatM

def ml(x): # moins l
    ml = 0
    alpha, beta, gamma, delta, sigma = x
    for n in range(N): # n varie de 0 à N-1
        t = T[n]
# calcul d'une somme utile dans mu_{c}(t) qui dépend de t, donc de n
        somme = 0
        for j in range(n): # j varie de 0 à n-1
            somme += (np.exp(delta*T[j+1]) - np.exp(delta*T[j])) * eta[j] / (gamma + eta[j])
# fin de ce calcul
# calcul de mu_c(t)
        muct = (mu0 - alpha/delta) * np.exp(-delta*t) + alpha/delta + beta*np.exp(-delta*t)/delta*somme
        for r in range(R):
# on prend g_{c, r}(t) dans le vecteur g_reforme_gene
            gcrt = g_reforme_gene[condition, r, n]
            if (np.isnan(gcrt)):
                continue
# à ce moment, on dispose de gcrt (et gcrt n'est pas un nan), et aussi de muct, on peut calculer la somme
            ml += (np.log(gcrt / muct) / sigma + sigma / 2)**2 / 2 + np.log(np.sqrt(2 * np.pi) * gcrt) + np.log(sigma)
    return ml

###################################
# initial guess pour l'optimisation
x0 = 2.*np.ones(G*5)
bnds = tuple((0., None) for x in x0)

# résultats de l'optimisation
res1 = optimize.minimize(ml, x0, bounds = bnds)
alpha, beta, gamma, delta, sigma = res1.x

## théorie
tempsvisu = np.arange(0, 30.5, 0.01) # à choisir -> faire un choix pertinent
resultatA = theorie(tempsvisu)

## mesures
observationA, observationM = mesure(T)

# 1. Create a figure of size 8x6 inches, 80 dots per inch
figure1 = plt.figure(figsize = (8, 6), dpi = 80)
figure1.canvas.set_window_title('My title')

plt.subplot(1,2,1)
# 2. Plot (T, resultatA)
plt.plot(tempsvisu, resultatA, '-', label = 'mu alginate')
plt.plot(T, observationA, 'o', label = 'g_gene_'+gene+' alginate')
plt.legend()
plt.grid(True)
"""
plt.subplot(1,2,2)
# 2. Plot (T, resultatM)
plt.plot(tempsvisu, resultatM, '-', label = 'mu maltose')
plt.plot(T, observationM, 'o', label = 'g_gene_'+gene+' maltose')
plt.legend()
plt.grid(True)
"""
# 3. Show result on screen
plt.show()
"""
# prepare data to write in a file
data = np.column_stack([tempsvisu, resultatA])

with open('fichier-test.txt', 'wb') as f: # here you open the ascii file
    np.savetxt(f, data, fmt = '%.2f')    # here the ascii file is populated
"""
