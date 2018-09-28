# Chargement des bibliothèques
import numpy as np

def theorie(tempsvisu, temps, eta, mu0, theta):
    # Docstring de la fonction
    """
    Renvoie l'expression de mu (équation juste après (9) dans le papier Khanin et al)
    Positional arguments:
    tempsvisu	--- vecteur de temps pour la visualisation
    temps		--- vecteur de subdivision du temps (instants de mesure)
    eta		--- vecteur des différents valeurs de eta dans la subdivision
    mu0		--- niveau d'expression à t = 0
    theta		--- theta = (alpha, beta, gamma, delta) vecteur des paramètres cinétiques du gène considéré
    Keyword arguments:
    """
    # Corps de la fonction
    alpha, beta, gamma, delta = theta
    Nvisu	= tempsvisu.size
    resultat	= np.zeros(Nvisu)
    for nvisu in range(Nvisu): # entrée dans la boucle sur les temps
        tvisu	= tempsvisu[nvisu]
# calcul de l'integrale
        tauvisu	  = tempsvisu[:nvisu+1]
        positions = np.searchsorted(temps, tauvisu, side='right') - 1
        etatau	  = eta[positions] / (gamma + eta[positions])
        integrale = np.trapz(np.exp(delta*tauvisu) * etatau, tauvisu)
# calcul de mu_c(t)
        resultat[nvisu] = (mu0 - alpha/delta) * np.exp(-delta*tvisu) + alpha/delta + beta*np.exp(-delta*tvisu)/delta*integrale
    return resultat

def mesure(g, gene_index):
    # Docstring de la fonction
    """
    Renvoie l'expression de g
    Positional arguments:
    gene_index	--- numéro du gène (ou indice de 0 à 16) # TODO : décider !
    Keyword arguments:
    """
    # Corps de la fonction
# on prend g_{0, r}(t) dans le vecteur g et on fait la moyenne sur les réplicats
    resultat = np.nanmean(g[0, :, :, gene_index], axis = 0)
    return resultat

def ml(x, C, R, N, G, c, g, mu0, temps): # moins l
    ml = 0
    eta		= np.zeros(temps.size - 1)
    eta[0]	= 1.
    eta[1 : ]	= x[-10 : -1]
    theta	= np.zeros((4, G))
    Sigma	= np.zeros(G)
    for gene_index in range(G): # gene_index varie de 0 à G-1
        theta[ : , gene_index]	= x[5*gene_index : 5*gene_index + 4]
        Sigma[gene_index]		= x[5*gene_index + 4]
        alpha = theta[0, gene_index]
        beta    = theta[1, gene_index]
        gamma = theta[2, gene_index]
        delta = theta[3, gene_index]
        sigma = Sigma[gene_index]
        for n in range(N): # n varie de 0 à N-1
            t = temps[n]
# calcul d'une somme utile dans mu_{c}(t) qui dépend de t, donc de n
            somme = 0
            for j in range(n): # j varie de 0 à n-1
                somme += (np.exp(delta*temps[j+1]) - np.exp(delta*temps[j])) * eta[j] / (gamma + eta[j])
# fin de ce calcul
# calcul de mu_c(t)
            muct = (mu0[gene_index] - alpha/delta) * np.exp(-delta*t) + alpha/delta + beta*np.exp(-delta*t)/delta*somme
            for r in range(R):
# on prend g_{c, r}(t) dans le vecteur g
                gcrt = g[c, r, n, gene_index]
                if (np.isnan(gcrt)):
                    continue
# à ce moment, on dispose de gcrt (et gcrt n'est pas un nan), et aussi de muct, on peut calculer la somme
                ml += (np.log(gcrt / muct) / sigma + sigma / 2)**2 / 2 + np.log(np.sqrt(2 * np.pi) * gcrt) + np.log(sigma)
    return ml
