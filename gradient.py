# Chargement des bibliothèques
import matplotlib.pyplot as plt
import numpy as np

def gradml(x):
    mu0 = g_reforme_gene[0, 0, 0]
    dmldalpha = 0
    dmldbeta  = 0
    dmldgamma = 0
    dmlddelta = 0
    dmldsigma = 0

#    alphaA, betaA, gammaA, deltaA, alphaM, betaM, gammaM, deltaM, sigma = x
    alphaA, betaA, gammaA, deltaA, sigma = x
#    alphaM, betaM, gammaM, deltaM, sigma = x

#    for c in range(C): # entrée dans la boucle sur les conditions
    for c in range(C-1): # entrée dans la boucle sur la condition alginate
#    for c in range(1, C): # entrée dans la boucle sur la condition maltose
        eta = eta[c, :]
        if c == 0: # faisons d'abord la condition alginate
            alpha, beta, gamma, delta = alphaA, betaA, gammaA, deltaA
        if c == 1:
            alpha, beta, gamma, delta = alphaM, betaM, gammaM, deltaM
        for n in range(N): # entrée dans la boucle sur les temps
            t = T[n]
            ### calcul d'une somme utile dans mu_{c}(t) qui dépend de t, donc de n
            somme = 0
            smudelta = 0
            coef = np.zeros(n)
            for j in range(n): # à faire : améliorer rectangle -> trapèze
                coef[j] = ( eta[j+1] + eta[j] ) / 2
#                somme += (np.exp(delta*T[j+1])-np.exp(delta*T[j])) * eta[j]
                somme += (np.exp(delta*T[j+1])-np.exp(delta*T[j])) * coef[j]
                # calcul de la somme pour la dérivée de mu par rapport à delta
                smudelta += (T[j+1]*np.exp(delta*T[j+1])-T[j]*np.exp(delta*T[j])) * coef[j]
            ### fin de ce calcul
            # calcul de mu_c(t)
            muct = (mu0 - alpha/delta) * np.exp(-delta*t) + alpha/delta + beta*np.exp(-delta*t)/delta*somme
            # calcul des dérivées de mu par rapport aux paramètres alpha, beta, delta et mu0
            dmudalpha = (1 - np.exp(-delta*t)) / delta
            dmudbeta  = np.exp(-delta*t) / delta * somme
            dmuddelta = alpha*np.exp(-delta*t)/delta**2 - (mu0-alpha/delta)*t*np.exp(-delta*t) - alpha/delta**2 - beta*np.exp(-delta*t)*(1+1/delta)*somme/delta + beta*np.exp(-delta*t)*smudelta/delta
            dmudmu0 = np.exp(-delta*t)
            for r in range(R):
                # on prend g_{c, r}(t) dans le vecteur g_reforme_gene
                gcrt = g_reforme_gene[c, r, n]
                # calcul des dérivées par rapport à sigma, alpha, beta, delta et mu0
                dmldsigma += (np.log(gcrt/muct) / sigma + sigma/2) * (-np.log(gcrt/muct) / (sigma**2) + 1/2) + 1/sigma
                dmldalpha += -dmudalpha * (np.log(gcrt/muct) / sigma + sigma/2) / muct
                dmldbeta  += -dmudbeta  * (np.log(gcrt/muct) / sigma + sigma/2) / muct
                dmlddelta += -dmuddelta * (np.log(gcrt/muct) / sigma + sigma/2) / muct
                dmldmu0   += -dmudmu0   * (np.log(gcrt/muct) / sigma + sigma/2) / muct

    return np.array([dalphaA, dbetaA, ddeltaA, dmu0A, dalphaM, dbetaM, ddeltaM, dmu0M, dsigma])
