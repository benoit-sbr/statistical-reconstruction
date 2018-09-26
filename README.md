# statistical-reconstruction

reconstruction-gene.py s'attache à reconstruire gène par gène
- alpha, beta, gamma, delta, sigma soit 5 paramètres
*en supposant eta connu* (disons que eta qui est en fait simplement le TFX sert de proxy pour TFA (cf paragraphe 4.2 de Khanin et al))
avec comme données
- les 4*11 mesures

reconstruction-tout.py s'attache à reconstruire
- 17*5 paramètres + 10 valeurs pour eta
*en NE supposant donc PAS eta connu*
avec comme données
- les 17x4x11 mesures
