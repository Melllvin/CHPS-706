from math import cos, pi, sin, sqrt, log


## Solution exacte du problème
## @param x, y : coordonnées du point
## @return : u(x, y)
def fct_u(x, y):
    return 1 + sin(x * pi / 2) + x * (x - 4) * cos(y*pi / 2 )

## Température extérieure sur le bord de référence 1
## @param x, y : coordonnées du point
## @return : uE(x, y)
def fct_uE(x, y):
    return 1

## Terme source f = -Δu
## @param x, y : coordonnées du point
## @return : f(x, y)
def fct_f(x, y):
    return ( pi**2 / 4 * sin(x*pi / 2) + (pi**2 / 4 * x**2 - pi**2 * x - 2) * cos(y * pi / 2))

## Conductivité thermique du matériau (constante ici)
## @param x, y : coordonnées du point
## @return : kappa(x, y)
def fct_kappa(x, y):
    return 1

## Coefficient de transfert (pénalisation de la condition de Dirichlet)
## @param x, y : coordonnées du point
## @return : alpha(x, y)
def fct_alpha(x, y):
    return 1e8


## Calcule la matrice de rigidité élémentaire du triangle Tl
## @param xl : abscisses des 3 sommets [x1, x2, x3]
## @param yl : ordonnées des 3 sommets [y1, y2, y3]
## @return : kl, matrice 3x3 sous forme de liste de listes
def coeffelem_P1_rigid(xl, yl):
    x1, x2, x3 = xl
    y1, y2, y3 = yl

    # Aire du triangle (déterminant, en valeur absolue)
    aire = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2

    kl = [[0.0] * 3 for _ in range(3)]

    # Conductivité au barycentre du triangle
    kappa = fct_kappa((x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3)

    # Termes diagonaux : longueur au carré du coté opposé au sommet
    kl[0][0] = kappa / (4 * aire) * ((x2 - x3) ** 2 + (y2 - y3) ** 2)
    kl[1][1] = kappa / (4 * aire) * ((x1 - x3) ** 2 + (y1 - y3) ** 2)
    kl[2][2] = kappa / (4 * aire) * ((x1 - x2) ** 2 + (y1 - y2) ** 2)

    # Termes hors diagonale : produit scalaire depuis le 3e sommet
    kl[0][1] = kappa / (4 * aire) * (-(x1 - x3) * (x2 - x3) - (y1 - y3) * (y2 - y3)) # 3e sommet = S3
    kl[0][2] = kappa / (4 * aire) * (-(x1 - x2) * (x3 - x2) - (y1 - y2) * (y3 - y2)) # 3e sommet = S2
    kl[1][2] = kappa / (4 * aire) * (-(x2 - x1) * (x3 - x1) - (y2 - y1) * (y3 - y1)) # 3e sommet = S1

    # Symétrie
    kl[1][0] = kl[0][1]
    kl[2][0] = kl[0][2]
    kl[2][1] = kl[1][2]

    return kl


## Calcule le vecteur source élémentaire du triangle Tl (quadrature au barycentre)
## @param xl : abscisses des 3 sommets [x1, x2, x3]
## @param yl : ordonnées des 3 sommets [y1, y2, y3]
## @return : fl, liste de 3 valeurs
def coeffelem_P1_source(xl, yl):
    x1, x2, x3 = xl
    y1, y2, y3 = yl

    # Aire du triangle
    aire = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2

    # Barycentre du triangle
    xG = (x1 + x2 + x3) / 3
    yG = (y1 + y2 + y3) / 3

    val_f = fct_f(xG, yG)

    fl = [aire / 3 * val_f] * 3

    return fl


## Calcule la matrice de poids élémentaire de l'arête de bord
## pa = longueur(Aa) / 6 * alpha * [[2, 1], [1, 2]]
## @param xa : abscisses des 2 sommets de l'arête [x1, x2]
## @param ya : ordonnées des 2 sommets de l'arête [y1, y2]
## @return : pa, matrice 2x2 sous forme de liste de listes
def coeffelem_P1_poids(xa, ya):
    x1, x2 = xa
    y1, y2 = ya

    # Longueur de l'arête (Pythagore)
    longueur = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5

    # Alpha au milieu de l'arête
    alpha = fct_alpha((x1 + x2) / 2, (y1 + y2) / 2)

    pa = [[0, 0], [0,0]]

    pa[0][0] = longueur / 6 * alpha * 2
    pa[1][1] = longueur / 6 * alpha * 2
    pa[0][1] = longueur / 6 * alpha * 1
    pa[1][0] = longueur / 6 * alpha * 1

    return pa


## Calcule le vecteur de transfert élémentaire de l'arête de bord
## ea_i = longueur(Aa) / 2 * alpha * uE(milieu)
## @param xa : abscisses des 2 sommets de l'arête [x1, x2]
## @param ya : ordonnées des 2 sommets de l'arête [y1, y2]
## @return : ea, liste de 2 valeurs
def coeffelem_P1_transf(xa, ya):
    x1, x2 = xa
    y1, y2 = ya

    # Longueur de l'arête
    longueur = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5

    # Milieu de l'arête
    xM = (x1 + x2) / 2
    yM = (y1 + y2) / 2

    alpha = fct_alpha(xM, yM)
    u_E = fct_uE(xM, yM)

    ea = [longueur / 2 * alpha * u_E * 1] * 2

    return ea



## Repris du TP1

## Lit le fichier .msh et retourne les données sous forme de liste
## @param chemin_relatif : chemin du fichier .msh
## @return : liste contenant les données du fichier .msh sous la forme [nbn,nbe,nba,coord,tri,ar,refn,reft,refa]
def lit_fichier_msh(chemin_relatif):
    with open(chemin_relatif, 'r') as fichier:
        lignes = fichier.readlines()

    nbn, nbe, nba = [int(x) for x in lignes[0].split()]

    # Coords + refn : lignes 1 à nbn
    coords = []
    refn = []
    debut = 1
    for i in range(nbn):
        valeurs = lignes[debut + i].split()
        coords.append([float(valeurs[0]), float(valeurs[1])])
        refn.append(int(valeurs[2]))

    # Tri + reft : lignes juste après les noeuds
    tri = []
    reft = []
    debut = 1 + nbn
    for i in range(nbe):
        valeurs = lignes[debut + i].split()
        tri.append([int(valeurs[0]), int(valeurs[1]), int(valeurs[2])])
        reft.append(int(valeurs[3]))

    # Ar + refa : lignes juste après les triangles
    ar = []
    refa = []
    debut = 1 + nbn + nbe
    for i in range(nba):
        valeurs = lignes[debut + i].split()
        ar.append([int(valeurs[0]), int(valeurs[1])])
        refa.append(int(valeurs[2]))

    resultat = [nbn, nbe, nba, coords, tri, ar, refn, reft, refa]

    return resultat


## Calcule le pas et la qualité du maillage
## @param nbe : nombre de triangles
## @param coord : coordonnées des noeuds [[x, y], ...]
## @param tri : sommets des triangles [[n1, n2, n3], ...] (numérotés à partir de 1)
## @return : [h, Q] avec h = max des hT et Q = max des QT
def pas_et_qualite_maillage(nbe, coord, tri):
    h = 0.0
    Q = 0.0

    for i in range(nbe):
        # Coordonnées des 3 sommets du triangle
        n1 = tri[i][0] - 1
        n2 = tri[i][1] - 1
        n3 = tri[i][2] - 1
        x1 = coord[n1][0]
        y1 = coord[n1][1]
        x2 = coord[n2][0]
        y2 = coord[n2][1]
        x3 = coord[n3][0]
        y3 = coord[n3][1]

        # Longueurs des 3 cotés
        a = sqrt((x2 - x3) ** 2 + (y2 - y3) ** 2)
        b = sqrt((x1 - x3) ** 2 + (y1 - y3) ** 2)
        c = sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

        # Aire du triangle (déterminant, en valeur absolue)
        aire = ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2
        if aire < 0:
            aire = -aire

        # Diamètre = plus grand coté
        hT = a
        if b > hT:
            hT = b
        if c > hT:
            hT = c

        rhoT = 2 * aire / (a + b + c)       # rayon du cercle inscrit
        QT = hT / (2 * sqrt(3) * rhoT)     # qualité du triangle

        # On garde le max sur tous les triangles
        if hT > h:
            h = hT
        if QT > Q:
            Q = QT

    return [h, Q]


## Assemble la matrice A et le second membre F du système A * Uh = F
## @param nbn, nbe, nba : nombre de noeuds, de triangles et d'aretes
## @param coord : coordonnées des noeuds [[x, y], ...]
## @param tri : sommets des triangles [[n1, n2, n3], ...] (numérotés à partir de 1)
## @param ar : sommets des aretes [[n1, n2], ...] (numérotés à partir de 1)
## @param refa : références des aretes
## @return : [A, F, K]
def assemblage_EF_P1(nbn, nbe, nba, coord, tri, ar, refa):
    A = [[0] * nbn for _ in range(nbn)]
    F = [0] * nbn

    for l in range(nbe):
        # Indice globaux
        I1 = tri[l][0] - 1
        I2 = tri[l][1] - 1
        I3 = tri[l][2] - 1

        xl = [coord[I1][0], coord[I2][0], coord[I3][0]]
        yl = [coord[I1][1], coord[I2][1], coord[I3][1]]

        kl = coeffelem_P1_rigid(xl, yl)
        fl = coeffelem_P1_source(xl, yl)

        # Ajoute kl et fl aux cases de A et F
        I = [I1, I2, I3]
        for i in range(3):
            for j in range(3):
                A[I[i]][I[j]] += kl[i][j]

            F[I[i]] += fl[i]

    # Copie de A avant la pénalisation
    K = [[0] * nbn for _ in range(nbn)]
    for i in range(nbn):
        for j in range(nbn):
            K[i][j] = A[i][j]

    for a in range(nba):
        if refa[a] == 1:
            # Indices globaux
            I1 = ar[a][0] - 1
            I2 = ar[a][1] - 1

            xa = [coord[I1][0], coord[I2][0]]
            ya = [coord[I1][1], coord[I2][1]]

            pa = coeffelem_P1_poids(xa, ya)
            ea = coeffelem_P1_transf(xa ,ya)

            # Ajoute pa et ea aux cases de A et F
            I = [I1, I2]
            for i in range(2):
                for j in range(2):
                    A[I[i]][I[j]] += pa[i][j]

                F[I[i]] += ea[i]

    return [A, F, K]


## Résouddre le système A * U = F avec Gauss
## @param A : matrice n x n (liste de listes)
## @param F : second membre (liste de n valeurs)
## @return : U, liste de n valeurs
def gauss(A, F):
    n = len(F)

    # Copies
    M = [ligne[:] for ligne in A]
    b = F[:]

    for k in range(n - 1):
        for i in range(k + 1, n):
            facteur = M[i][k] / M[k][k]

            # Rien à faire si la case est déjà nulle
            if facteur == 0:
                continue

            for j in range(k, n):
                M[i][j] -= facteur * M[k][j]

            b[i] -= facteur * b[k]

    # Remontée
    U = [0] * n
    for i in range(n - 1, -1, -1):
        somme = 0
        for j in range(i + 1, n):
            somme += M[i][j] * U[j]

        U[i] = (b[i] - somme) / M[i][i]

    return U


## Calcule l'erreur  (e_h = racine( (U - Uh)^T K (U - Uh) ))
## @param nbn : nombre de noeuds
## @param coord : coordonnées des noeuds [[x, y], ...]
## @param K : matrice de rigidité sans pénalisation
## @param Uh : solution calculée aux noeuds
## @return : e_h
def erreur_H1(nbn, coord, K, Uh):

    # Solution exacte aux noeuds
    U = [0.0] * nbn
    for i in range(nbn):
        U[i] = fct_u(coord[i][0], coord[i][1])

    # Ecart noeud par noeud
    E = [0.0] * nbn
    for i in range(nbn):
        E[i] = U[i] - Uh[i]

    # Produit matrice vecteur KE = K * E
    KE = [0.0] * nbn
    for i in range(nbn):
        for j in range(nbn):
            KE[i] += K[i][j] * E[j]

    # Produit scalaire E^T * KE
    total = 0.0
    for i in range(nbn):
        total += E[i] * KE[i]

    e_h = total ** 0.5

    return e_h


## Lit un maillage, résout le problème et calcule h, Q et l'erreur
## @param fichier_msh : chemin du fichier .msh
## @return : [h, Q, e_h]
def resout_maillage(fichier_msh):
    nbn, nbe, nba, coord, tri, ar, refn, reft, refa = lit_fichier_msh(fichier_msh)

    h, Q = pas_et_qualite_maillage(nbe, coord, tri)
    A, F, K = assemblage_EF_P1(nbn, nbe, nba, coord, tri, ar, refa)
    Uh = gauss(A, F)
    e_h = erreur_H1(nbn, coord, K, Uh)

    return [h, Q, e_h]

maillages = ["m00", "m1", "m2", "m3", "m4"]
liste_h = []
liste_eh = []

print("maillage |         h          |         Q          |        e_h")
for nom in maillages:
    h, Q, e_h = resout_maillage("./Maillages_M1-CHPS/" + nom + ".msh")
    liste_h.append(h)
    liste_eh.append(e_h)
    print(nom, "      |", h, "|", Q, "|", e_h)

print("")
for i in range(1, len(maillages)):
    p = log(liste_eh[i - 1] / liste_eh[i]) / log(2)
    print("p entre", maillages[i - 1], "et", maillages[i], "=", p)
