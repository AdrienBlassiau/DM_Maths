from numpy import*


def verifcompletregul(M):
    '''Fonction qui vérifie si une matrice est complètement régulière ou non'''
    
    t=M.shape[0] # On prend la taille de la matrice M.
    
    for i in range(1,t+1):
        deter=linalg.det(M[:i,:i]) # On teste l'inversibilité des matrices extraites de M de taille (i,i) avec 1<=i<=t.
        if deter==0.0: # Si le déterminant de la matrice extraite est nul, elle n'est pas inversible est M n'est pas complétement régulière.
            return False
            
    return True  # La matrice M a passé le test, elle est validée !
    
    

def inverse(M):
    '''Fonction qui permet de calculer l'inverse de toute matrice complétement régulière'''

    # Etape 0 : on prend la taille de la matrice M.
    
    t=M.shape[0]
    
    # Etape 1 : on vérifie que l'inverse n'est pas trivial.
    
    if t==1 or t==2: # Si la matrice est de taille 1 ou 2, cas de base.
        return linalg.inv(M)  
        
        
    # Etape 2 : on découpe la matrice M en 4 blocs : A (Nord-Ouest),B (Nord-Est), C (Sud-Ouest) et D(Sud-Est).
    
    l=int(t//2)   # On définit le plan de coupe de manière à avoir 4 blocs sensiblement de même taille.


    A=M[:l,:l] # Bloc NORD-OUEST de taille (l,l).
    B=M[:l,l:] # Bloc NORD-EST de taille (l,t-l).
    C=M[l:,:l] # Bloc SUD-OUEST de taille (t-l,l).
    D=M[l:,l:] # Bloc SUD-EST de taille (t-l,t-l).

    
    # Etape 3 : on inverse la matrice A.
    
    InvA=inverse(A)
    
    # Etape 4 : on définit le complément de Schur et on l'inverse.
    
    S=dot(dot(-C,InvA),B)+D
    InvS=inverse(S)

    # Etape 5 : on définit les 4 blocs de la matrice inversée en utilisant la formule du Dm trouvée en 3. c).
    
    BlocNO=InvA+dot(dot(InvA,B),dot(dot(InvS,C),InvA)) # Bloc NORD-OUEST de taille (l,l).
    BlocNE=dot(dot(-InvA,B),InvS) # Bloc NORD-EST de la matrice M inversée de taille (l,t-l).
    BlocSO=dot(dot(-InvS,C),InvA) # Bloc SUD-OUEST de la matrice M inversée de taille (t-l,l).
    BlocSE=InvS # Bloc SUD-EST de la matrice M inversée de taille (t-l,t-l).
    
    # Etape 6 : on réunit les 4 blocs pour avoir l'inverse de M.
    
    InvM=hstack((vstack((BlocNO,BlocSO)),vstack((BlocNE,BlocSE))))
    
    return InvM
   
   
   
   
## QUESTION 2   
    

M=random.randint(100,size=(13,13)) # On génère aléatoirement une matrice .   

print(M)

test=verifcompletregul(M)  # On vérifie que M est complètement régulière :

if test==True:# Si oui on passe au calcul de l'inverse ...
    print("\n\n","La matrice M est complètement régulière.\n\n")
    
    InvMRecur=inverse(M) # On calcule l'inverse avec la méthode proposée.
    print("L'inverse de M par récursivité est :\n\n",InvMRecur,"\n\n")
    
    InvM=linalg.inv(M)   # On calcule l'inverse avec la fonction dédiée pour comparer les résultats.
    print("L'inverse de M par linalg.inv(M) est :\n\n",InvM,"\n\n")
    
else: # ...sinon le problème s'arrête là.
    print("La matrice M n'est pas complètement régulière\n\n")
    

## QUESTION 3   



#On vérifie que M est inversible et que la transposée de M multipliée par M est complètement régulière :
test1=linalg.det(M)
test2=verifcompletregul(dot(M.T,M))  

if test1!=0 and test2==True:# Si oui on passe au calcul de l'inverse ...
    print("La transposée de M multipliée par M est complètement régulière.\n\n")
    
    InvMTMRecur=inverse(dot(M.T,M))  # On calcule l'inverse avec la méthode proposée.
    InvMRecu=dot(InvMTMRecur,M.T) # Pour obtenir M on multiplie l'inverse de (transposée de M)*M  par la transposée de M. (voir copie).
    print("\n\n","L'inverse de M par récursivité est :\n\n",InvMRecur,"\n\n")
    
else: # ...sinon le problème s'arrête là.
    print("La matrice M n'est pas complètement régulière")
    







