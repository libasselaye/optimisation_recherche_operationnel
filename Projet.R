# ALGORITHME 1
Arbre_Couvrant = function(X,U){
  #i est le sommet que je l'ai s√©l√©ctionn√© ici le 6 
  i=6
  #Boolean ou je vais l'utiliser pour stoper m'as boucle .
  fin=TRUE
  #X' re√ßoit le sommet s√©lectionn√©e 
  X_prime <- c(i)
  #Matrice √©gale a la matrice d'adjacence.
  A_prime= matrix(100,nrow =7,ncol= 7)
  #U' re√ßoit un tab vide
  U_prime<-c()
  #tant que la longeur de x est diff de x' cela veut dire que j'ai pas encore visit√© tout les sommets
  while (length(X_prime)!=length(X))
  {
    while (fin) {
      for(j in 1:nrow(A)){
        #si le sommet exisite dans x' je cherche son voisin dans la matrice d'adjacence
        if((is.element(i,X_prime)) & !(is.element(j,X_prime)) ){
          #si le sommet j est un voisin 
          if(A[i,j]>0){
            #je sauvgarde cette valeur dans la matrice que je l'ai cr√©er 
            A_prime[i,j]=A[i,j]
          }
        }
      }
      #je cherche le minimum dans la matrice 
      x<-which(A_prime==min(A_prime),arr.ind=T)
      #si cette element n'est pas dans x'
      if(!(is.element(x[1,2],X_prime)))
      {
        #je sauvgarde ces coordonn√©es. dans U' et X'
        U_prime<-c(U_prime,x[1,])
        X_prime<-c(X_prime,x[1,2])
        #apr√®s le i devient le sommet s√©l√©ctionn√©
        i=x[1,2]
        
        A_prime[x[1],x[2]]=100
      }
      #Condition d'arret de ma boucle si j'ai visit√© tout les sommets de ce graphe.
      if(length(X)==length(X_prime)){
        fin=FALSE
      }
    }
  }
  #Seulement pour voir mes r√©sultat j'ai pas retourner j'ai afficher .
  print("U PRIME")
  print(U_prime)
  print("X PRIME")
  print(X_prime)
}

#ALGORITHME 2
Ford_Belleman = function(X,A,s){
  #Algorithme de plus court chemin de Moore Dijkstra (poids non nÈgatifs)
  #INPUT
  #X est l'ensemble des sommets
  #A est la matrice d'adjacence pondÈrÈe
  #s le sommet initial
  #OUTPUT
  #pi vecteur donnant la longeur du ppc entre s et les autres sommets
  
  #Longueur du vecteur X 
  n=length(X)
  pi=rep(0,n)
  pi[s]=0
  Sb=setdiff(X,s)#initialisation de Sb
  
  for (i in Sb)
  {
    pi[i] = Inf
  }
  
  #initialisation de la variable permettant de stocker 
  #la derniÈre valeur de pi
  stock = NULL
  
  #boucle repeat until 
  repeat{
    for (i in Sb)
    {
      #calcul du plus court chemin
      pi[i] = min(pi[i],pi[which(A[,i]!=0)]+A[which(A[,i]!=0),i])
    }
    
    #on sort de la boucle si la derniÈre valeur de pi ne change pas 
    if (setequal(pi,stock)==TRUE) break
    stock = pi
  }
  return(pi)
}

#ALGORITHME 3
FordFulkerson = function(X,A,s,p)
{
  
  phi = 0 #flot realisable
  #Preparation des triplets de chaque noeuds
  capacites = List(X,s)
  chaine = paste("m",s,sep="")
  capacites[chaine,2] = 9999
  capacites[chaine,3] = "+"
  #Initialisation d'une matrice P qui va contenir le flot qui passe entre deux noeuds
  P = matrix(0, nrow = length(X), ncol = length(X))
  P[s,s]=9999 #flot qui part de la source, vaut :  +l'infini
  while(TRUE)
  {
    S = c(s)#Liste des noeuds marquÈs, s marquÈ dËs le dÈpart de chaque itÈration
    Sb = setdiff(X,s)#Liste des noeuds non marquÈs
    tmp = A-P > 0 #Test Cij - PHIij > 0 de l'algorithme
    tmp2 = t(P) >0#Recuperation des PHI non nuls
    C= tmp | tmp2
    index = which(matrix(C[S,Sb] ==TRUE,nrow=length(S), ncol=length(Sb)),arr.ind = TRUE)
    index2 = which(matrix(P[Sb,S] > 0, nrow=length(Sb), ncol = length(S)), arr.ind = TRUE)
    
    while(length(index)>0 || length(index2)>0)#Ligne 2 algo
    {
      if(length(index)>0)
      {#ligne 3
        i = index[1,1]
        i = S[i]
        j = index[1,2]
        j = Sb[j]
        #Construction de la chaine de caracteres pour correspondre a la ligne du triplet mi et mj
        chaine = paste("m",i,sep="")
        chaine2 = paste("m",j,sep="")
        capacites[chaine2,1] = i
        capacites[chaine2,2] = min(capacites[chaine,2],A[i,j]-P[i,j])
        capacites[chaine2,3] = "+"
        #mj=(i,min(capacites[chaine,2],A[i,j]-P[i,j]),"+")
      }
      else
      {#ligne 5
        if(length(index2)>0)
        {
          i = index2[1,1]
          i = S[i]
          j = index2[1,2]
          j = Sb[j]
          #On met a jour le triplet mj dans le dataframe
          chaine=paste("m",i,sep="")
          chaine2 = paste("m",j,sep="")
          capacites[chaine2,1] = i
          capacites[chaine2,2] = min(capacites[chaine,2],P[j,i])
          capacites[chaine2,3] = "-"
          #mj=(i,min(capacites[chaine,2],P[j,i]),"-")
        }
      }
      #On met dans S le sommet qui vient d'etre mis a jour(on le marque)
      S = cbind(S,j)#ligne 8
      Sb = setdiff(X,S)
      #Est ce que j est le puit?
      if(j==p)
      {#oui, on met a jour la valeur du flot maximal et on arrete la boucle courrante
        phi = phi + capacites[chaine2,2]
        phiP = capacites[chaine2,2]
        break
      }
      #On met a jour l'index
      tmp = A-P > 0
      tmp2 = t(P)
      C= tmp | tmp2
      index = which(matrix(C[S,Sb] ==TRUE,nrow=length(S), ncol=length(Sb)),arr.ind = TRUE)
      index2 = which(matrix(P[Sb,S] > 0, nrow=length(Sb), ncol = length(S)), arr.ind = TRUE)
    }
    #Le puit est il marque?
    if(is.element(p,S))
    {#ligne 14
      repeat
      {#ligne 15
        if(j==s)
        {
          break
        }
        chaine = paste("m",j,sep="")
        if(capacites[chaine,3]=="+")
        {#ligne 16
          i = capacites[chaine,1]
          P[i,j] = P[i,j]+phiP#Mise a jour de P aux indices des sommets qui ont permis d'atteindre le puit en additionnant la valeur du flot
        }
        else
        {
          if(capacites[chaine,3]=="-")
          {
            i = capacites[chaine,1]
            P[j,i] = P[j,i]-phiP#Mise a jour de P aux indices des sommets qui ont permis d'atteindre le puit en soustrayant la valeur du flot
          }
        }
        j = capacites[chaine,1]#ligne 21
      }
    }
    else
    {
      break
    }
    
  }
  return(phi)
}