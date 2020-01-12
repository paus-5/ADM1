def load_parameters():
    with open("parameters.txt",'r') as file:

        data = file.read()
    data = data.split ('\n')

    dico = {}
    for element in data:

        elem = element.split (' ') # affiche le contenu de chaque ligne
        var = elem [0] # affiche le contenu de premier champs de la ligne 
       
        if len(elem)>1:
            valeur = elem [2] 
            valeur= float(valeur)
            dico [var] = valeur
        else:
            break

    return(dico)
dico=load_parameters()
#print(dico['k_dis'])ou (dico["k_dis"]) # pour afficher la valeur de k_dis
