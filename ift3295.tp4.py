# -*- coding: utf-8 -*-
import math
import argparse
import sys
from node import Node

def build_sequence_dict(seqFile):   
    sequence, label = '',''
    sequenceDict = {}
    first = True
    fastaFile = open(seqFile,'r')
    for line in fastaFile:
        if line.startswith(">"):
            if first:
                label = line.strip('>').strip('\n')
                first = False
            if line is not '' and label is not '':
                sequenceDict[label] = sequence
                label = line.strip('>').strip('\n')
                sequence =''
        else:
            sequence+=line.strip("\n")
    sequenceDict[label] = sequence
    return sequenceDict

# La methode parse_newick_string prend une string representant un arbre binaire sous le format newick ainsi qu'un dictionnaire des 
# sequence correspondant aux feuilles de l'arbre. Les cles du dictionnaire et les etiquettes des feuilles de l'arbre 
# doivent concorder. La methode retourne un noeud qui est en fait la racine de l'arbre correspond a la sequence parsee
def parse_newick_string(s, seqDict):
    newickTree, leafLabel = s.strip(';'), ''
    nodeID, first = 0, True
    for character in newickTree:
        if first:
            root = Node(0)
            currentNode = root
            first = False
        elif character == '(':
            nodeID+=1
            node = Node(nodeID)
            node.set_parent(currentNode)
            currentNode.add_child(node)
            currentNode = node
        elif character == ',':
            if leafLabel != '':
                nodeID+=1
                newNode = Node(leafLabel.strip())
                newNode.set_node_sequence(seqDict[leafLabel.strip()])
                currentNode.add_child(newNode)
                leafLabel = ''
            if currentNode.get_right_child() is not None and currentNode.get_left_child() is not None:
                currentNode = currentNode.get_parent()
        elif character == ')':           
            if leafLabel != '':
                nodeID+=1
                newNode = Node(leafLabel.strip(' '))
                newNode.set_node_sequence(seqDict[leafLabel.strip(' ')])
                currentNode.add_child(newNode)
                if currentNode is 0:
                    return 0
                leafLabel = ''
            if currentNode.get_right_child()    is not None and currentNode.get_left_child() is not None:
                currentNode = currentNode.get_parent()
        else:
            leafLabel+=character
    return root

# La methode build_mutation_dict construit un dictionnaire qui represente une matrice de cout de substitution entre 
# acides aminés. La matrice doit etre contenue dans le fichier mutfile et avoir le format suivant:
#
#                         A B C
#                       A 0 2 3
#                       B 2 0 4
#                       C 3 4 0
#
# Le dictionnaire construit aura la forme suivante:
#   {
#         'A': 
#         {
#             'A': 0,
#             'B': 2,
#             'C': 3
#         },
#         ...
#         ,
#         'C':
#         {
#             'A': 3,
#             'B': 4,
#             'C': 0
#         }
# }
#
#
# Ainsi, selon un dictionnaire d, le cout de substitution entre X et Y est donné par d['X']['Y']
def build_mutation_dict(mutfile):           
    mutDict, count, first = {}, 0, True
    f = open(mutfile,'r')
    for lines in f:
        if first:
            first = False
            alphabet = lines.strip('\n').split()
            for character in alphabet:
                mutDict[character] = {}
        else:
            line = lines.strip('\n').split()
            for i in range(1, len(line)):
                mutDict[alphabet[count]][alphabet[i-1]] = int(line[i])
            count+=1
    return mutDict 

# La methode recursive get_tree_score retourne le score de la racine du sous arbre ou chaque feuille porte le ie caractere de la sequence
# Vu que la recursion de Sankoff determine le score d'un noeud en fonction du score de ses enfants, c'est un parcours post-fixe qui nous 
# permet d'identifier le score d'un arbre pour une position de la sequence donnee  
def get_tree_score(node, i, mutDict):
    alphabet, score = mutDict.keys(), 0
    if node is not None:
        get_tree_score(node.get_left_child(), i, mutDict)
        get_tree_score(node.get_right_child(), i, mutDict)        
        potentialScore = (get_min_score(node, i, mutDict))
        if node.get_parent() is None:
            return potentialScore
    
# La methode get_min_score permet de d'obtenir le score d'un noeud en fonction de ses enfants. Cette methode doit etre appelee dans 
# un parcours post-fixe car le score des enfants d'un noeud doivent etre fixe avant d'etablir son propre score car il depend de ces valeurs 
def get_min_score(node, pos, mutations):
    alphabet, score, discard= mutations.keys(), 0, False
    # Si les deux enfants d'un noeud sont null alors ce noeud est une feuille
    if node.get_left_child() is None and node.get_right_child() is None:
        x = node.get_character_in_sequence(pos)
        if x == '-':
            discard = True
        for i in range(0, len(alphabet)):
            if x == alphabet[i]:
                node.set_node_score(0, i )
            else:
                node.set_node_score(float('inf'), i)
        return
    else:

        # Remplir la matrice de score pour chacune des lettres de l'alphabet
        # Ensuite, choisir le minimum
        for j in range(0, len(alphabet)):
            minLeft = min(node.get_left_child().get_node_alphabet_scores())
            minRight = min(node.get_right_child().get_node_alphabet_scores())

            # Si un des caracteres des sequences est un '-' alors on ne considere par l'arbre sous-jacent
            if math.isinf(minRight) or math.isinf(minRight):
                discard = True
                break
            penaltyLeft = int(mutations[alphabet[node.get_left_child().get_alphabet_index_from_score(minLeft)]][alphabet[j]])
            penaltyRight = int(mutations[alphabet[node.get_right_child().get_alphabet_index_from_score(minRight)]][alphabet[j]])
            substitutionCost = minLeft + minRight + penaltyLeft + penaltyRight
            node.set_node_score(substitutionCost ,j)
        score=min(node.get_node_alphabet_scores())
    node.get_right_child().set_node_alphabet_scores([float('inf')]*23)   
    node.get_left_child().set_node_alphabet_scores([float('inf')]*23)
    if discard:
        return float('inf')
    else:
        return score

# Methode pour imprimer tout les noeuds de l'arbre. L'ensemble des enfants gauche d'un noeud est imprimer, ensuite le noeud 
# est imprimer puis les enfants droits du noeud sont imprimer
def print_tree(tree):
    if tree is not None:
        print_tree(tree.get_left_child())
        print tree.get_label()
        print_tree(tree.get_right_child())


# Methode d'execution des questions du TP
def main(treeFile, seqFile, mutFile):
    treeScores, count = [], 0
    seqDict, mutationDict = build_sequence_dict(seqFile), build_mutation_dict(mutFile)
    trees = open(treeFile,'r')
    for tree in trees:
        score, discardedTrees = 0, 0
        count += 1
        newickTree = tree.strip('\n')
        root = parse_newick_string(newickTree,seqDict)
        if root is 0:
            print 'L\'arbre %d n\'est pas en format binaire!!'%(count)
            continue
        sequence = seqDict[seqDict.keys()[0]]
        for i in range(0, len(sequence)):
            charTreeScore = get_tree_score(root, i, mutationDict)
            if not math.isinf(charTreeScore):
                score += charTreeScore
            else:
                discardedTrees+=1
        treeScores.append(score)
        print 'Arbre ' + str(count) + ' : ' + str(score)
    print '%d Arbres rejetés à cause de gaps'%(discardedTrees)
    print 'L\'Arbre le plus parcimonieux est le numéro ' + str(treeScores.index(min(treeScores))+1) + '\n avec un score de ' + str(min(treeScores))

# La méthode main n'est executée que si le programme est appelé directement et non importé en tant que package 
# argparse permet la gestion organisée des arguments
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Ce programme converti des arbres binaires en format newick en objets Node où un arbre est représenté par son noeud racine. '
        + 'Après avoir assigné chaque sequences aux feuilles correspondantes, l\'algorithme de Sankoff est appliqué sur l\'arbre. '
        + 'La sortie du programme montre le score de chaque arbre, le nombre d\'arbre(s) rejeté(s) par arbre et l\'arbre le plus parcimonieux. '
        + 'Sankoff calcule le score pour chaque caractère des séquences. Un arbre de caractère est rejeté lorsqu\'il contient un gap. '
        + 'Le même nombre d\'arbre(s) de caractère est rejeté pour l\'ensemble des arbres newick car ils représentent le même alignement. '
        )
    parser.add_argument(
                            '-t',
                            '--newick_file', 
                            help='Le nom du fichier qui contient les arbres binaire en format newick\n\tex:((A,B));', 
                            required=True
                        )
    parser.add_argument(
                            '-s',
                            '--sequence_file', 
                            help='Le nom du fichier en format fasta qui contient l\'ensemble des séquences et leurs'
                            + ' étiquettes correspondants aux feuilles des arbres du fichier newick', 
                            required=True
                        )
    parser.add_argument(
                            '-m',
                            '--mutation_file', 
                            help='Le nom du fichier contenant la matrice de couts de substitution pour les différents acides aminés', 
                            required=True
                        )
    # Si aucun argument n'est présent sauf le nom du programme, imprimer le message d'aide et indiquer la nature obligatoire des arguments
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    args=vars(parser.parse_args())
    main(args['newick_file'], args['sequence_file'], args['mutation_file'])

