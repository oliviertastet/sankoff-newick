# -*- coding: utf-8 -*-
import math
import argparse
import sys

# La classe Node instancie les noeuds d'un arbre binaire. Il est binaire par définition car un noeud ne peux avoir que 2 enfants 
# La méthode add_child permet d'ajouter un enfant à un noeud de manière à ce que l'arbre soit balancé. Si on demande l'ajout de 3 
# enfants, l'assignation de l'enfant ne se fait pas et un message d'erreur est retourné. La classe contient également plusieurs 
# méthodes d'accession et de fixation
class Node():
    def __init__(self,label):
        self.right = None
        self.left = None 
        self.parent = None
        self.label = label
        self.sequence = None
        self.score= 0
        self.scores = [float('inf')]*23

    def add_child(self, node):
        if self.left == None:
            self.left = node
            node.parent = self
        elif self.right == None:
            self.right = node
            node.parent = self
        else:
            return 0

    def get_label(self):
        return self.label

    def set_node_sequence(self, sequence):
        self.sequence = sequence

    def get_character_in_sequence(self, i):
        return self.sequence[i]

    def set_node_score(self, score, i):
        self.scores[i] = score

    def get_node_alphabet_scores(self):
        return self.scores

    def set_node_alphabet_scores(self, scores):
        self.scores = scores

    def get_alphabet_index_from_score(self, value):
        return self.scores.index(value)

    def get_right_child(self):
        return self.right

    def get_left_child(self):
        return self.left

    def get_parent(self):
        return self.parent

    def set_parent(self, node):
        self.parent = node
