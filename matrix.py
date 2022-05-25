import pandas as pd
from seq import Seq, BASES, CONTEXTS


class Matrix(Seq):

    def __init__(self, coding):
        Seq.__init__(self)
        self.coding = coding
        self.trinucleotides_occurences = pd.DataFrame(0, index=BASES, columns=CONTEXTS)
        self.mutations = pd.DataFrame(0, index=BASES, columns=CONTEXTS)

    def add_tri_occ_noncoding(self, pentamer):
        '''all posiible mutation locations are added'''
        for i in range(3):
            codon = pentamer[i:i+3]
            for mutant_base in ['A', 'C', 'G', 'T']:
                self.trinucleotides_occurences[codon][mutant_base] += 1

    def add_tri_occ_coding(self, pentamer):
        '''
        gets a pentamer, adds the trinucleotides occurences to neutral matrix
        based on their possibility to mutate synonymsly
        '''
        codon = pentamer[1:-1]
        count = 0
        for variant in self.syn_variants_dict[codon]:
            if self.get_hamming_distance(codon, variant) == 1: #don't consider double mutants
                count = 1
                mutant_base, mutant_index = self.get_mutant_base_and_its_index(codon, variant)
                context = pentamer[mutant_index:mutant_index+3]
                self.trinucleotides_occurences[context][mutant_base] += 1

    def add_mutation(self, ref_codon, alg_codon, pentamer):
        """
        docstring
        """
        if (not self.coding) or self.get_mutation_type(ref_codon, alg_codon) == "coding-synon":

            if self.get_hamming_distance(ref_codon, alg_codon) != 1: return #don't consider double mutants

            base, idx = self.get_mutant_base_and_its_index(ref_codon, alg_codon)
            trin = pentamer[idx: idx+3]
            self.mutations[trin][base] += 1

    def get_mutational_matrix(self):
        matrix = self.mutations/self.trinucleotides_occurences
        return matrix.fillna(0)
