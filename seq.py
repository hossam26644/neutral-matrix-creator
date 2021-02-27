import sys
from importer import Importer


class Seq(object):
    """
    docstring
    """
    def __init__(self):
        self.codon_table = Importer.parse_json("Auxiliary/codon_table.json")
        self.mutations_count = {"coding-synon":0, "nonsense":0, "missense":0}
        self.possible_syn_counts_table = self.generate_possible_syn_mutations_counts_table()
        self.syn_variants_dict = self.get_syn_variants_dict()
        self.bases = ["A", "C", "G", "T"]

    @classmethod
    def get_reverse_complement(cls, sequence):
        complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        rev_comp = ""
        for base in reversed(sequence):
            rev_comp += complement[base]

        return rev_comp

    def get_mutation_type(self, ref_codon, alg_codon):

        if self.codon_table[alg_codon] == "Stop" or  self.codon_table[ref_codon] == "Stop":
            return "nonsense"

        elif self.codon_table[ref_codon] == self.codon_table[alg_codon]:
            return "coding-synon"

        return "missense"

    def get_syn_variants_dict(self):
        """
        returns a dict of codons as keys, and codons that encode the same as
        amino acid as values
        """
        syn_variants = {}
        for codon, amino_acid in self.codon_table.items():
            syn_codons = []
            for target_codon, target_amino_acid in self.codon_table.items():
                if codon != target_codon and \
                amino_acid == target_amino_acid and \
                amino_acid != "Stop":
                    syn_codons.append(target_codon)

            syn_variants[codon] = syn_codons

        return syn_variants

    @staticmethod
    def print_progres(current_progress, max_progress):
        """
        docstring
        """
        if current_progress%1000 == 0:
            sys.stderr.write("%.f%% done.\r" %(100.*current_progress/max_progress))

    def generate_possible_syn_mutations_counts_table(self):
        """
        docstring
        """
        nuclutides = ["A", "C", "G", "T"]
        possible_syn_counts_table = {}

        for codon in self.codon_table.keys():
            possible_syn_counts_table[codon] = [0, 0, 0]

            for idx, base in enumerate(codon):
                for mut_base in nuclutides:
                    mut_codon = list(codon)
                    mut_codon[idx] = mut_base
                    mut_codon = "".join(mut_codon)
                    is_syn = self.get_mutation_type(codon, mut_codon) == "coding-synon"
                    if is_syn and base != mut_base:
                        possible_syn_counts_table[codon][idx] += 1
        return possible_syn_counts_table

    @staticmethod
    def get_hamming_distance(frst_seq, scnd_seq):
        return sum(base1 != base2 for base1, base2 in zip(frst_seq, scnd_seq))

    @staticmethod
    def get_mutant_base_and_its_index(first_seq, second_seq):
        """
        gets two sequences, returns the first mutant base and it's indexed
        """
        first_seq = list(first_seq)
        second_seq = list(second_seq)

        for idx, base in enumerate(first_seq):
            if base != second_seq[idx]:
                return second_seq[idx], idx

        raise ValueError("no mutation")

    @staticmethod
    def get_triplets_contexts_from_pentamer(pentamer):
        """
        gets a pentamer returns a list of three contexts that can be
        made by this pentamer
        """
        contexts = []
        for i in range(len(pentamer)-2):
            contexts.append(pentamer[i:i+3])

        return contexts

    @staticmethod
    def get_lines_from_file(filename):
        fin = open(filename)
        lines = fin.readlines()
        fin.close()
        return lines

    @staticmethod
    def get_row_col_by_pos(pos=int):
        row = 1 + (pos//50)
        col = (pos%50) - 1
        if col == -1:
            row -= 1
            col = 49
        return row, col

    @staticmethod
    def get_base_by_pos(ref, pos=int):
        row, col = Seq.get_row_col_by_pos(pos)
        return ref[row][col]

    @staticmethod
    def get_context_by_pos(ref, pos=int):
        base0 = Seq.get_base_by_pos(ref, pos-1).upper()
        base1 = Seq.get_base_by_pos(ref, pos).upper()
        base2 = Seq.get_base_by_pos(ref, pos+1).upper()
        return ([base0, base1, base2])
