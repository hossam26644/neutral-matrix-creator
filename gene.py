import pandas as pd
from seq import Seq, BASES, CONTEXTS

PENTS = [x+y+c for x in BASES for y in BASES for c in CONTEXTS]
MUTTYPES = ['coding-synon', 'missense', 'nonsense']

class Gene(Seq):

    def __init__(self, name, coding=True):
        Seq.__init__(self)
        self.length = 0
        self.name = name
        self.coding = coding
        self.pentamer_occ = pd.Series(0, index=PENTS)
        '''pentamers are not all the pentamers, they are the codons with a base before
           and a base after'''
        self.exp_ops = pd.DataFrame(0, index=MUTTYPES, columns=['exp', 'obs'])
        self.seq = ''

    def add_pentamer(self, pentamer, alg_codon):
        self.pentamer_occ[pentamer] += 1
        self.length += 3
        ref_codon = pentamer[1:4]
        self.seq += self.codon_table[alg_codon]
        if ref_codon != alg_codon:
            try: self.change_in_dtime.append(self.dtime[ref_codon][0]-self.dtime[alg_codon][0])
            except: pass
            mutation_type = self.get_mutation_type(ref_codon, alg_codon) if self.coding else "missense"
            self.exp_ops['obs'][mutation_type] += 1

    def calculate_exp_mutations(self, neutral_matrix, pent_matrix=None):
        pent_matrix = pent_matrix if pent_matrix is not None else self.expected_mutations_per_pentamer(neutral_matrix)
        for mut_type in MUTTYPES:
            self.exp_ops['exp'][mut_type] += sum(self.pentamer_occ * pent_matrix[mut_type])

    def expected_mutations_per_pentamer(self, neutral_matrix):
        ''' gets a tri neutral matrix, returns the equiv of mutations by type for pentamers
            ex from AAA->C:0.05 => AAAAA->missense:5, nonsense:3
        '''
        pent_mutations = pd.DataFrame(0, index=PENTS, columns=MUTTYPES)

        for codon in CONTEXTS:
            for i in range(3):
                for mut_b in BASES:
                    mut_codon = codon[:i] + mut_b + codon[i+1:]
                    if codon == mut_codon: continue
                    pentamers = [pent for pent in PENTS if pent[1:4]==codon]
                    contexts = {pentamer[i:i+3] for pentamer in pentamers}
                    contexts_pentamers = [(pentamer[i:i+3], pentamer) for pentamer in pentamers]
                    mutation_type = "missense" if not self.coding else self.get_mutation_type(codon, mut_codon)
                    for contexts_pentamer in contexts_pentamers:
                        context, pentamer = contexts_pentamer
                        pent_mutations.at[pentamer, mutation_type] += neutral_matrix[context][mut_b]

        return pent_mutations

class GenesList():

    def __init__(self):
        self.genes = {}

    def has_gene(self, gene_name):
        return gene_name in self.genes

    def add_gene(self, gene_name):
        if self.has_gene(gene_name): raise KeyError('Gene exists')
        self.genes[gene_name] = Gene(gene_name)

    def get_gene(self, gene_name):
        if not self.has_gene(gene_name): raise KeyError('Gene does not exist')
        return self.genes[gene_name]

    def add_pentamer(self, gene_name, pentamer, alg_codon):
        try:
            gene = self.get_gene(gene_name)
        except KeyError:
            self.add_gene(gene_name)
            gene = self.get_gene(gene_name)
        gene.add_pentamer(pentamer, alg_codon)

    def calculate_exp_mutations(self, neutral_matrix):
        pent_matrix = None
        for gene in self.genes.values():
            pent_matrix = pent_matrix if pent_matrix is not None else gene.expected_mutations_per_pentamer(neutral_matrix)
            gene.calculate_exp_mutations(neutral_matrix, pent_matrix)

    def clear_occ(self):
        for gene in self.genes.values():
            gene.pentamer_occ = pd.Series(0, index=PENTS)


    def __iter__(self):
        yield list(self.genes.values()())
