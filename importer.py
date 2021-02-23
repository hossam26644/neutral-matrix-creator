''' used by the parameters class to import the static data
'''
import json
import sys
from blist import blist


class Importer(object):

    def __init__(self, filepath):
        self.filepath = "./" + filepath + "/"

    def import_quintuplets(self):
        ''' quiruples are the pentamers context
        '''
        filename = self.filepath + "quintuplets.txt"
        fin = open(filename)
        lines = fin.readlines()
        fin.close()
        return [line.strip().split()[0] for line in lines]

    def import_special_genes(self, filename):
        ''' special genes are cancer genes, essential genes and zero genes
            the assumption is that they won't have neutral background
        '''
        filename = self.filepath + filename
        fin = open(filename)
        lines = fin.readlines()
        fin.close()

        c_genes = set()
        for line in lines:
            c_genes.add(line.strip().split()[0])
        return c_genes

    def import_context_freqs(self, filename):
        ''' the precalculated triplets occurence frequencies
        '''
        filename = self.filepath + filename
        fin = open(filename)
        lines = fin.readlines()
        fin.close()
        occs = []

        for line in lines:
            field = line.strip().split()
            occs.append(float(field[1]))

        return occs

    def import_maf_data(self, filename, context_mode, triplets_map):
        '''	(1)	Import mutation annotation file (maf) including header line, "context" is 0-based.
            Format: ["gene", "muttype", "mutbase", "context"]

            gets the input filename (.cbase), context: 0for triplets 1 for pentamers, and the map between the new and old
            triplets code

            returns the mutation array: array of dictionaries, one per each input file line
            keys are the names of the fileds, attributes are gene name, mutation type, and the context
        '''
        fin = open("./" + filename)
        lines = fin.readlines()
        fin.close()
        mut_array = []
        for line in lines[1:]:
            field = line.strip().split("\t")
            if len(field) != 4:
                sys.stderr.write("Number of columns in maf file not as expected (=4): %i.\n" %len(field))
                sys.exit()
            # context is 0-based; self.triplets mapped to legacy
            if context_mode == 0:
                mut_array.append({"gene": field[0], "muttype": field[1], "mutbase": field[2].upper(), "context": triplets_map[int(field[3])]})
            else:
                mut_array.append({"gene": field[0], "muttype": field[1], "mutbase": field[2].upper(), "context": int(field[3])})

        return mut_array

    def import_codons_by_gene(self, filename):
        ''' the input file has the structure of: headers => ex: gene   A3GALT2
            and a context => ex: 15 50 46
            context is the trinuclutide context of each of the three bases of the codon

            returns a list of dictionaries, each has the key-value pairs of:
            "gene":gene name and "context_triplets":list of lists, each is a three values one,
            representing the context of each of the three bases of the codon

        '''
        filename = self.filepath + filename
        sys.stderr.write("Importing codons...\n")
        fin = open(filename)
        lines = blist(fin.readlines())
        fin.close()
        c_genes = blist(); cur_gene = "bla"; cur_cods = blist()

        for idx, line in enumerate(lines):
            field = line.strip().split("\t")
            if len(field) == 2: #header
                c_genes.append({"gene": cur_gene, "context_triplets": cur_cods})
                cur_cods = []
                cur_gene = field[1]
            else:
                cur_cods.append([int(el) for el in field])

            if idx%1000 == 0:
                sys.stderr.write("%.f%% done.\r" %(100.*idx/len(lines)))

        c_genes.append({"gene": cur_gene, "context_triplets": cur_cods})
        sys.stderr.write("\n")
        return c_genes[1:]

    def import_muttype_index(self):
        filename = self.filepath + "muttype_index.json"
        muttype_index = self.parse_json(filename)
        muttype_index = {int(k):v for k, v in muttype_index.items()}

        return muttype_index

    @staticmethod
    def get_lines_from_file(filename):
        fin = open(filename)
        lines = fin.readlines()
        fin.close()
        return lines

    @staticmethod
    def parse_json(filename):
        with open(filename, 'r') as fp:
            payload = json.load(fp)
        return payload
