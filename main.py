from __future__ import division
import sys

from contig import Contig
from exon import Exon, ExonsList, NoMoreExons
from seq import Seq
from importer import Importer

class Interpreter(Seq):

    def __init__(self, maf_file, annotations_file, regions_names, coding=True):

        Seq.__init__(self)

        self.coding = coding
        annotations = Importer.get_lines_from_file(annotations_file)
        maf_file = Importer.get_lines_from_file(maf_file)


        for region_name in regions_names:
            print(region_name)
            self.tri_num = 0
            self.mutation_number = 0

            self.trinucleotides_occurences = self.generate_empty_params()
            self.neutral_matrix = self.generate_empty_params()

            regions = self.extract_regions_from_annotations(annotations, region_name)

            self.analyse_maf_file(maf_file, regions)

            self.export_tri_occ_and_mutations(region_name)
            self.normalize_neutral_matrix()
            self.export_results(region_name)


    def extract_regions_from_annotations(self, annotation_lines, region_name):
        ''' gets the annotation file lines
            returns a dictionary of exons positions, keys are the start position
            value is a list [start_pos, end_pos, chromosome,
            strand(- or + ), and the name of the gene]
            also returns the sorted keys of this dictionary (sorted start positions)
        '''
        exons = ExonsList()

        for idx, line in enumerate(annotation_lines):
            line = line.split('\t')
            if line[2] == region_name:
                exon = Exon(line)
                exons.add_exon(exon)

            self.print_progres(idx, len(annotation_lines))

        sys.stderr.write("\n")
        return exons

    def analyse_maf_file(self, maf_file, exons):
        ''' gets the maf_file lines, the dictionary of exons, and the exons keys 9start positions)
            list in order

            analysis the maf_file line by line, according to the exons from the annotations file
            analyse means generating counting the observed mutations per, and calculating the
            expected mutations per gene

            expected mutations are calculated by multiplying the number of occurence for each
            triblet and multiplying it by the chance of mutation from the neutral matrix
        '''
        file_len = len(maf_file)
        contig = Contig("hg38", "mrca") #creat an initial empty contig
        #exon = exons[start_pos_list.pop(0)] #start with the first exon

        for progress, line in enumerate(maf_file):

            if Contig.new_contig(line):
                contig = Contig("hg38", "mrca")
            else:
                contig.add_line(line)

            if contig.got_both_seq() and exons.not_empty(contig.chrom):

                exon = exons.get_first_exon(contig.chrom)

                if contig.is_before_exon(exon):
                    contig = Contig("hg38", "mrca")
                    continue #go for the next contig

                while contig.is_after_exon(exon):
                    try:
                        exon = exons.get_next_exon(contig.chrom)
                    except NoMoreExons:
                        break


                if contig.has_exon(exon):
                    self.analyse_seq(contig, exon)

            Interpreter.print_progres(progress, file_len)

    def analyse_seq(self, contig, exon):

        pointer = self.get_contig_pointer_by_exon(contig, exon)

        stop_point = min((exon.end_position - contig.strt_pos - 4),
                         (contig.length - 4))

        while pointer < stop_point:
            ref_codon = contig.ref_seq[pointer:pointer+3]
            pentamer = contig.ref_seq[pointer-1:pointer+4]
            alg_codon = contig.alg_seq[pointer:pointer+3]

            if "-" in ref_codon: #has insertion, unrechable now as we remove gaps first
                gap_index = pointer + ref_codon.index("-")
                contig.ref_seq.pop(gap_index)
                contig.alg_seq.pop(gap_index)
                stop_point -= 1
                continue

            else:
                base_ref_index = pointer + contig.strt_pos
                self.add_codon(ref_codon, alg_codon, pentamer, exon, base_ref_index)
                pointer += 3

    def add_codon(self, ref_codon, alg_codon, pentamer, exon, base_ref_index):

        if self.check_region_not_available(ref_codon, alg_codon, pentamer):
            return

        if exon.negative_strand: #negative strand
            ref_codon = self.get_reverse_complement(ref_codon)
            alg_codon = self.get_reverse_complement(alg_codon)
            pentamer = self.get_reverse_complement(pentamer)

        self.add_trinucleotides_occurences(pentamer)

        if alg_codon != ref_codon: #mutation
            self.add_codon_to_neutral_matrix(ref_codon, alg_codon, pentamer, exon)


    @classmethod
    def get_contig_pointer_by_exon(cls, contig, exon):

        if  exon.start_position > contig.strt_pos:
            start_pos = exon.start_position
        else:
            start_pos = contig.strt_pos + 3 #+3 to be changed

        if not exon.negative_strand: #positive strand
            while (start_pos - exon.start_position)%3 != 0: #so we don't start middle of a codon
                start_pos += 1
        else:
            while (exon.end_position - start_pos)%3 != 0: #so we don't start  middle of a codon
                start_pos += 1

        pointer = start_pos - contig.strt_pos #pointer of the seq bases one by one


        return pointer

    @classmethod
    def generate_empty_params(cls):
        nuclutides = ["A", "C", "G", "T"]
        neutral_matrix = {}
        for i in nuclutides:
            for j in nuclutides:
                for k in nuclutides:
                    neutral_matrix[i+j+k] = {"A":0, "C":0, "G":0, "T":0}
        return neutral_matrix

    def normalize_neutral_matrix(self):
        '''
        divides the number of mutations (context => base) by the number of occurences
        of this context where it was possible for that context to mutate
        '''
        for context, bases in self.neutral_matrix.items():
            for base, _value in bases.items():
                if self.neutral_matrix[context][base] != 0: #to avoid zero div
                    self.neutral_matrix[context][base] /= \
                    float(self.trinucleotides_occurences[context][base])

    def add_trinucleotides_occurences(self, pentamer):
        ''' calls the coding or the noncoding method to calculate trin occ'''
        if self.coding:
            self.add_trinucleotides_occurences_coding(pentamer)
        else:
            self.add_trinucleotides_occurences_noncoding(pentamer)

    def add_trinucleotides_occurences_noncoding(self, pentamer):
        '''all posiible mutation locations are added'''
        for i in range(3):
            codon = pentamer[i:i+3]
            for mutant_base in ['A', 'C', 'G', 'T']:
                self.trinucleotides_occurences[codon][mutant_base] += 1
                self.tri_num += 1


    def add_trinucleotides_occurences_coding(self, pentamer):
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
        self.tri_num += count

    def add_codon_to_neutral_matrix(self, ref_codon, alg_codon, pentamer, exon):
        """
        docstring
        """
        if (not self.coding) or self.get_mutation_type(ref_codon, alg_codon) == "coding-synon":
            self.mutation_number += 1

            if self.get_hamming_distance(ref_codon, alg_codon) == 1: #don't consider double mutants
                base, idx = self.get_mutant_base_and_its_index(ref_codon, alg_codon)
                trin = pentamer[idx: idx+3]
                self.neutral_matrix[trin][base] += 1

    def export_results(self, region_name):
        matrix_file_name = "neutral_matrices/" + region_name + "_normalized_matrix.txt"
        Importer.export_dict_to_tsv(matrix_file_name, self.neutral_matrix)


        logs_file_name = "neutral_matrices/" + region_name + "_logs.txt"
        logs = 'mutability:' + '\t' + str(self.mutation_number/float(self.tri_num)) + '\n'
        logs += 'mutations number:' + '\t' + str(self.mutation_number) + '\n'
        logs += 'trin number:' + '\t' + str(self.tri_num) + '\n'
        Importer.export_text(logs_file_name, logs)

    def export_tri_occ_and_mutations(self, region_name):
        prefix = "neutral_matrices/" + region_name + '_'
        Importer.export_dict_to_tsv(prefix + "tri_occ.txt", self.trinucleotides_occurences)
        Importer.export_dict_to_tsv(prefix + "mutations.txt", self.neutral_matrix)

    def check_region_not_available(self, *args):
        for arg in args:
            if "-" in arg or "X" in arg or "N" in arg:
                return True
        return False

Interpreter("mrca_mult.maf", "hg38Regionsannotations.gtf", ["Intergenic"], False)
