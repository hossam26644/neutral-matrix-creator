from __future__ import division
import sys
import os

from contig import Contig
from exon import Exon, ExonsList, NoMoreExons #TODO change name to region
from seq import Seq
from importer import Importer
from excludedregions import ExcludedRegions
from matrix import Matrix
from gene import GenesList

class Interpreter(Seq):

    def __init__(self, maf_file, annotations_file, regions_names, coding=True, excluded_regions=None):

        Seq.__init__(self)
        self.coding = coding
        self.genes = GenesList()
        if excluded_regions: self.check_excluded_regions(excluded_regions)

        for region_name in regions_names:
            print(region_name)
            self.matrix = Matrix(self.coding)

            print('extracting annotations')
            regions = self.extract_regions_from_annotations(annotations_file, region_name)

            print('analysing sequence')
            self.analyse_maf_file(maf_file, regions)

            self.matrix.get_mutational_matrix().T.to_csv('results/'+region_name+"_matrix.csv", sep='\t')
            neutral_matrix = self.matrix.get_mutational_matrix()

            self.genes.calculate_exp_mutations(neutral_matrix)

        #TODO self.genes.export_mks('mks.csv')


    def extract_regions_from_annotations(self, annotation_file, region_name):
        ''' gets the annotation file lines
            returns a dictionary of exons positions, keys are the start position
            value is a list [start_pos, end_pos, chromosome,
            strand(- or + ), and the name of the gene]
            also returns the sorted keys of this dictionary (sorted start positions)
        '''
        fileSize = os.path.getsize(annotation_file)
        progress = 0
        exons = ExonsList()
        with open(annotation_file) as f:
            for _idx, line in enumerate(f):
                progress = progress + len(line)
                line = line.split('\t')
                if line[2] in region_name:
                    exon = Exon(line)
                    exons.add_exon(exon)

                Interpreter.print_progres((1.0*progress), fileSize)

        sys.stderr.write("\n")
        return exons

    def check_excluded_regions(self, excluded_regions):
        self.excluded_regions = ExcludedRegions()

        for region_name, filename in excluded_regions.items():
            lines = Importer.get_lines_from_file(filename)
            regions = self.extract_regions_from_annotations(lines, region_name)
            self.excluded_regions.add_regions(regions)

    def analyse_maf_file(self, maf_file, exons):
        ''' gets the maf_file lines, the dictionary of exons, and the exons keys 9start positions)
            list in order
            analysis the maf_file line by line, according to the exons from the annotations file
            analyse means generating counting the observed mutations per, and calculating the
            expected mutations per gene
            expected mutations are calculated by multiplying the number of occurence for each
            triblet and multiplying it by the chance of mutation from the neutral matrix
        '''
        contig = Contig("HCLCA", "hg38") #creat an initial empty contig
        #exon = exons[start_pos_list.pop(0)] #start with the first exon
        last_contig_pos = {}

        fileSize = os.path.getsize(maf_file)
        progress = 0

        with open (maf_file, "r") as f:
            for _idx, line in enumerate(f):
                progress = progress + len(line)

                if Contig.new_contig(line):
                    contig = Contig("HCLCA", "hg38")
                else:
                    contig.add_line(line)

                #if contig.got_both_seq():
                #    last_contig_pos = self.check_order(contig, last_contig_pos)

                if contig.got_both_seq() and exons.not_empty(contig.chrom):

                    exon = exons.get_first_exon(contig.chrom)

                    if contig.is_before_exon(exon):
                        contig = Contig("HCLCA", "hg38")
                        continue #go for the next contig

                    while contig.is_after_exon(exon):
                        try:
                            exon = exons.get_next_exon(contig.chrom)
                        except NoMoreExons:
                            break

                    while contig.has_exon(exon):
                        self.analyse_seq(contig, exon)
                        if contig.end_before_exon(exon):
                            contig = Contig("HCLCA", "hg38")
                            break

                        try:
                            exon = exons.get_next_exon(contig.chrom)
                        except NoMoreExons:
                            break

                Interpreter.print_progres((1.0*progress), fileSize)

    def analyse_seq(self, contig, exon):

        pointer = self.get_contig_pointer_by_exon(contig, exon)

        stop_point = min((exon.end_position - contig.strt_pos - 4),
                         (contig.length - 4))

        while pointer < stop_point:
            ref_codon = contig.ref_seq[pointer:pointer+3]
            pentamer = contig.ref_seq[pointer-1:pointer+4]
            alg_codon = contig.alg_seq[pointer:pointer+3]

            if "-" in ref_codon: #has insertion, unrechable now as we remove gaps first
                #gap_index = pointer + ref_codon.index("-")
                #contig.ref_seq = contig.ref_seq[:gap_index] + contig.ref_seq[gap_index+1:]
                #contig.alg_seq = contig.alg_seq[:gap_index] + contig.alg_seq[gap_index+1:]
                #stop_point -= 1

                pointer += 3
                continue

            else:
                base_ref_index = pointer + contig.strt_pos
                self.add_codon(ref_codon, alg_codon, pentamer, exon.negative_strand, exon.gene_name)
                pointer += 3

    def add_codon(self, ref_codon, alg_codon, pentamer, on_neg_strand, gene_name):
        if self.check_region_not_available(ref_codon, alg_codon, pentamer): return

        if on_neg_strand: #negative strand
            ref_codon = self.get_reverse_complement(ref_codon)
            alg_codon = self.get_reverse_complement(alg_codon)
            pentamer = self.get_reverse_complement(pentamer)

        self.matrix.add_tri_occ_coding(pentamer) if self.coding else 1+1# else self.matrix.add_tri_occ_noncoding(pentamer)
        self.genes.add_pentamer(gene_name, pentamer, alg_codon) if self.coding else self.genes.add_intron_occ(gene_name, pentamer)

        if alg_codon != ref_codon and self.coding: #mutation
            self.matrix.add_mutation(ref_codon, alg_codon, pentamer)

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

if __name__ == "__main__":

    ref_file = "../Data/refs/5LCAs_noSNVs_removed.maf"
    annotations = "hg38Regionsannotations.gtf"
    annotations = 'small_annotations.gtf'
    ref_file = 'trial.maf'

    Interpreter(ref_file, annotations, ["first_coding_exon",
                                        "internal_exon",
                                        "last_coding_exon"], True)
