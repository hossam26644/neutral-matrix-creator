
class Contig(object):

    def __init__(self, ref_spc, alg_spc):
        self.ref_spc = ref_spc
        self.alg_spc = alg_spc

        self.got_ref = False
        self.got_alg = False

        self.ref_seq = None
        self.alg_seq = None
        self.chrom = None
        self.strt_pos = None
        self.end_pos = None
        self.length = 0
        self.strand = None

    def got_both_seq(self):
        idx = 0
        if self.got_ref and self.got_alg:

            while idx < len(self.ref_seq):
                if self.ref_seq[idx] == "-":
                    self.ref_seq.pop(idx)
                    self.alg_seq.pop(idx)
                else:
                    idx += 1

            self.ref_seq = ("".join(self.ref_seq)).upper()
            self.alg_seq = ("".join(self.alg_seq)).upper()
            return True

        return False
        #return self.got_ref and self.got_alg

    def add_line(self, line):
        line = line.split()
        if line[0] == "s":
            species, chrom = line[1].split(".")

            if species == self.ref_spc:
                self.got_ref = True
                self.ref_seq = list(line[6])
                self.chrom = chrom
                self.strt_pos = int(line[2])+1
                self.end_pos = self.strt_pos + int(line[3])
                self.length = self.end_pos - self.strt_pos
                self.strand = line[4]

            elif species == self.alg_spc:
                self.got_alg = True
                self.alg_seq = list(line[6])

    def is_before_exon(self, exon):
        return self.end_pos < exon.start_position

    def is_after_exon(self, exon):
        return self.strt_pos > exon.end_position

    def has_exon(self, exon):
        no_exon = self.is_before_exon(exon) or self.is_after_exon(exon)
        return not no_exon

    @staticmethod
    def new_contig(line):
        line = line.split()
        return line == []
