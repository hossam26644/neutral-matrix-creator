import bisect
import re

class Exon(object):
    """
    docstring
    """
    def __init__(self, line):


        self.gene_name = self._get_gene_name(line[8])
        self.chrom = line[0]
        self.replication_dir = self.get_repli_dir(line[8])

        direction = line[6]
        frame = self.get_frame(line[7])

        if direction == "+":
            self.def_strt_end_postn_pstv(line[3], line[4], frame)
            self.negative_strand = False

        elif direction == "-":
            self.def_strt_end_postn_ngtv(line[3], line[4], frame)
            self.negative_strand = True

        #self.gene_name = self.gene_name + str(self.start_position)

    def get_frame(self, frame):
        if frame in ['0', '1', '2']:
            return frame
        return '0'

    def def_strt_end_postn_pstv(self, strt_pos, end_pos, frame):
        """
        docstring
        """

        self.start_position = int(strt_pos) + int(frame) #start by the first codon
        self.end_position = int(end_pos)

    def def_strt_end_postn_ngtv(self, strt_pos, end_pos, frame):
        """
        docstring
        """
        self.start_position = int(strt_pos)
        self.end_position = int(end_pos) - int(frame) - 2

    def _get_gene_name(self, line):

        attrs = line.split(";")[:-1]
        for attr in attrs:
            attr = attr.split()
            if attr[0] == "gene_name":
                return re.sub('\"', '', attr[1])

        return 'Not_a_gene'

    def get_repli_dir(self, line):

        attrs = line.split(";")[:-1]
        for attr in attrs:
            attr = attr.split()
            if attr[0:1] == "replication strand":
                return re.sub('\"', '', attr[-1])

        return 'no dir'

    def __eq__(self, other):
        """
        docstring
        """
        if self.chrom == other.chrom:
            if self.start_position == other.start_position:
                if self.end_position == other.end_position:
                    return True

        return False

class ExonsList(object):
    """
    docstring
    """
    def __init__(self):
        self.exons_by_chrom = {}
        self.strt_position_by_chrom = {}
        self.end_position_by_chrom = {}

    def add_exon(self, exon):

        if exon.chrom in self.exons_by_chrom:
            if exon.start_position not in self.strt_position_by_chrom[exon.chrom]:
                self.exons_by_chrom[exon.chrom][exon.start_position] = exon
                bisect.insort(self.strt_position_by_chrom[exon.chrom], exon.start_position)
                bisect.insort(self.end_position_by_chrom[exon.chrom], exon.end_position)
            else:
                return
        else:
            self.exons_by_chrom[exon.chrom] = {exon.start_position:exon}
            self.strt_position_by_chrom[exon.chrom] = [exon.start_position]
            self.end_position_by_chrom[exon.chrom] = [exon.end_position]

    def get_first_exon(self, chrom):
        frst_strt_pos = self.strt_position_by_chrom[chrom][0]
        return self.exons_by_chrom[chrom][frst_strt_pos]

    def get_next_exon(self, chrom):
        self.strt_position_by_chrom[chrom].pop(0)
        if self.not_empty(chrom):
            return self.get_first_exon(chrom)
        else:
            raise NoMoreExons

    def not_empty(self, chrom):
        """
        docstring
        """
        if chrom in self.strt_position_by_chrom and self.strt_position_by_chrom[chrom]:
            return True

        return False

class NoMoreExons(ValueError):
    pass
