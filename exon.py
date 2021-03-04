import bisect


class Exon(object):
    """
    docstring
    """
    def __init__(self, line):


        self.gene_name = self._get_gene_name(line[8])
        self.chrom = line[0]

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

    def _get_gene_name(self, line8):
        """
        docstring
        """
        try:
            gene_name = line8.split()[1]
            gene_name = gene_name[1:-2] #remove semicolon and qutation marks
        except IndexError:
            gene_name = 'Not_a_gene'

        return gene_name

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
        self.strt_position_by_chrom2 = {}

    def add_exon(self, exon):

        if exon.chrom in self.exons_by_chrom:
            if exon.start_position not in self.strt_position_by_chrom[exon.chrom]:
                self.exons_by_chrom[exon.chrom][exon.start_position] = exon
                bisect.insort(self.strt_position_by_chrom[exon.chrom], exon.start_position)
                bisect.insort(self.strt_position_by_chrom2[exon.chrom], exon.start_position)
            else:
                return
        else:
            self.exons_by_chrom[exon.chrom] = {exon.start_position:exon}
            self.strt_position_by_chrom[exon.chrom] = [exon.start_position]
            self.strt_position_by_chrom2[exon.chrom] = [exon.start_position]

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
