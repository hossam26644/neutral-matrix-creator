import bisect

class ExcludedRegions(object):
    def __init__(self):
        self.start_positions_by_chrom = {}
        self.end_positions_by_chrom = {}
        self.current_excluded_region = ('bla', range(0, 0))
        self.current_included_region = ('bla', range(0, 0))


    def add_regions(self, regions):
        if not self.start_positions_by_chrom:
            self.start_positions_by_chrom = regions.strt_position_by_chrom
            self.end_positions_by_chrom = regions.end_position_by_chrom
            return

        for chrom, start_positions in regions.strt_position_by_chrom.items():
            if not self.start_positions_by_chrom:
                self.start_positions_by_chrom[chrom] = start_positions
                self.end_positions_by_chrom[chrom] = regions.end_position_by_chrom[chrom]
            else:
                #TODO implement insorting the new regions
                pass

    def is_excluded(self, chrom, pos):
        if (self.current_excluded_region[0] == chrom) and (pos in self.current_excluded_region[1]):
            return True
        elif (self.current_included_region[0] == chrom) and (pos in self.current_included_region[1]):
            return False
        else:
            return self.chec_included_excluded(chrom, pos)

    def chec_included_excluded(self, chrom, pos):
        if not self.start_positions_by_chrom:
            return False

        try:
            idx = bisect.bisect_left(self.start_positions_by_chrom[chrom], pos) - 1
            if self.start_positions_by_chrom[chrom][idx] < pos < self.end_positions_by_chrom[chrom][idx]:

                self.current_excluded_region = (chrom, range(self.start_positions_by_chrom[chrom][idx],
                                                            self.end_positions_by_chrom[chrom][idx]))
                return True
            else:
                try:
                    self.current_included_region = (chrom,
                                                    range(pos, self.start_positions_by_chrom[chrom][idx+1]))
                except IndexError:
                    pass
            return False

        except KeyError:
            return False
