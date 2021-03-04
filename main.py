from extractor import Extractor
from seq import Seq
from importer import Importer

class MatrixCreator(object):

    def __init__(self, regions):
        for region in regions:
            excluded_regions = {'CpGisland': 'CpGislands.gtf'}
            extractor = Extractor("MRCA_mult_full.maf", "hg38Regionsannotations.gtf")
            results = extractor.extract_mutational_matrix(region, False, excluded_regions)

            Seq.export_matrix(region, 'neutral_matrix', results[0])
            Seq.export_matrix(region, 'mutations', results[1])
            Seq.export_matrix(region, 'occurences', results[2])
            Importer.export_text("neutral_matrices/" + region +'_logs.txt', extractor.logs)


if __name__ == '__main__':
    MatrixCreator(["5UTR", "3UTR", "Intergenic", "intron"])
