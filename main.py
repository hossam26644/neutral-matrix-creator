from extractor import Extractor

extractor = Extractor("mrca_mult.maf", "hg38Regionsannotations.gtf")
results = extractor.extract_mutational_matrix("intron", False)

print(results)
