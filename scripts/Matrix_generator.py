from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import sys

sample=sys.argv[1]
path=sys.argv[2]
matrices = matGen.SigProfilerMatrixGeneratorFunc(sample, "GRCh37", path, exome=False, bed_file=None, chrom_based=False, plot=True, tsb_stat=False)
