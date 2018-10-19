#!/usr/bin/env python3
import sys
import argparse
import re
import logging
import bisect
import gzip
import pysam
from vase.var_by_region import VarByRegion
from vase.gnomad_filter import GnomadFilter

logger = logging.getLogger("gnomAD Assoc")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter(
       '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)

class CovAnalyzer(object):

    """
        Infer allele numbers at given depth thresholds from gnomAD
        coverage files
    """

    def __init__(self, coverage_directory, dp_cutoff=10, pops_file=None):
        """
            Find valid coverage files in coverage_directory and determine
            what coverage cutoff to use (bisect available coverage
            threshold columns left on dp_cutoff). Create tabix
        """
        self.logger = logging.getLogger("CovAnalyzer")
        self.logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter(
               '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        self.file_to_dp = dict()    #cov column to use for each coverage file
        self.file_to_tabix = dict() #tabix iters for searching by coordinate
        self._collect_cov_files(coverage_directory, dp_cutoff)
        for f in self.file_to_dp:
            self.file_to_tbx[f] = pysam.TabixFile(f, parser=pysam.asTuple())
        #TODO Read population file and get number of samples per population

    def search_coordinates(contig, pos):
        """
            Return fraction of samples above threshold at this position.
        """
        for f, tbx in self.file_to_tbx.items():
            for row in tbx.fetch(contig, pos -1 , pos):
                #first (and only?) hit should be the correct one
                return float(row[self.file_to_dp[f]])
        return None

    def _columns_from_header(header):
        return dict((c, n) for (n, c) in enumerate(header.split()))

    def _get_next_dp(self, dp, thresholds):
        thresholds.sort()
        try:
            return thresholds[bisect.bisect_left(thresholds, dp)]
        except IndexError:
            return thresholds[-1]

    def _get_dp_threshold(self, cols, dp):
        thresholds = [int(x) for x in cols if x.isdigit()]
        if not thresholds:
            return None
        if dp in thresholds:
            return dp
        return get_next_dp(dp, thresholds)

    def _check_header(self, header, f, dp):
        expected_cols = set(['#chrom', "pos", "mean", "median"])
        cols = columns_from_header(header)
        if not expected_cols.issubset(set(cols)):
            self.logger.warn("Did not find expected column headers for file " +
                        "{} - skipping".format(f))
            return False
        t = get_dp_threshold(cols, dp)
        if t is None:
            self.logger.warn("Did not find any depth thresholds headers for " +
                        "file {} - skipping".format(f))
            return False
        if dp != t:
            self.logger.warn("Depth threshold column {} not ".format(dp) +
                             "found in {}".format(f) + " - will use a " +
                             "threshold of {}".format(t))
        self.file_to_dp[f] = cols[str(t)]
        return True

    def _cov_file_ok(self, f, dp):
        if os.path.isfile(f):
            o_func = gzip.open if f.endswith('.gz') else open
            with o_func(f, 'rt') as infile:
                for line in infile:
                    if line.startswith('##'):
                        next
                    return check_header(line, f, dp)

    def _collect_cov_files(self, cov_dir, dp):
        '''
            Set self.file_to_dp dictionary - keys are all valid coverage
            files, values are the columns of the closest depth threshold
            to that given to __init__.
        '''
        cov_files = [os.path.join(cov_dir, f) for f in os.listdir(cov_dir) if
                     self._cov_file_ok(os.path.join(cov_dir, f), dp)]
        if not cov_files:
            sys.exit("No valid coverage files found - exiting")



def get_options():
    parser = argparse.ArgumentParser(description='''
            Do a crude association test against a gnomAD VCF file.''')
    parser.add_argument("vcf", help='''Input VCF file''')
    parser.add_argument("gnomad", help='''gnomAD VCF file''')
    parser.add_argument("-c", "--coverage_dir", help='''Coverage directory as
                        downloaded from gnomAD containing one or more coverage
                        files. Coverage files must be sorted, bgzip compressed
                        and tabix indexed. These will be used to infer the
                        number of alleles at a given site above the coverage
                        threshold if not present in the gnomAD VCF file. If No
                        coverage files are provided p-values will only be
                        calculated for variants present in the gnomAD VCF. Note
                        that these coverage files should correspond to the same
                        data used in the gnomAD VCF file (e.g. exome or WGS).''')
    parser.add_argument("-d", "--dp_cutoff", type=int, default=10,
                        help='''Calculate allele numbers using coverage files
                        at this depth threshold. Default=10''')
    parser.add_argument("-p", "--pops", nargs='+', help='''One or more gnomAD
                        populations to test against.''')
    parser.add_argument("-s", "--samples", nargs='+', help='''One or more
                        samples to process. Defaults to all samples in input
                        file.''')
    parser.add_argument("-b", "--bed", help='''Restrict analysis to regions in
                        this BED format file.''')
    parser.add_argument("--p_value", type=float, default=0.05,
                        help='''Only output variants with a p-value equal to or
                        lower than this value. Default=0.05''')
    return parser

def get_gnomad_pops(vcf):
    pop_ac_re = re.compile(r'''^AC_([A-Z]+)$''')
    pops = []
    for f in vcf.metadata['INFO']:
        match = pop_ac_re.match(f)
        if match:
            p = match.group(1)
            an = 'AN_' + p
            ac_num = vcf.metadata['INFO'][f][-1]['Number']
            ac_typ = vcf.metadata['INFO'][f][-1]['Type']
            an_num = vcf.metadata['INFO'][an][-1]['Number']
            an_typ = vcf.metadata['INFO'][an][-1]['Type']
            if (ac_num == 'A' and ac_typ == 'Integer' and
                an_num == '1' and an_typ == 'Integer'):
                pops.append(p)
    if not pops:
        raise RuntimeError("No gnomAD populations found for VCF input!")
    return set(pops)

def process_variant(record, gnomad_Filter, p_value, cov_analyzer=None):
    pass

def main(vcf, gnomad, pops=None, samples=None, bed=None, p_value=0.05,
         dp_cutoff=10, coverage_dir=None):
    vcfreader = VcfReader(vcf)
    if bed is None:
        varstream = vcfreader
    else:
        logger.info("Reading, sorting and merging intervals in " +
                    "{}".format(bed))
        varstream = VarByRegion(vcfreader, bed=bed)
    gnomad_filter = GnomadFilter(vcf=gnomad, prefix=gnomad_assoc)
    avail_pops = get_gnomad_popts(gnomad_filter.vcf)
    if not avail_pops:
        sys.exit("ERROR: No gnomAD populations found in {}".format(pops))
    if pops:
        if not avail_pops.issuperset(set(pops)):
            sys.exit("ERROR - the following specified populations were not" +
                     " found in your gnomAD file ({}): ".format(gnomad) +
                     ", ".join(x for x in set(pops).difference(avail_pops)))
    else:
        pops = avail_pops
    cov_analyzer = None
    if coverage_dir is not None:
        cov_analyzer = CovAnalyzer(coverage_dir, dp_cutoff)
    pops = sorted(list(pops))
    header = ["#chrom", "pos", "id", "ref", "alt", "cases_ref", "cases_alt"]
    header.extend([x + "_ref\t" + x + "_alt" for x in pops])
    print("\t".join(header))
    for record in varstream:
        process_variant(record, gnomad_filter, p_value, cov_analyzer)

if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))
