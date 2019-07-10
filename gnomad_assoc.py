#!/usr/bin/env python3
import sys
import argparse
import re
import os
import logging
import bisect
import gzip
import pysam
import scipy.stats as stats
from collections import defaultdict
from parse_vcf import VcfReader
from vase.var_by_region import VarByRegion
from vase.gnomad_filter import GnomadFilter
from vase.sample_filter import GtFilter
from vase.ped_file import PedFile
from Bio import bgzf

logger = logging.getLogger("gnomAD Assoc")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter(
       '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
prog_string = ''
g_columns = ['alt', 'non_alt', 'P', 'OR',]

class CovAnalyzer(object):
    """
        Infer allele numbers at given depth thresholds from gnomAD
        coverage files
    """

    def __init__(self, coverage_directory, dp_cutoff=10, pops_file=None,
                 cohort="Exomes"):
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
        self.file_to_tbx = dict() #tabix iters for searching by coordinate
        self._collect_cov_files(coverage_directory, dp_cutoff)
        for f in self.file_to_dp:
            self.file_to_tbx[f] = pysam.TabixFile(f, parser=pysam.asTuple())
        self.pop_counts = self._read_pops(pops_file)
        valid_types = ["Exomes", "Genomes", "Total"]
        if cohort not in valid_types:
            sys.exit("Invalid cohort '{}' for CovAnalyzer".format(cohort)+
                     "valid types are {}.".format(", ".join(valid_types)))
        self.cohort = cohort

    def samples_at_site(self, contig, pos):
        frac = self.search_coordinates(contig, pos)
        if frac is None:
            return {None: 0}
        return dict((k, int(frac * v[self.cohort])) for (k, v) in
                    self.pop_counts.items())

    def search_coordinates(self, contig, pos):
        """
            Return fraction of samples above threshold at this position.
        """
        #TODO check whether coverage files have chr prefix - currently we
        # assume that coverage files will NOT have chr prefix
        contig = contig.lstrip("chr")
        for f, tbx in self.file_to_tbx.items():
            try:
                for row in tbx.fetch(contig, pos -1 , pos):
                    #first (and only?) hit should be the correct one
                    return float(row[self.file_to_dp[f]])
            except ValueError:
                pass
        return None

    def _columns_from_header(self, header):
        return dict((c, n) for (n, c) in enumerate(header.split()))

    def _read_pops(self, pop_file):
        counts = defaultdict(dict)
        if pop_file is None:
            pop_file = os.path.join(os.path.dirname(__file__),
                                    "data",
                                    "pops.tsv")
        with open(pop_file, 'rt') as infile:
            header = next(infile)
            cols = self._columns_from_header(header)
            for row in (line.split() for line in infile):
                pop = row[cols['Population']]
                for x in ["Exomes", "Genomes", "Total"]:
                    counts[pop][x] = int(row[cols[x]])
        return counts

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
        return self._get_next_dp(dp, thresholds)

    def _check_header(self, header, f, dp):
        expected_cols = set(['#chrom', "pos", "mean", "median"])
        cols = self._columns_from_header(header)
        if not expected_cols.issubset(set(cols)):
            self.logger.warn("Did not find expected column headers for file " +
                        "{} - skipping".format(f))
            return False
        t = self._get_dp_threshold(cols, dp)
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
            try:
                with o_func(f, 'rt') as infile:
                    for line in infile:
                        if line.startswith('##'):
                            next
                        return self._check_header(line, f, dp)
            except ValueError:
                self.logger.warn("Error parsing {} - skipping".format(f))
                return False

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
    parser.add_argument("-i", "--vcf_input", "--vcf", required=True,
                        help='Input VCF file')
    parser.add_argument("-g", "--gnomad", required=True, help='gnomAD VCF file'
                        )
    parser.add_argument("-t", "--table_output", help='''Table output of
                        variants, counts and p-values. Default=STDOUT''')
    parser.add_argument("-v", "--vcf_output", help='''VCF output. If provided
                        all variants will be written to this file with counts
                        and p-values annotated. If the filename ends with '.gz'
                        it will automatically compressed with bgzip.''')
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
    parser.add_argument("--cohort", action='store_true', default="Exomes",
                        help='''gnomAD cohort to use for inferring sample
                        number at sites without a variant. Population counts
                        will be calculated from coverage data provided to the
                        --coverage_dir argument. Valid values are "Exomes",
                        "Genomes" or "Total". Default=Exomes.''')
    parser.add_argument("-s", "--samples", nargs='+', help='''One or more
                        samples to process. Defaults to all samples in input
                        file.''')
    parser.add_argument("-b", "--bed", help='''Restrict analysis to regions in
                        this BED format file.''')
    parser.add_argument("--p_value", type=float, default=0.05,
                        help='''Only output variants with a p-value equal to or
                        lower than this value. Default=0.05''')
    parser.add_argument("--require_all_p_values", action='store_true',
                        help='''Require the p-value to be under threshold for
                        all populations rather than just one.''')
    parser.add_argument("--max_one_per_sample", action='store_true',
                        help='''Only count one allele per sample irrespective
                        of whether they are homozygous or heterozygous.''')
    parser.add_argument("--ped", help='''PED file detailing family membership
                        and affected status of individuals. Only individuals
                        with an affected status (2) will be counted and only a
                        single member will be counted per family. For variants
                        where any affected family member is confidently called
                        as homozygous reference a count of 0 will be given even
                        if other members carry the variant allele.''')
    parser.add_argument('--log_progress', action='store_true',
                        help='''Use logging output for progress rather than
                        wiping progress line after each update.''')
    parser.add_argument('--progress_interval', type=int, default=1000, metavar='N',
                        help='''Report progress information every N variants.
                        Default=1000.''')
    return parser

def get_gnomad_pops(vcf):
    vreader = VcfReader(vcf)
    pop_ac_re = re.compile(r'''^AC_([A-Z]+)$''')
    pops = []
    for f in vreader.metadata['INFO']:
        match = pop_ac_re.match(f)
        if match:
            p = match.group(1)
            an = 'AN_' + p
            ac_num = vreader.metadata['INFO'][f][-1]['Number']
            ac_typ = vreader.metadata['INFO'][f][-1]['Type']
            an_num = vreader.metadata['INFO'][an][-1]['Number']
            an_typ = vreader.metadata['INFO'][an][-1]['Type']
            if (ac_num == 'A' and ac_typ == 'Integer' and
                an_num == '1' and an_typ == 'Integer'):
                pops.append(p)
    if not pops:
        raise RuntimeError("No gnomAD populations found for VCF input!")
    return set(pops)

def _one_individual_per_fam(gts, gt_filter, allele, families):
    '''
        Return a single affected individual from each family on the
        condition that all carry the ALT allele (or are not confident
        calls).

        Args:
            gts:    parsed_gts from VcfRecord

            gt_filter:
                    GtFilter from vase.sample_filter

            allele: index of ALT allele to test

            families:
                    dict of family ID to affected members
    '''
    indvs = []
    for f,members in families.items():
        if len(members) == 1:
            indvs.append(members[0])
        else:
            a_counts = [gts['GT'][s].count(allele) if
                        gt_filter.gt_is_ok(gts, s, allele) else -1 for s in
                        members]
            if all(a_counts) and any(x > 0 for x in a_counts):
                #all carry variant or are no-call, but at least one is not a
                # no-call. If homs and hets present add the first het
                i = a_counts.index(min((x for x in a_counts if x > 0)))
                indvs.append(members[i])
    return indvs

def write_record(fh, record, results, pops):
    '''
        Annotate VcfRecord with counts and Fisher's test results from
        comparison with gnomAD.

        Args:
            fh: output filehandle

            record:
                input VcfRecord

            results:
                list with results (one per-allele) from process_variant
                method. Each per-allele result should contain the ALT
                allele counts from input, REF allele counts from input,
                and for each population the ALT allele counts, REF allele
                counts, p-value and odds-ratio.

            pops:
                The name of each population in the same order as it
                appears in each result.
    '''
    inf = defaultdict(list)
    for i in range(len(record.ALLELES)-1):
        inf['gassoc_cohort_alt'].append(results[i][0])
        inf['gassoc_cohort_non_alt'].append(results[i][1])
        for j in range(len(pops)):
            s = j * 4 + 2
            for f, r in zip(["gassoc_" + pops[j] + "_" + x for x in g_columns],
                            results[i][s:s+4]):
                inf[f].append(r)
    record.add_info_fields(inf)
    fh.write(str(record) + "\n")

def write_vcf_header(vcf, fh, pops):
    vcf.header.add_header_field(name="gnomad_assoc",
                               string='"' + str.join(" ", sys.argv) + '"')
    inf = {'gassoc_cohort_alt':     {'Number': 'A', 'Type': 'Integer',
                                     'Description':
                                     '"ALT allele counts in cohort"'},
           'gassoc_cohort_non_alt': {'Number': 'A', 'Type': 'Integer',
                                     'Description':
                                     '"non-ALT allele counts in cohort"'},
           }
    for p in pops:
        for f in g_columns:
            ftype = 'Float' if f in ['P', 'OR'] else 'Integer'
            field_name = "gassoc_" + p + "_" + f
            if f == 'alt':
                desc = '"ALT allele counts for {} population"'.format(p)
            elif f == 'non_alt':
                desc = '"non-ALT allele counts for {} population"'.format(p)
            elif f == 'P':
                desc = ('"Fisher\'s P-value for {} population vs '.format(p) +
                        'cohort"')
            elif f == 'OR':
                desc = ('"Odds ratio from Fisher\'s test for {} '.format(p) +
                       'population vs cohort"')
            inf[field_name] = {'Number': 'A', 'Type': ftype,
                               'Description': desc }
    for f,d in inf.items():
        vcf.header.add_header_field(name=f, dictionary=d, field_type='INFO')
    fh.write(str(vcf.header))

def update_progress(n, record, log=False):
    n_prog_string = "{:,} variants processed, at {}:{}".format(n, record.CHROM,
                                                               record.POS)
    global prog_string
    if log:
        logger.info(n_prog_string)
    else:
        n_prog_string = '\r' + n_prog_string
        if len(prog_string) > len(n_prog_string):
            sys.stderr.write('\r' + ' ' * len(prog_string) )
        sys.stderr.write(prog_string)
    prog_string = n_prog_string

def process_variant(record, gnomad_filter, p_value, pops, gts, gt_filter,
                    table_out, vcf_out, cov_analyzer=None,
                    max_one_per_sample=False, families=None,
                    require_all_p_values=False):
    if cov_analyzer:
        covered = cov_analyzer.samples_at_site(record.CHROM, record.POS)
        if sum(covered.values()) == 0:
            return #no population covered sufficiently
    overlapping = gnomad_filter.get_overlapping_records(record)
    p_check = all if require_all_p_values else any
    i = 0
    per_allele_results = []
    for allele in record.DECOMPOSED_ALLELES:
        if allele.ALT == '*':
            if vcf_out:
                per_allele_results.append(['.'] * (4 * len(pops) + 2))
            continue
        i += 1
        filt,keep,matched,annot = gnomad_filter._compare_var_values(allele,
                                                                    overlapping
                                                                   )
        if families is not None:
            indvs = _one_individual_per_fam(gts, gt_filter, i, families)
        else:
            indvs = [x for x in gts['GT']]
        if max_one_per_sample:
            alts = sum(1 for s in indvs if gt_filter.gt_is_ok(gts, s, i)
                       and i in gts['GT'][s])
        else:
            alts= sum((gts['GT'][s].count(i)) for s in indvs if
                       gt_filter.gt_is_ok(gts, s, i) and i in gts['GT'][s])
        #TODO rethink for sex chromosomes
        chroms = sum(2 for s in gts['GT'] if gt_filter.gt_is_ok(gts, s, i))
        refs = chroms - alts
        results = [allele.CHROM, allele.POS, record.ID,  allele.REF,
                   allele.ALT, alts, refs]
        all_pvals = []
        if not annot: #no matching variant in gnomad VCF
            if cov_analyzer:
                for p in pops:
                    an = covered[p]
                    odds, pval = stats.fisher_exact([(alts, refs), (0, an)],
                                                    alternative='greater')
                    all_pvals.append(pval)
                    results.extend([0, an, pval, odds])
        else:
            total_ac = 0
            total_an = 0
            for p in pops:
                ac = int(annot['AC_' + p])
                an = int(annot['AN_' + p])
                odds, pval = stats.fisher_exact([(alts, refs), (ac, an-ac)],
                                                alternative='greater')
                all_pvals.append(pval)
                results.extend([ac, an - ac, pval, odds])
                total_ac += ac
                total_an += an
            odds, pval = stats.fisher_exact([(alts, refs), (total_ac,
                                                            total_an-total_ac)],
                                            alternative='greater')
            results.extend([total_ac, total_an - total_ac, pval, odds])
        if vcf_out:
            per_allele_results.append(results[5:])
        if all_pvals and p_check(x <= p_value for x in all_pvals):
            table_out.write("\t".join((str(x) for x in results)) + "\n")
    if vcf_out:
        write_record(vcf_out, record, per_allele_results, pops)

def main(vcf_input, gnomad, table_output=None, vcf_output=None, pops=None,
         samples=None, bed=None, p_value=0.05, dp_cutoff=10, coverage_dir=None,
         cohort="Exomes", gq=20, dp=5, max_dp=250, het_ab=0.25, hom_ab=0.95,
         ped=None, require_all_p_values=False, max_one_per_sample=False,
         progress_interval=1000, log_progress=False):
    vcfreader = VcfReader(vcf_input)
    out_fh = sys.stdout if table_output is None else open(table_output, 'w')
    vcf_writer = None
    if vcf_output is not None:
        if vcf_output.endswith(('.gz', '.bgz')):
            vcf_writer = bgzf.BgzfWriter(vcf_output)
        else:
            vcf_writer = open(vcf_output, 'w')
        write_vcf_header(vcfreader, vcf_writer, pops)
    if samples is None:
        samples = vcfreader.header.samples
    else:
        missing = [x for x in samples if x not in vcfreader.header.samples]
        if missing:
            sys.exit("Missing {} user specified samples".format(len(missing)) +
                     " in VCF ({}).".format(", ".join(missing)))
    families = None
    if ped is not None:
        pedfile = PedFile(ped)
        samples = [aff for aff in pedfile.get_affected() if aff in samples]
        if not samples:
            sys.exit("No affected samples from PED in VCF\n")
        families = defaultdict(list)
        for fid, fam in pedfile.families.items():
            families[fid].extend(x for x in fam.get_affected() if x in samples)
    gt_filter = GtFilter(vcfreader, gq=gq, dp=dp, max_dp=max_dp, het_ab=het_ab,
                         hom_ab=hom_ab)
    gt_fields = gt_filter.fields
    if bed is None:
        varstream = vcfreader
    else:
        logger.info("Reading, sorting and merging intervals in " +
                    "{}".format(bed))
        varstream = VarByRegion(vcfreader, bed=bed)
    avail_pops = get_gnomad_pops(gnomad)
    if not avail_pops:
        sys.exit("ERROR: No gnomAD populations found in {}".format(pops))
    if pops:
        if not avail_pops.issuperset(set(pops)):
            sys.exit("ERROR - the following specified populations were not" +
                     " found in your gnomAD file ({}): ".format(gnomad) +
                     ", ".join(x for x in set(pops).difference(avail_pops)))
    else:
        pops = avail_pops
    gnomad_filter = GnomadFilter(vcf=gnomad, prefix="gnomad_assoc", pops=pops)
    cov_analyzer = None
    if coverage_dir is not None:
        cov_analyzer = CovAnalyzer(coverage_dir, dp_cutoff=dp_cutoff,
                                   cohort=cohort)
    pops = sorted(list(pops))
    header = ["#chrom", "pos", "id", "ref", "alt", "cases_alt", "cases_ref"]
    header.extend([x + "_alt\t" + x + "_ref\t" + x + "_p\t"  + x + "_odds" for
                   x in pops])
    out_fh.write("\t".join(header) + "\n")
    n = 0
    for record in varstream:
        gts = record.parsed_gts(fields=gt_fields, samples=samples)
        process_variant(record=record, gnomad_filter=gnomad_filter, pops=pops,
                        table_out=out_fh, vcf_out=vcf_writer, p_value=p_value,
                        cov_analyzer=cov_analyzer, gts=gts,
                        gt_filter=gt_filter, families=families,
                        require_all_p_values=require_all_p_values,
                        max_one_per_sample=max_one_per_sample)
        n += 1
        if n % progress_interval == 0:
            update_progress(n, record, log_progress)
    if out_fh is not sys.stdout:
        out_fh.close()

if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))
