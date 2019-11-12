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
logger.addHandler(ch)
prog_string = ''

g_columns = ['alt', 'non_alt', 'P', 'OR',]

sex_chr_re = re.compile(r'^(chr)?([XY])')

class CovAnalyzer(object):
    """
        Infer allele numbers at given depth thresholds from gnomAD
        coverage files
    """

    def __init__(self, coverage_files=[], coverage_directory=None,
                 dp_cutoff=10, pops_file=None, cohort="Exomes",
                 gnomad_version='3.0', genders_file=None):
        """
            Find valid coverage files in coverage_directory and determine
            what coverage cutoff to use (bisect available coverage
            threshold columns left on dp_cutoff). Create tabix
        """
        if not coverage_files and not coverage_directory:
            raise ValueError("One of 'coverage_files' or 'coverage_directory'"+
                             " args must be given.")
        self.logger = logging.getLogger("CovAnalyzer")
        self.logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter(
               '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        self.gnomad_version = gnomad_version
        self.file_to_dp = dict()    #cov column to use for each coverage file
        self.file_to_tbx = dict() #tabix iters for searching by coordinate
        self._collect_cov_files(coverage_files, coverage_directory, dp_cutoff)
        for f in self.file_to_dp:
            self.file_to_tbx[f] = pysam.TabixFile(f, parser=pysam.asTuple())
        self.pop_counts = self._read_pops(pops_file)
        valid_types = ["Exomes", "Genomes", "Total"]
        if cohort not in valid_types:
            sys.exit("Invalid cohort '{}' for CovAnalyzer".format(cohort)+
                     "valid types are {}.".format(", ".join(valid_types)))
        self.cohort = cohort
        self.fraction_xy = self._get_xy_fraction(genders_file)

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
        # assume that coverage files will match the input VCF
        #contig = contig.lstrip("chr")
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
                                    "pops.{}.tsv".format(self.gnomad_version))
        with open(pop_file, 'rt') as infile:
            header = next(infile)
            cols = self._columns_from_header(header)
            for row in (line.split() for line in infile):
                pop = row[cols['Population']]
                for x in ["Exomes", "Genomes", "Total"]:
                    counts[pop][x] = int(row[cols[x]])
        return counts

    def _read_genders(self, gender_file):
        if gender_file is None:
            gender_file = os.path.join(os.path.dirname(__file__),
                                       "data",
                                       "genders.2.1.tsv")
        return self._read_pops(gender_file)

    def _get_xy_fraction(self, gender_file):
        gender_counts = self._read_genders(gender_file)
        xx = gender_counts['female'][self.cohort]
        xy = gender_counts['male'][self.cohort]
        return float(xy)/(xx + xy)

    def _get_next_dp(self, dp, thresholds):
        thresholds.sort()
        try:
            return thresholds[bisect.bisect_left(thresholds, dp)]
        except IndexError:
            return thresholds[-1]

    def _get_dp_threshold(self, cols, dp):
        thresholds = [int(x.lstrip("over_")) for x in cols if
                      x.lstrip("over_").isdigit()]
        if not thresholds:
            return None
        if dp in thresholds:
            return dp
        return self._get_next_dp(dp, thresholds)

    def _check_header(self, header, f, dp):
        expected_cols = set(['chrom', "pos", "mean", "median"])
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
        self.file_to_dp[f] = cols["over_" + str(t)]
        return True

    def _cov_file_ok(self, f, dp):
        if os.path.isfile(f):
            o_func = gzip.open if f.endswith(('.gz', '.bgz')) else open
            try:
                with o_func(f, 'rt') as infile:
                    for line in infile:
                        if line.startswith('##'):
                            next
                        return self._check_header(line, f, dp)
            except ValueError:
                self.logger.warn("Error parsing {} - skipping".format(f))
                return False

    def _collect_cov_files(self, cov_files, cov_dir, dp):
        '''
            Set self.file_to_dp dictionary - keys are all valid coverage
            files, values are the columns of the closest depth threshold
            to that given to __init__. Setting of self.file_to_dp
            dictionary is handled within _cov_file_ok method.
        '''
        if cov_dir:
            cov_files.extend(os.path.join(cov_dir, f) for f in
                             os.listdir(cov_dir))

        valid_cov_files = [f for f in cov_files if self._cov_file_ok(f, dp)]
        if not valid_cov_files:
            sys.exit("No valid coverage files found - exiting")


def get_options():
    parser = argparse.ArgumentParser(description='''
            Do a crude association test against gnomAD VCF file(s).''')
    comgrp = parser.add_mutually_exclusive_group()
    parser.add_argument("-i", "--vcf_input", "--vcf", required=True,
                        help='Input VCF file')
    parser.add_argument("-e", "--exomes", help='gnomAD exomes VCF file')
    parser.add_argument("-g", "--genomes", help='gnomAD genomes VCF file')
    parser.add_argument("-t", "--table_output", help='''Table output of
                        variants, counts and p-values. Default=STDOUT''')
    parser.add_argument("-v", "--vcf_output", help='''VCF output. If provided
                        all variants will be written to this file with counts
                        and p-values annotated. If the filename ends with '.gz'
                        it will automatically compressed with bgzip.''')
    parser.add_argument("--pop_counts", metavar='FILE', help='''Tab delimited
                        file giving the number of exome and genome samples per
                        population. This is used with exome and genome coverage
                        files to infer the number of samples with likely
                        reference genotypes at sites without a variant called
                        in gnomAD VCF files. See the ".tsv" files in the data/
                        directory for examples.''')
    parser.add_argument("--gnomad_version", metavar='VERSION', default="3.0",
                        help='''gnomAD data version. This setting simply
                        chooses the default file used for sample counts per
                        population if --pop_counts is not set and exome/genome
                        coverage data is being used.''')
    parser.add_argument("--exome_coverage_files", nargs='+', help='''One or
                        more coverage summary files for Exome cohorts. Coverage
                        files must be sorted, bgzip compressed and tabix
                        indexed. These will be used to infer the number of
                        alleles at a given site above the coverage threshold if
                        a variant is not present in the gnomAD VCF file. If No
                        coverage files are provided p-values will only be
                        calculated for variants present in the gnomAD VCFs.''')
    parser.add_argument("--genome_coverage_files", nargs='+', help='''One or
                        more coverage summary files for whole genome
                        cohorts.''')
    parser.add_argument("--exome_coverage_dir", help='''Coverage directory as
                        downloaded from gnomAD containing one or more coverage
                        summary files for Exome cohorts. Coverage files must be
                        sorted, bgzip compressed and tabix indexed. These will
                        be used to infer the number of alleles at a given site
                        above the coverage threshold if a variant is not
                        present in the gnomAD VCF file. If No coverage files
                        are provided p-values will only be calculated for
                        variants present in the gnomAD VCF.''')
    parser.add_argument("--genome_coverage_dir", help='''Coverage directory as
                        above but for whole genome sequenced samples.''')
    parser.add_argument("-d", "--dp_cutoff", type=int, default=10,
                        help='''Calculate allele numbers using coverage files
                        at this depth threshold. Default=10''')
    parser.add_argument("-p", "--pops", nargs='+', help='''One or more gnomAD
                        populations to test against. Default to all 3-letter
                        population codes identified in VCF header excluding
                        ASJ/asj and OTH/oth.''')
    parser.add_argument("-s", "--samples", nargs='+', help='''One or more
                        samples to process. Defaults to all samples in input
                        file.''')
    parser.add_argument("--count_no_calls", action='store_true', help='''Count
                        chromosomes of samples even when they do not have a
                        called genotype for a variant or if they fail quality
                        filters. Default behaviour is to only count the number
                        of samples with a called genotype when calculating the
                        total number of chromosomes for a variant. This is
                        meant as a workaround where samples in your VCF were
                        not joint-called but instead merged together after
                        variant calling.''')
    parser.add_argument("-b", "--bed", help='''Restrict analysis to regions in
                        this BED format file.''')
    parser.add_argument("--p_value", type=float, default=0.05,
                        help='''Only output variants with a p-value equal to or
                        lower than this value. Default=0.05''')
    comgrp.add_argument("--require_all_p_values", action='store_true',
                        help='''Require the p-value to be under threshold for
                        all populations rather than just one.''')
    parser.add_argument("--test_combined_pops", action='store_true',
                        help='''Check --p-value threshold against combined
                        counts for all populations (see --pops) in addition to
                        individual populations.''')
    comgrp.add_argument("--test_combined_pops_only", action='store_true',
                        help='''Check --p-value threshold against combined
                        counts for all populations (see --pops) only.''')
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
    pop_ac_re = re.compile(r'''^AC_([A-Za-z]{3})$''')
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
    for p in pops + ['all']:
        for f in g_columns:
            ftype = 'Float' if f in ['P', 'OR'] else 'Integer'
            field_name = "gassoc_" + p + "_" + f
            if f == 'alt':
                desc = '"ALT allele counts for {} populations"'.format(p)
            elif f == 'non_alt':
                desc = '"non-ALT allele counts for {} populations"'.format(p)
            elif f == 'P':
                desc = ('"Fisher\'s P-value for {} populations vs '.format(p) +
                        'cohort"')
            elif f == 'OR':
                desc = ('"Odds ratio from Fisher\'s test for {} '.format(p) +
                       'populations vs cohort"')
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

def count_alleles(i, gts, gt_filter, indvs, max_one_per_sample, chromosome,
                  xx_samples, xy_samples, no_gender, count_no_calls=False):
    if count_no_calls:
        chrom_check = lambda x,y,z: True
    else:
        chrom_check = lambda x,y,z: (gt_filter.gt_is_ok(x, y, z) and
                                     gts['GT'][y] != (None,) *
                                     len(gts['GT'][y]))
    sex_match = sex_chr_re.match(chromosome)
    if max_one_per_sample:
        alts = sum(1 for s in indvs if gt_filter.gt_is_ok(gts, s, i)
                   and i in gts['GT'][s])
    elif sex_match:
        alts = sum(1 for s in no_gender if gt_filter.gt_is_ok(gts, s, i)
                   and i in gts['GT'][s])
        alts += sum(1 for s in xy_samples if gt_filter.gt_is_ok(gts, s, i)
                   and i in gts['GT'][s])
        if sex_match.group(2) == 'X':
            alts += sum((gts['GT'][s].count(i)) for s in xx_samples if
                       gt_filter.gt_is_ok(gts, s, i) and i in gts['GT'][s])
    else:
        alts = sum((gts['GT'][s].count(i)) for s in indvs if
                   gt_filter.gt_is_ok(gts, s, i) and i in gts['GT'][s])
    if sex_match:
        #TODO - rethink this fudge for samples without gender
        chroms = sum(2 for s in no_gender if chrom_check(gts, s, i))
        chroms += sum(1 for s in xy_samples if chrom_check(gts, s, i))
        if sex_match.group(2) == 'X':
            chroms += sum(2 for s in xx_samples if chrom_check(gts, s, i))
    else:
        chroms = sum(2 for s in indvs if chrom_check(gts, s, i))
    refs = chroms - alts
    return refs, alts

def process_variant(record, gnomad_filters, p_value, pops, gts, gt_filter,
                    table_out, vcf_out, cov_analyzers=dict(), pop_ids=None,
                    max_one_per_sample=False, families=None,
                    test_combined_pops=False, test_combined_pops_only=False,
                    require_all_p_values=False, xx_samples=None,
                    xy_samples=None, no_gender=None, count_no_calls=False):
    cohort_covered = dict()
    cov_ok = False
    alt_ref_counts = []
    if not pop_ids:
        pop_ids = sorted(set(p for x in pops.values() for p in x))
    for i in range(1, len(record.ALLELES)):
        if record.ALLELES[i] == '*':
            alt_ref_counts.append([-1, -1])
            continue
        if families is not None:
            indvs = _one_individual_per_fam(gts, gt_filter, i, families)
        else:
            indvs = [x for x in gts['GT']]
        refs, alts = count_alleles(i, gts, gt_filter, indvs, max_one_per_sample,
                                   record.CHROM, xx_samples, xy_samples,
                                   no_gender, count_no_calls=count_no_calls)
        alt_ref_counts.append([alts, refs])
    g_cohort_to_counts = defaultdict(list) #dict of pops -> [AC,AN]
    for c in ['Exomes', 'Genomes']:
        if c not in gnomad_filters:
            continue
        covered = dict()
        if c in cov_analyzers:
            covered = cov_analyzers[c].samples_at_site(record.CHROM,record.POS)
            if sum(covered.values()) == 0:
                continue  #no population covered sufficiently - skip cohort
        overlapping = gnomad_filters[c].get_overlapping_records(record)
        for i in range(len(record.DECOMPOSED_ALLELES)):
            if record.DECOMPOSED_ALLELES[i].ALT == '*':
                for p in pops[c]:
                    g_cohort_to_counts[p].append([-1, -1])
                continue
            filt,keep,matched,annot = gnomad_filters[c]._compare_var_values(
                record.DECOMPOSED_ALLELES[i], overlapping)
            for p in pops[c]:
                if len(g_cohort_to_counts[p]) <= i:
                    g_cohort_to_counts[p].append([0, 0])
                ac,an = 0,0
                if not annot: #no matching variant in gnomad VCF
                    if p in covered:
                        xy_match = sex_chr_re.match(record.CHROM)
                        if xy_match:
                            an = int(cov_analyzers[c].fraction_xy * covered[p])
                            if xy_match.group(2) == 'X':
                                xx = covered[p] - an
                                an += xx * 2
                        else:
                            an = 2 * covered[p]
                else:
                    ac = int(annot['AC_' + p])
                    an = int(annot['AN_' + p])
                g_cohort_to_counts[p][i][0] += ac
                g_cohort_to_counts[p][i][1] += an
    if not g_cohort_to_counts and not vcf_out: #no gnomAD data collected
        return
    #got counts for each cohort type for each allele - compute p-values
    p_check = all if require_all_p_values else any
    per_allele_results = []
    for i in range(len(record.DECOMPOSED_ALLELES)):
        allele = record.DECOMPOSED_ALLELES[i]
        if allele.ALT == '*':
            per_allele_results.append(['.'] * (4 * len(pop_ids) + 6))
            continue
        alts,refs = alt_ref_counts[i]
        results = [allele.CHROM, allele.POS, record.ID,  allele.REF,
                   allele.ALT, alts, refs]
        all_pvals = []
        total_ac, total_an = 0,0
        for p in pop_ids:
            if p in g_cohort_to_counts:
                ac, an = g_cohort_to_counts[p][i]
                odds, pval = stats.fisher_exact([(alts, refs), (ac, an-ac)],
                                                alternative='greater')
                all_pvals.append(pval)
                results.extend([ac, an - ac, pval, odds])
                total_ac += ac
                total_an += an
            else:
                results.extend(['.'] * 4)
        odds, pval = stats.fisher_exact([(alts, refs), (total_ac,
                                                        total_an-total_ac)],
                                        alternative='greater')
        results.extend([total_ac, total_an - total_ac, pval, odds])
        if test_combined_pops:
            all_pvals.append(pval)
        if vcf_out:
            per_allele_results.append(results[5:])
        if test_combined_pops_only:
            if pval <= p_value:
                table_out.write("\t".join((str(x) for x in results)) + "\n")
        else:
            if all_pvals and p_check(x <= p_value for x in all_pvals):
                table_out.write("\t".join((str(x) for x in results)) + "\n")
    if vcf_out:
        write_record(vcf_out, record, per_allele_results, pop_ids + ['all'])

def main(vcf_input, genomes=None, exomes=None, table_output=None,
         vcf_output=None, pops=None, samples=None, bed=None, p_value=0.05,
         dp_cutoff=10, exome_coverage_dir=None, genome_coverage_dir=None,
         exome_coverage_files=None, genome_coverage_files=None,
         gnomad_version='3.0', pop_counts=None, gq=20, dp=5, max_dp=250,
         het_ab=0.25, hom_ab=0.95, ped=None, test_combined_pops=False,
         test_combined_pops_only=False, require_all_p_values=False,
         max_one_per_sample=False, count_no_calls=False,
         progress_interval=1000, log_progress=False):
    if genomes is None and exomes is None:
        sys.exit('''At least one of --genomes or --exomes arguments must be
                 provided.''')
    vcfreader = VcfReader(vcf_input)
    out_fh = sys.stdout if table_output is None else open(table_output, 'w')
    vcf_writer = None
    if samples is None:
        samples = vcfreader.header.samples
    else:
        missing = [x for x in samples if x not in vcfreader.header.samples]
        if missing:
            sys.exit("Missing {} user specified samples".format(len(missing)) +
                     " in VCF ({}).".format(", ".join(missing)))
    families = None
    xx_samples = []
    xy_samples = []
    if ped is not None:
        pedfile = PedFile(ped)
        samples = [aff for aff in pedfile.get_affected() if aff in samples]
        if not samples:
            sys.exit("No affected samples from PED in VCF\n")
        families = defaultdict(list)
        for fid, fam in pedfile.families.items():
            families[fid].extend(x for x in fam.get_affected() if x in samples)
        xy_samples = [x for x in pedfile.get_males() if x in samples]
        xx_samples = [x for x in pedfile.get_females() if x in samples]
    no_gender = [x for x in samples if x not in xx_samples and x not in
                 xy_samples]
    gt_filter = GtFilter(vcfreader, gq=gq, dp=dp, max_dp=max_dp, het_ab=het_ab,
                         hom_ab=hom_ab)
    gt_fields = gt_filter.fields
    if bed is None:
        varstream = vcfreader
    else:
        logger.info("Reading, sorting and merging intervals in " +
                    "{}".format(bed))
        varstream = VarByRegion(vcfreader, bed=bed)
    gnomad_filters = dict()
    cov_analyzers = dict()
    cohort_pops = dict()
    for gnomad, cohort, cov_dir, cov_files in [
        (exomes, 'Exomes', exome_coverage_dir, exome_coverage_files),
        (genomes, 'Genomes', genome_coverage_dir, genome_coverage_files)
    ]:
        if gnomad is None:
            continue
        avail_pops = get_gnomad_pops(gnomad)
        if not avail_pops:
            sys.exit("ERROR: No gnomAD populations found in {}".format(gnomad))
        if pops:
            if not avail_pops.issuperset(set(pops)):
                logger.warn("The following specified populations were not " +
                            "found in your gnomAD file ({}): ".format(gnomad) +
                            ", ".join(set(pops).difference(avail_pops)))
            these_pops = [p for p in pops if p in avail_pops]
            if not these_pops:
                sys.exit("ERROR: No specified populations available in " +
                         "gnomAD file '{}'".format(gnomad))
        else:
            pop_blacklist = ['raw', 'asj', 'ASJ', 'oth', 'OTH']
            these_pops = [x for x in avail_pops if x not in pop_blacklist]
        gnomad_filters[cohort] = GnomadFilter(vcf=gnomad,
                                              prefix="gnomad_assoc",
                                              pops=these_pops)
        if cov_dir is not None or cov_files is not None:
            cov_analyzers[cohort] = CovAnalyzer(coverage_files=cov_files,
                                                coverage_directory=cov_dir,
                                                dp_cutoff=dp_cutoff,
                                                gnomad_version=gnomad_version,
                                                pops_file=pop_counts,
                                                cohort=cohort)
        cohort_pops[cohort] = sorted(list(these_pops))
    all_pops = sorted(set(p for x in cohort_pops.values() for p in x))
    if vcf_output is not None:
        if vcf_output.endswith(('.gz', '.bgz')):
            vcf_writer = bgzf.BgzfWriter(vcf_output)
        else:
            vcf_writer = open(vcf_output, 'w')
        write_vcf_header(vcfreader, vcf_writer, all_pops)
    header = ["#chrom", "pos", "id", "ref", "alt", "cases_alt", "cases_ref"]
    header.extend([x + "_alt\t" + x + "_ref\t" + x + "_p\t"  + x + "_odds" for
                   x in all_pops + ['total']])
    out_fh.write("\t".join(header) + "\n")
    n = 0
    for record in varstream:
        gts = record.parsed_gts(fields=gt_fields, samples=samples)
        process_variant(record=record, gnomad_filters=gnomad_filters,
                        pops=cohort_pops, pop_ids=all_pops,
                        cov_analyzers=cov_analyzers, gts=gts, table_out=out_fh,
                        vcf_out=vcf_writer, p_value=p_value,
                        gt_filter=gt_filter,
                        families=families,
                        test_combined_pops=test_combined_pops,
                        test_combined_pops_only=test_combined_pops_only,
                        require_all_p_values=require_all_p_values,
                        max_one_per_sample=max_one_per_sample,
                        xx_samples=xx_samples, xy_samples=xy_samples,
                        no_gender=no_gender, count_no_calls=count_no_calls)
        n += 1
        if n % progress_interval == 0:
            update_progress(n, record, log_progress)
    if out_fh is not sys.stdout:
        out_fh.close()
    if vcf_writer is not None:
        vcf_writer.close()

if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))
