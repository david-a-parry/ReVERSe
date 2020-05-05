import logging
import bisect
import os
import gzip
import pysam
from collections import defaultdict


class CovAnalyzer(object):
    """
        Infer allele numbers at given depth thresholds from gnomAD
        coverage files
    """

    def __init__(self, coverage_files=[], coverage_directory=None,
                 dp_cutoff=10, pops_file=None, cohort="Exomes",
                 gnomad_version='2.1', genders_file=None):
        """
            Find valid coverage files in coverage_directory and determine
            what coverage cutoff to use (bisect available coverage
            threshold columns left on dp_cutoff). Create tabix
        """
        if not coverage_files and not coverage_directory:
            raise ValueError("One of 'coverage_files' or " +
                             "'coverage_directory' args must be given.")
        if coverage_files is None:
            coverage_files = []
        self.logger = logging.getLogger("CovAnalyzer")
        self.logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter(
               '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        self.gnomad_version = gnomad_version
        self.file_to_dp = dict()   # cov column to use for each coverage file
        self.file_to_tbx = dict()  # tabix iters for searching by coordinate
        self._collect_cov_files(coverage_files, coverage_directory, dp_cutoff)
        for f in self.file_to_dp:
            self.file_to_tbx[f] = pysam.TabixFile(f, parser=pysam.asTuple())
        self.pop_counts = self._read_pops(pops_file)
        valid_types = ["Exomes", "Genomes", "Total"]
        if cohort not in valid_types:
            raise ValueError("Invalid cohort '{}' for ".format(cohort) +
                             "CovAnalyzer valid types are {}.".format(
                                 ", ".join(valid_types)))
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
        # TODO check whether coverage files have chr prefix - currently we
        # assume that coverage files will match the input VCF
        # contig = contig.lstrip("chr")
        for f, tbx in self.file_to_tbx.items():
            try:
                for row in tbx.fetch(contig, pos - 1, pos):
                    # first (and only?) hit should be the correct one
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
            raise ValueError("No valid coverage files found - exiting")
