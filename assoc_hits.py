#!/usr/bin/env python3
import sys
import argparse
import re
import os
import logging
import gzip
from collections import defaultdict
from parse_vcf import VcfReader
from vase.ped_file import PedFile
from vase.vep_filter import VepFilter
from vase.family_filter import FamilyFilter, RecessiveFilter
from vase.vase_runner import VariantCache
from Bio import bgzf

logger = logging.getLogger("gnomAD Assoc")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter(
       '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
prog_string = ''
progress_interval = 1000

def main(args):
    log_progress = False
    variant_cache = VariantCache()
    vcfreader = VcfReader(args.vcf_input)
    assoc_fields = get_assoc_fields(vcfreader)
    freq_fields = []
    if args.freq:
        freq_fields = get_vase_freq_annots(vcfreader)
    pedfile = PedFile(args.ped)
    family_filter = FamilyFilter(ped=pedfile, vcf=vcfreader,
                                 logging_level=logger.level)
    recessive_filter = RecessiveFilter(family_filter,
                                       gq=args.gq,
                                       dp=args.dp,
                                       max_dp=args.max_dp,
                                       het_ab=args.het_ab)
    #VepFilter for non-association hits (i.e. functional second hits)
    csq_filter = VepFilter(vcf=vcfreader,
                           csq=args.csq,
                           impact=args.impact,
                           canonical=args.canonical,
                           biotypes=args.biotypes,
                           in_silico=args.missense_filters,
                           filter_unpredicted=args.filter_unpredicted,
                           keep_any_damaging=args.keep_if_any_damaging,
                           loftee=args.loftee,
                           splice_in_silico=args.splice_filters,
                           splice_filter_unpredicted=args.splice_filter_unpredicted,
                           splice_keep_any_damaging=args.splice_keep_if_any_damaging,
                           retain_labels=args.retain_labels,
                           filter_flagged_features=args.flagged_features,
                           freq=args.freq,
                           afs=args.vep_af,
                           blacklist=args.feature_blacklist,
                           pathogenic=args.pathogenic,
                           no_conflicted=args.no_conflicted,
                           logging_level=logger.level)
    #VepFilter for association hits - check biotype, frequency etc. but not csq
    no_csq_filter = VepFilter(vcf=vcfreader,
                              csq=['all'],
                              canonical=args.canonical,
                              biotypes=args.biotypes,
                              filter_flagged_features=args.flagged_features,
                              freq=args.freq,
                              afs=args.vep_af,
                              blacklist=args.feature_blacklist,
                              logging_level=logger.level)
    if args.output is None:
        vcf_writer = sys.stdout
    elif args.output.endswith(('.gz', '.bgz')):
        vcf_writer = bgzf.BgzfWriter(args.output)
    else:
        vcf_writer = open(args.output, 'w')
    write_vcf_header(vcfreader, vcf_writer)
    n = 0
    for record in vcfreader:
        pr = process_record(record, args.p_value, args.min_alleles,
                            assoc_fields, csq_filter, no_csq_filter,
                            recessive_filter, args.freq, freq_fields)
        if pr:
            variant_cache.add_record(record)
        else:
            variant_cache.check_record(record)
        if variant_cache.output_ready:
            process_cache(variant_cache, recessive_filter, args.p_value,
                          args.min_alleles, assoc_fields, vcf_writer)
        n += 1
        if n % progress_interval == 0:
            update_progress(n, record, log_progress)
    process_cache(variant_cache, recessive_filter, args.p_value,
                  args.min_alleles, assoc_fields, vcf_writer, final=True)
    if vcf_writer is not sys.stdout:
        vcf_writer.close()

def process_cache(cache, recessive_filter, pval, min_alleles, assoc_fields,
                  outfh, final=False):
    #first check if any set of alleles segregates as recessive
    var_id_to_seg = recessive_filter.process_potential_recessives(final=final)
    #OrderedDict([('chr21:46357169-G/A', [<vase.family_filter.PotentialSegregant object at 0x7f2780d8f3f0>]), 
    if final:
        variant_cache.add_cache_to_output_ready()
    

def process_record(record, pval, min_alleles, assoc_fields, csq_filter,
                   hit_vep_filter, recessive_filter, freq, freq_fields):
    remove_alleles  = [False] * (len(record.ALLELES) -1)
    if freq_fields and freq:
        r = filter_on_existing_freq(record, freq, freq_fields)
        set_true_if_true(remove_alleles, r)
    is_hit = is_assoc_hit(record, pval, min_alleles)
    if is_hit:
        r_alts, remove_csq = hit_vep_filter.filter(record)
    else:
        r_alts, remove_csq = csq_filter.filter(record)
    set_true_if_true(remove_alleles, r_alts)
    if all(remove_alleles) or all(remove_csq):
        return
    return recessive_filter.process_record(record, remove_alleles, remove_csq)

def is_assoc_hit(record, pval, min_alleles):
    pass #TODO

def set_true_if_true(a, b):
    for i in range(len(a)):
        if b[i]:
            a[i] = True

def filter_on_existing_freq(record, freq, freq_fields):
    remove  = [False] * (len(record.ALLELES) -1)
    parsed = record.parsed_info_fields(fields=freq_fields)
    for annot in parsed:
        if parsed[annot] is None:
            continue
        for i in range(len(remove)):
            if parsed[annot][i] is not None:
                if parsed[annot][i] >= freq:
                    remove[i] = True
    return remove

def get_vase_freq_annots(vcf):
    frqs = list()
    for annot in vcf.metadata['INFO']:
        match = re.search('^VASE_dbSNP|gnomAD(_\d+)?_(CAF|AF)(_\w+)?', annot)
        if match:
            if (vcf.metadata['INFO'][annot][-1]['Number'] == 'A' and
                vcf.metadata['INFO'][annot][-1]['Type'] == 'Float'):
                logger.info("Found previous allele frequency annotation " +
                            "'{}'".format(annot))
                frqs.append(annot)
    return frqs

def get_assoc_fields(vcf):
    for f in ['gassoc_cohort_alt', 'gassoc_cohort_non_alt']:
        if f not in vcf.metadata['INFO']:
            sys.exit("Missing required annotation field '{}' in ".format(f) +
                     "VCF header. The gnomad_assoc.py script must be used to" +
                     " produce your input VCF before running this program.\n")
    annots = set()
    for info in vcf.metadata['INFO']:
        match = re.search('^gassoc_\w+_P', info)
        if match:
            annots.add(info)
            logger.debug("Identified association INFO field '{}'".format(info))
    if not annots:
        sys.exit("No P-value annotations found in VCF header. Exiting.\n")
    return annots

def update_progress(n, record, log=False):
    n_prog_string = "{:,} variants processed, at {}:{}".format(n,
                                                               record.CHROM,
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


def write_vcf_header(vcf, fh):
    vcf.header.add_header_field(name="assoc_hits",
                               string='"' + str.join(" ", sys.argv) + '"')
    inf = dict()
    #TODO - get new INFO fields
    for f,d in inf.items():
        vcf.header.add_header_field(name=f, dictionary=d, field_type='INFO')
    fh.write(str(vcf.header))

def get_options():
    parser = argparse.ArgumentParser(description='''
            Find genes with an association hit on one allele and an allele
            matching consequence requirements on the other.''')
    parser.add_argument("-i", "--vcf_input", "--vcf", required=True,
                        help='Input VCF file')
    parser.add_argument("-p", "--ped", required=True,
                        help='Ped file detailing family relationships')
    parser.add_argument("-o", "--output", help='Output VCF file')
    parser.add_argument("-v", "--p_value", type=float, default=1e-4,
                        help='''Minimum association test p-value.
                        Default=1e-4''')
    parser.add_argument("-a", "--min_alleles", type=int, default=2, help=
                        '''Minimum observed alleles from association test.
                        Default=2.''')
    parser.add_argument('--freq', type=float, default=0.0, help=
                        '''Allele frequency cutoff. Frequency information will
                        be read from existing VASE and VEP annotations.
                        Variant alleles with an allele frequency equal to or
                        greater than this value will be filtered.''')
    parser.add_argument('--gq', type=int, default=20, help=
                         '''Minimum genotype quality score threshold. Sample
                         genotype calls with a score lower than this threshold
                         will be treated as no-calls. Default = 20.''')
    parser.add_argument('-dp', '--dp', type=int, default=0, help=
                        '''Minimum genotype depth threshold. Sample genotype
                        calls with a read depth lower than this threshold
                        will be treated as no-calls. Default = 0.''')
    parser.add_argument('-max_dp', '--max_dp', type=int, default=0, help=
                        '''Maximum genotype depth threshold. Sample genotype
                        calls with a read depth higher than this threshold
                        will be treated as no-calls. Default = 0 (i.e. not
                        used).''')
    parser.add_argument('-het_ab', '--het_ab', type=float, default=0.,
                         metavar='AB', help='''Minimum genotype allele balance
                         for heterozygous genotypes. Heterozygous sample
                         genotype calls with a ratio of the alternate allele vs
                         total depth lower than this threshold will be treated
                         as no-calls. Default = 0.''')
    parser.add_argument('--csq', nargs='+', help=
                        '''One or more VEP consequence classes to retain.
                        Variants which do not result in one of these VEP
                        consequence classes will be filtered.
                        Default=['TFBS_ablation', 'TFBS_amplification',
                        'inframe_deletion', 'inframe_insertion',
                        'frameshift_variant', 'initiator_codon_variant',
                        'missense_variant', 'protein_altering_variant',
                        'regulatory_region_ablation',
                        'regulatory_region_amplification',
                        'splice_acceptor_variant', 'splice_donor_variant',
                        'start_lost', 'stop_gained', 'stop_lost',
                        'transcript_ablation', 'transcript_amplification']''')
    parser.add_argument('--impact', nargs='+', help='''One or more VEP 'IMPACT'
                        types to retain. Valid values are 'HIGH', 'MODERATE',
                        'LOW' and 'MODIFIER'. Any consequence classes specified
                        by the '--csq' argument will still be retained
                        irrespective of values specified here.''')
    parser.add_argument('--canonical', action='store_true', help=
                        '''When used in conjunction with --csq argument,
                        ignore consequences for non-canonical transcripts.''')
    parser.add_argument('--flagged_features', action='store_true', help=
                        '''Ignore consequences for flagged transcripts/features
                        (i.e. with a non-empty 'FLAGS' CSQ field).''')
    parser.add_argument('--biotypes',  nargs='+', default=[],
                        metavar='BIOTYPE', help='''Ignore consequences in
                        biotypes other than those specified here. Default =
                        ['3prime_overlapping_ncrna', 'antisense',
                        'CTCF_binding_site', 'enhancer', 'IG_C_gene',
                        'IG_D_gene', 'IG_J_gene', 'IG_V_gene', 'lincRNA',
                        'miRNA', 'misc_RNA', 'Mt_rRNA', 'Mt_tRNA',
                        'open_chromatin_region', 'polymorphic_pseudogene',
                        'processed_transcript', 'promoter',
                        'promoter_flanking_region', 'protein_coding', 'rRNA',
                        'sense_intronic', 'sense_overlapping', 'snoRNA',
                        'snRNA', 'TF_binding_site',
                        'translated_processed_pseudogene', 'TR_C_gene',
                        'TR_D_gene', 'TR_J_gene', 'TR_V_gene']''')
    parser.add_argument('--feature_blacklist', '--blacklist', help=
                        '''A file containing a list of Features (e.g. Ensembl
                        transcript IDs) to ignore. These must correspond
                        to the IDs in the 'Feature' field annotated by VEP.''')
    parser.add_argument('--loftee', default=False, action='store_true', help=
                        '''Retain LoF (stop_gained, frameshift_variant,
                        splice_acceptor_variant and splice_donor_variant)
                        classes only if the LoF annotation from loftee is
                        'HC'.''')
    parser.add_argument('-m', '--missense_filters', default=[], nargs='+',
                        help='''A list of in silico prediction programs to
                        use for filtering missense variants (must be used
                        in conjunction with --csq argument). The programs
                        provided here must have been annotated on the
                        input VCF file either directly by VEP or via the
                        dbNSFP VEP plugin. Recognised program names and
                        default 'damaging' values are provided in the
                        "data/vep_insilico_pred.tsv" file.

                         You may optionally specify score criteria for
                         filtering (e.g. FATHMM_pred=D) or you may just provide
                         the program names and the default 'damaging'
                         prediction values will be used.

                         By default, a missense consequence is filtered unless
                         each of the programs listed here have an appropriate
                         or missing prediction/score. This behaviour can be
                         changed using the --filter_unpredicted or
                         --keep_if_any_damaging flags.''')

    parser.add_argument('--filter_unpredicted', action='store_true',
                        default=False, help='''For use in conjunction with
                        --missense_filters. The default behaviour when using
                        --missense_filters is to ignore a program if there is
                        no prediction given (i.e. the score/pred is empty).
                        That is, if there are no predictions for any of the
                        programs annotating a missense consequence, it will not
                        be filtered, while if predictions are missing for only
                        some, filtering will proceed as normal with the other
                        programs. If this option is given, missense variants
                        will be filtered if any program does not have a
                        prediction/score.''')
    parser.add_argument('--keep_if_any_damaging', action='store_true',
                        default=False, help= '''For use in conjunction with
                        --missense_filters. If this option is provided, a
                        missense consequence is only filtered if ALL of the
                        programs provided to --missense_filters do not have an
                        appropriate prediction/score - that is, the missense
                        consequence will be retained if ANY of the given
                        programs has an appropriate value for the
                        prediction/score. This behaviour is overridden by
                        '--filter_unpredicted' when a prediction/score is
                        missing for any program.''')

    parser.add_argument('--splice_filters', nargs='+', help='''Similar to
                        --missense_filters except only splice consequences
                        (splice_donor_variant, splice_acceptor_variant and
                        splice_region_variant) are checked versus the given in
                        silico prediction programs. Currently only dbscSNV,
                        (rf_score and ada_score), MaxEntScan and SpliceDistance
                        (https://github.com/david-a-parry/SpliceDistance)
                        plugins are supported. For example '--splice_filters
                        ada_score' will filter splice region variants with a
                        dbscSNV ada_score cutoff below the default value (0.7).
                        Alternatively, '--splice_filters ada_score=0.9' would
                        filter on a higher threshold of 0.9 or above.''')

    parser.add_argument('--splice_filter_unpredicted', action='store_true',
                        default=False, help='''Same as --filter_unpredicted but
                        for --splice_filters only.''')

    parser.add_argument('--splice_keep_if_any_damaging', action='store_true',
                        default=False, help= '''Same as --keep_if_any_damaging
                        but for --splice_filters only.''')
    parser.add_argument('--retain_labels', metavar='Label=Value', nargs='+',
                        help='''Retain consequence annotations if there is a
                        matching annotation for the given label. For example,
                        to retain any consequence where there is a VEP
                        annotation for 'FOO' matching 'BAR' use
                        "--retain_labels FOO=BAR".''')
    parser.add_argument('--no_vep_freq', '-no_vep_freq', action='store_true',
                        help='''Use this option if you want to ignore VEP
                        annotated allele frequencies when using the --freq
                        option.''')
    parser.add_argument('--vep_af', '-vep_af', nargs='+', default=[], help=
                        '''One or more VEP allele frequency annotations to
                        use for frequency filtering. Default is to use all
                        standard VEP AF annotations.''')

    parser.add_argument('--pathogenic', action='store_true', help='''When used
                        in conjunction with --csq argument, retain variants
                        flagged as pathogenic by either 'CLIN_SIG' or
                        'clinvar_clnsig' VEP annotations even if the
                        consequence class is not included in those selected
                        using the --csq argument. Note that this only alters
                        filtering as specified by --csq and --missense_filters
                        options; frequency, canonical transcript,
                        flagged_features and biotype filtering will still occur
                        as normal.''')
    parser.add_argument('--no_conflicted', action='store_true', help='''When
                        used in conjunction with --pathogenic argument,
                        variants labelled as pathogenic will only be retained
                        if there are no conflicting 'benign' or 'likely benign'
                        assertions.''')
    return parser

if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(args)
