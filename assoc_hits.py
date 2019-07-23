#!/usr/bin/env python3
import sys
import argparse
import re
import os
import logging
import gzip
from collections import defaultdict, OrderedDict
from parse_vcf import VcfReader
from vase.ped_file import PedFile
from vase.vep_filter import VepFilter
from vase.family_filter import FamilyFilter, RecessiveFilter
from vase.family_filter import SegregatingVariant
from vase.vase_runner import VariantCache
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

class AssocSegregator(RecessiveFilter):
    '''
        Look for recessive combinations including at least one designated
        association 'hit'

        Arguments are same as for RecessiveFilter except:

                max_incidentals:
                        Maximum number of families with an individual
                        that may carries an allele of a biallelic variant
                        where affected individuals in that family do not
                        have a second qualifying variant. Useful for
                        filtering out implausibly common alleles that by
                        chance segregate in a subset of families.
                        Default=0 (i.e. not applied).

    '''

    def __init__(self, family_filter, gq=0, dp=0, max_dp=0, het_ab=0.,
                 hom_ab=0., min_control_dp=None, max_control_dp=None,
                 min_control_gq=None, control_het_ab=None, control_hom_ab=None,
                 con_ref_ab=None, sv_gq=0, sv_dp=0, sv_max_dp=0,
                 sv_het_ab=0., sv_hom_ab=0., sv_min_control_dp=None,
                 sv_max_control_dp=0, sv_min_control_gq=None,
                 sv_control_het_ab=None, sv_control_hom_ab=None,
                 sv_con_ref_ab=None, min_families=1, strict=False,
                 exclude_denovo=False, report_file=None, max_incidentals=0):
        super().__init__(family_filter, gq=gq, dp=dp, max_dp=max_dp,
                         het_ab=het_ab, hom_ab=hom_ab,
                         min_control_dp=min_control_dp,
                         max_control_dp=max_control_dp,
                         min_control_gq=min_control_gq,
                         control_het_ab=control_het_ab,
                         control_hom_ab=control_hom_ab, con_ref_ab=con_ref_ab,
                         sv_gq=sv_gq, sv_dp=sv_dp, sv_max_dp=sv_max_dp,
                         sv_het_ab=sv_het_ab, sv_hom_ab=sv_hom_ab,
                         sv_min_control_dp=sv_min_control_dp,
                         sv_max_control_dp=sv_max_control_dp,
                         sv_min_control_gq=sv_min_control_gq,
                         sv_control_het_ab=sv_control_het_ab,
                         sv_control_hom_ab=sv_control_hom_ab,
                         sv_con_ref_ab=sv_con_ref_ab,
                         min_families=min_families,
                         report_file=report_file,)
        self.prefix = "gassoc_biallelic"
        self.header_fields = [("gassoc_biallelic_homozygous",
               '"Samples that carry homozygous biallelic changes ' +
               ' parsed by {}"' .format(type(self).__name__)),
               ("gassoc_biallelic_compound_het",
               '"Samples that carry compound heterozygous biallelic changes ' +
               'parsed by {}"'.format(type(self).__name__)),
               ("gassoc_biallelic_de_novo",
               '"Samples that carry biallelic alleles that appear to have ' +
               'arisen de novo"'),
                ('gassoc_biallelic_families',
                '"Family IDs for gassoc_biallelic alleles"'),
               ("gassoc_biallelic_features",
               '"Features (e.g. transcripts) that contain qualifying ' +
               'biallelic variants parsed by {}"' .format(
                type(self).__name__)),]
        self.annot_fields = ('homozygous', 'compound_het', 'de_novo',
                            'families', 'features')
        self.max_incidentals = max_incidentals

    def process_potential_recessives(self, assoc_alleles, final=False):
        '''
            Check whether stored PotentialSegregant alleles make up
            biallelic variation in the same transcript for affected
            individuals/families where at least one allele is in the
            assoc_alleles. Adds labels to INFO fields of VCF records and
            returns an OrderedDict of 'var_ids' to lists of
            PotentialSegregant objects that appear to segregate
            consistent with recessive inheritance.

            Clears the cache of stored PotentialSegregant alleles.

            Args:
                assoc_alleles:
                        list of alt_ids (matching the 'alt_id' property
                        of PotentialSegregant objects) that must be
                        present as at least one allele of biallelics.

                final:  if True, all cached variants will be processed.
                        Normal behaviour is not to process cached
                        variants for features present in the last cached
                        variant to ensure genes/transcripts are only
                        processed once all variants in that feature have
                        been encountered.

        '''
        segregating = OrderedDict() #keys are alt_ids, values are SegregatingBiallelic
        comp_hets = defaultdict(dict) # store compound hets per feature per
                                      # family for checking max_incidentals
        for feat, prs in self._potential_recessives.items():
            if not final and feat in self._last_added:
                continue
            feat_segregating = [] #list of tuples of values for creating SegregatingBiallelic
            un_hets = defaultdict(list)  #store het alleles carried by each unaffected
            aff_hets = defaultdict(list) #store het alleles carried by each affected
            biallelics = defaultdict(list)  #store biallelic combinations for affecteds
            for pid,p in prs.items():
                for un in self.unaffected:
                    if p.allele_counts[un] == 1:  #already checked for homs when adding
                        #store allele carried in this unaffected
                        un_hets[un].append(pid)
                for aff in (x for x in self.affected
                            if self.ped.fid_from_iid(x) in p.families):
                    if p.allele_counts[aff] == 1:
                        aff_hets[aff].append(pid)
                    elif p.allele_counts[aff] == 2:
                        if pid in assoc_alleles:
                            biallelics[aff].append(tuple([pid]))
            incompatibles = [] #create a list of sets of incompatible hets
            for hets in un_hets.values():
                if len(hets):
                    incompatibles.append(set(hets))
            for aff,hets in aff_hets.items():
                for i in range(len(hets)):
                    for j in range(i+1, len(hets)):
                        incomp = False
                        for iset in incompatibles:
                            if iset.issuperset([hets[i], hets[j]]):
                                incomp = True
                                break
                        if incomp:
                            continue
                        if prs[hets[i]].record.in_cis_with(
                            sample=aff,
                            allele=prs[hets[i]].allele,
                            other=prs[hets[j]].record,
                            other_allele=prs[hets[j]].allele):
                            #phase groups indicate alleles in cis
                            continue
                        if (hets[i] in assoc_alleles or hets[j] in
                                assoc_alleles):
                                biallelics[aff].append(tuple([hets[i],hets[j]]))
            if not biallelics:
                continue
            #see if all affecteds in the same family share the same biallelics
            for fid,affs in self._fam_to_aff.items():
                b_affs = set(x for x in affs if x in biallelics)
                if len(b_affs) == 0 or b_affs != affs:
                    continue
                affs = list(affs)
                absent_in_aff = False
                for i in range(len(affs)):
                    for bi in biallelics[affs[i]]:
                        for j in range(i+1, len(affs)):
                            if bi not in biallelics[affs[j]]:
                                absent_in_aff = True
                                break
                        if not absent_in_aff:
                            segs,de_novo = self._check_parents(feat, bi, affs)
                            if not segs:
                                continue
                            if len(bi) == 1:
                                model = 'homozygous'
                            else:
                                model = 'compound_het'
                                if fid not in comp_hets[feat]:
                                    comp_hets[feat][fid] = list()
                                comp_hets[feat][fid].append(bi)
                            for bi_pr in (prs[x] for x in bi):
                                feat_segregating.append((bi_pr, affs, [fid],
                                                         model, [feat],
                                                         de_novo[bi_pr.alt_id],
                                                         self.prefix))
            if self.max_incidentals:
                feat_segregating = self._filter_incidentals(feat_segregating,
                                                            comp_hets)
            fam_count = len(set([fam for tup in feat_segregating for fam in
                                 tup[2]]))
            if fam_count >= self.min_families:
                for tp in feat_segregating:
                    if tp[0] in segregating:
                        segregating[tp[0]].add_samples(*tp[1:6])
                    else:
                        segregating[tp[0]] = SegregatingVariant(*tp)
        var_to_segregants = OrderedDict()
        for sb in segregating.values():
            sb.annotate_record(self.report_file, self.annot_fields)
            if sb.segregant.var_id in var_to_segregants:
                var_to_segregants[sb.segregant.var_id].append(sb.segregant)
            else:
                var_to_segregants[sb.segregant.var_id] = [sb.segregant]
        #clear the cache except for the last entry which will be a new gene
        self._potential_recessives = self._last_added
        self._last_added = dict()
        return var_to_segregants

    def _filter_incidentals(self, feat_segregating, comp_hets):
        filtered_segregating = []
        filtered = False
        alt2affs = defaultdict(list)
        alt2fams = defaultdict(list)
        for ft in feat_segregating:
            #each alt may occur more than once under different model or family
            alt2affs[ft[0].alt_id].extend(ft[1])
            alt2fams[ft[0].alt_id].extend(ft[2])
        for ft in feat_segregating:
            incidental_affs = (k for k,v in ft[0].allele_counts.items() if k
                               not in alt2affs[ft[0].alt_id] and k in
                               self.ped.individuals and v)
            incidentals = set(self.ped.fid_from_iid(x) for x in incidental_affs
                              if self.ped.fid_from_iid(x) not in
                              alt2fams[ft[0].alt_id])
            if len(incidentals) < self.max_incidentals:
                filtered_segregating.append(ft)
            else:
                filtered = True
        if not filtered:
            return feat_segregating
        #reassess compound hets after filtering
        still_segregating = []
        hets = defaultdict(list)
        for ft in filtered_segregating:
            if ft[3] == 'homozygous':
                still_segregating.append(ft)
            else:
                for fam in ft[2]:
                    hets[fam].append(ft)
        for fam,het_list in hets.items():
            c_het_indices = set()
            for i in range(len(het_list)):
                combis = []
                for feat in het_list[i][4]:
                    combis.extend(comp_hets[feat][fam])
                for j in range(i+1, len(het_list)):
                    a_i = het_list[i][0].alt_id
                    a_j = het_list[j][0].alt_id
                    if (a_i, a_j) in combis or (a_j, a_i) in combis:
                        c_het_indices.add(i)
                        c_het_indices.add(j)
            for k in sorted(c_het_indices):
                still_segregating.append(het_list[k])
        return still_segregating


def main(args):
    variant_cache = VariantCache()
    vcfreader = VcfReader(args.vcf_input)
    avail_assoc_fields = get_assoc_fields(vcfreader)
    if args.pops:
        assoc_fields = ['gassoc_{}_P'.format(x) for x in args.pops]
        if not set(avail_assoc_fields).issuperset(set(assoc_fields)):
            sys.exit("ERROR: The following P-value annotations populations " +
                     "were not found in your VCF: " + ", ".join(
                         set(assoc_fields).difference(avail_assoc_fields)) +
                     ". Check populations provided to --pops.")
    else:
        assoc_fields = avail_assoc_fields
    freq_fields = []
    if args.freq:
        freq_fields = get_vase_freq_annots(vcfreader)
    allele_filters = [] #lambdas for filtering alleles for non-hit ALTs
    if args.cadd_phred:
        if 'CADD_PHRED_score' not in vcfreader.metadata['INFO']:
            sys.exit("No pre-exisiting CADD_PHRED_score INFO annotation in " +
                     "header - please annotate CADD scores using vase if you" +
                     " want to use --cadd_phred cutoffs.")
        cadd_filter = lambda x: filter_on_existing_cadd(x, args.cadd_phred)
        allele_filters.append(cadd_filter)
    pedfile = PedFile(args.ped)
    family_filter = FamilyFilter(ped=pedfile, vcf=vcfreader,
                                 logging_level=logger.level)
    assoc_segregator = AssocSegregator(family_filter,
                                       gq=args.gq,
                                       dp=args.dp,
                                       max_dp=args.max_dp,
                                       min_families=args.min_families,
                                       het_ab=args.het_ab,
                                       max_incidentals=args.max_incidentals)
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
    any_or_all = any if args.allow_any_pop else all
    write_vcf_header(vcfreader, vcf_writer, assoc_segregator)
    n,w = 0,0
    assoc_alts = list()
    for record in vcfreader:
        assocs = process_record(record, variant_cache=variant_cache,
                                pval=args.p_value,
                                min_alleles=args.min_alleles,
                                assoc_fields=assoc_fields,
                                csq_filter=csq_filter,
                                hit_vep_filter=no_csq_filter,
                                assoc_segregator=assoc_segregator,
                                freq=args.freq, freq_fields=freq_fields,
                                any_or_all=any_or_all,
                                allele_filters=allele_filters)
        assoc_alts.extend(assocs)
        if variant_cache.output_ready:
            w += process_cache(variant_cache, assoc_alts, assoc_segregator,
                               args.p_value, args.min_alleles, assoc_fields,
                               vcf_writer)
            variant_cache.output_ready.clear()
            assoc_alts = assoc_alts[-1:]
        n += 1
        if n % args.progress_interval == 0:
            update_progress(n, w, record, args.log_progress)
    w += process_cache(variant_cache, assoc_alts, assoc_segregator,
                       args.p_value, args.min_alleles, assoc_fields,
                       vcf_writer, final=True)
    variant_cache.output_ready.clear()
    update_progress(n, w, record, args.log_progress)
    if vcf_writer is not sys.stdout:
        vcf_writer.close()

def process_cache(cache, assoc_alts, assoc_segregator, pval, min_alleles,
                  assoc_fields, outfh, final=False):
    written = 0
    if final:
        cache.add_cache_to_output_ready()
    if not assoc_alts:
        return 0
    var_id_to_seg = assoc_segregator.process_potential_recessives(assoc_alts,
                                                                  final=final)
    for var in cache.output_ready:
        if var.can_output or var.var_id in var_id_to_seg:
            outfh.write(str(var.record) + '\n')
            written += 1
    return written

def process_record(record, variant_cache, pval, min_alleles, assoc_fields,
                   csq_filter, hit_vep_filter, assoc_segregator, freq,
                   freq_fields, any_or_all, allele_filters=[]):
    remove_alleles = [False] * (len(record.ALLELES) -1)
    for i in range(1, len(record.ALLELES)):
        if record.ALLELES[i] == '*':
            remove_alleles[i-1] = True
    if freq_fields and freq:
        r = filter_on_existing_freq(record, freq, freq_fields)
        set_true_if_true(remove_alleles, r)
    is_hit = is_assoc_hit(record, pval, min_alleles, assoc_fields, any_or_all)
    if any(is_hit):
        r_alts, h_remove_csq = hit_vep_filter.filter(record)
        set_true_if_true(remove_alleles, r_alts)
    if all(is_hit):
        remove_csq = h_remove_csq
    else:
        r_alts, remove_csq = csq_filter.filter(record)
        if any(is_hit):
            remove_alleles = [r_alts[i] if not is_hit[i] else remove_alleles[i]
                              for i in range(len(r_alts))]
            remove_csq = [remove_csq[i] if not
                          is_hit[record.CSQ[i]['alt_index']-1]
                          else h_remove_csq[i] for i in range(len(remove_csq))]
        else:
            set_true_if_true(remove_alleles, r_alts)
        for afilter in allele_filters:
            r_alts = afilter(record)
            if any(is_hit):
                remove_alleles = [r_alts[i] if not is_hit[i] else
                                  remove_alleles[i] for i in
                                  range(len(r_alts))]
            else:
                set_true_if_true(remove_alleles, r_alts)
    if all(remove_alleles) or all(remove_csq):
        return []
    segs = assoc_segregator.process_record(record, remove_alleles, remove_csq)
    # cache is checked here because it MUST only be checked after running
    # assoc_segregator.process_record or the two caches go out of sync
    assoc_alts = list()
    if segs:
        variant_cache.add_record(record)
        for i in (x for x in range(len(is_hit)) if is_hit[x]):
            alt_id = "{}:{}-{}/{}".format(record.CHROM, record.POS,
                                          record.REF, record.ALLELES[i+1])
            assoc_alts.append(alt_id)
    else:
        variant_cache.check_record(record)
    return assoc_alts

def is_assoc_hit(record, pval, min_alleles, assoc_fields, any_or_all=all):
    '''
        For each ALT allele check if they meet the requirement for an
        enrichment hit
    '''
    allele_is_hit = [False] * (len(record.ALLELES) -1)
    info = record.parsed_info_fields(fields=['gassoc_cohort_alt'] + assoc_fields)
    for i in range(len(record.ALLELES) - 1):
        if (info['gassoc_cohort_alt'][i] is None or
                info['gassoc_cohort_alt'][i] < min_alleles):
            continue
        if any_or_all(info[f][i] is not None and info[f][i] <= pval for f in
                      assoc_fields):
            allele_is_hit[i] = True
    return allele_is_hit

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

def filter_on_existing_cadd(record, score):
    remove  = [False] * (len(record.ALLELES) -1)
    phreds = record.parsed_info_fields(fields=['CADD_PHRED_score'])
    if 'CADD_PHRED_score' in phreds:
        for i in range(len(remove)):
            if (phreds['CADD_PHRED_score'][i] is not None and
                phreds['CADD_PHRED_score'][i] < score):
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
    annots = list()
    for info in vcf.metadata['INFO']:
        match = re.search('^gassoc_\w+_P', info)
        if match:
            annots.append(info)
            logger.debug("Identified association INFO field '{}'".format(info))
    if not annots:
        sys.exit("No P-value annotations found in VCF header. Exiting.\n")
    return annots

def update_progress(n, w, record, log=False):
    n_prog_string = ("{:,} variants processed, {:,} written ".format(n, w) +
                     "at {}:{}".format(record.CHROM, record.POS))
    global prog_string
    if log:
        logger.info(n_prog_string)
    else:
        n_prog_string = '\r' + n_prog_string
        if len(prog_string) > len(n_prog_string):
            sys.stderr.write('\r' + ' ' * len(prog_string) )
        sys.stderr.write(prog_string)
    prog_string = n_prog_string


def write_vcf_header(vcf, fh, assoc_seg):
    vcf.header.add_header_field(name="assoc_hits",
                               string='"' + str.join(" ", sys.argv) + '"')
    inf = dict()
    for f in assoc_seg.header_fields:
        inf[f[0]] = {'Number' : 'A', 'Type' : 'String',
                     'Description' : f[1] }
    for f,d in inf.items():
        vcf.header.add_header_field(name=f, dictionary=d, field_type='INFO')
    fh.write(str(vcf.header))

def get_options():
    biotype_default = ['3prime_overlapping_ncrna', 'antisense',
                        'CTCF_binding_site', 'enhancer', 'lincRNA',
                        'miRNA', 'misc_RNA', 'Mt_rRNA', 'Mt_tRNA',
                        'open_chromatin_region', 'polymorphic_pseudogene',
                        'processed_transcript', 'promoter',
                        'promoter_flanking_region', 'protein_coding', 'rRNA',
                        'sense_intronic', 'sense_overlapping', 'snoRNA',
                        'snRNA', 'TF_binding_site',
                        'translated_processed_pseudogene', 'TR_C_gene',
                        'TR_D_gene', 'TR_J_gene', 'TR_V_gene']
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
    parser.add_argument("--min_families", type=int, default=2, help=
                        '''Minimum number of families with biallelic hits in a
                        transcript. Default=2.''')
    parser.add_argument("--allow_any_pop", action='store_true',
                        help='''Require the p-value to be under threshold for
                        any population instead of requiring all analyzed
                        populations to have a p-value lower than the
                        threshold.''')
    parser.add_argument("--pops", nargs='+', help='''One or more gnomAD
                        populations to test against. Default to all
                        annotated p-values from gnomad_assoc.py.''')
    parser.add_argument('--freq', type=float, default=0.0, help=
                        '''Allele frequency cutoff. Frequency information will
                        be read from existing VASE and VEP annotations.
                        Variant alleles with an allele frequency equal to or
                        greater than this value will be filtered.''')
    parser.add_argument('--max_incidentals', type=int, default=0, help=
                        '''Maximum number of families with an individual that
                        carries one allele of a biallelic variant without
                        having a second qualifying variant. Useful for
                        filtering out implausibly common alleles that by chance
                        segregate in a subset of families. Default=0 (i.e. not
                        applied).''')
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
                        Variants that are not below the P-value cut-off which
                        do not result in one of these VEP consequence classes
                        will be filtered. Default=['TFBS_ablation',
                        'TFBS_amplification', 'inframe_deletion',
                        'inframe_insertion', 'frameshift_variant',
                        'initiator_codon_variant', 'missense_variant',
                        'protein_altering_variant',
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
    parser.add_argument('--biotypes',  nargs='+', default=biotype_default,
                        metavar='BIOTYPE', help='''Ignore consequences in
                        biotypes other than those specified here. Default =
                        {}'''.format(biotype_default))
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
    parser.add_argument('--cadd_phred', type=float, help='''CADD Phred score
                        cutoff for variants that are not below the P-value
                        threshold. Requires CADD scores to have been annotated
                        on your VCF using vase.''')
    parser.add_argument('--log_progress', action='store_true',
                        help='''Use logging output for progress rather than
                        wiping progress line after each update.''')
    parser.add_argument('--progress_interval', type=int, default=1000, metavar='N',
                        help='''Report progress information every N variants.
                        Default=1000.''')
    return parser

if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(args)
