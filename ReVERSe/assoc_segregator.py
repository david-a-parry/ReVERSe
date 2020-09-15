from collections import defaultdict, OrderedDict
from vase.family_filter import RecessiveFilter
from vase.family_filter import SegregatingVariant


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

    def __init__(self, family_filter, gt_args, min_families=1, strict=False,
                 exclude_denovo=False, report_file=None, max_incidentals=0):
        super().__init__(family_filter, gt_args, min_families=min_families,
                         report_file=report_file,)
        self.prefix = "ReVERSe_biallelic"
        self.header_fields = [
            ("ReVERSe_biallelic_homozygous",
             'Samples that carry homozygous biallelic changes ' +
             ' parsed by {}' .format(type(self).__name__)),
            ("ReVERSe_biallelic_compound_het",
             'Samples that carry compound heterozygous biallelic changes ' +
             'parsed by {}'.format(type(self).__name__)),
            ("ReVERSe_biallelic_de_novo",
             'Samples that carry biallelic alleles that appear to have ' +
             'arisen de novo'),
            ('ReVERSe_biallelic_families',
             'Family IDs for ReVERSe_biallelic alleles'),
            ("ReVERSe_biallelic_features",
             'Features (e.g. transcripts) that contain qualifying ' +
             'biallelic variants parsed by {}' .format(
                 type(self).__name__))]
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
        segregating = OrderedDict()
        # keys are alt_ids, values are SegregatingBiallelic
        comp_hets = defaultdict(dict)
        # compound hets per feature per family for checking max_incidentals
        bi_and_carrier_fams = dict()
        # keys are feats, values are dicts indicating biallelic/carrier fams
        for feat, prs in self._potential_recessives.items():
            if not final and feat in self._current_features:
                continue
            feat_segregating = []  # tuples of values for SegregatingBiallelic
            un_hets = defaultdict(list)  # het alleles carried by unaffected
            aff_hets = defaultdict(list)  # het alleles carried by affected
            biallelics = defaultdict(list)  # biallelic combinations for affs
            phase_known_biallelics = set()
            phase_unknown_biallelics = set()
            carrier_fams = set()
            for pid, p in prs.items():
                for un in self.unaffected:
                    if p.allele_counts[un] == 1:
                        # checked for homs when adding store allele carried in
                        # this unaffected
                        un_hets[un].append(pid)
                        carrier_fams.add(self.ped.fid_from_iid(un))
                for aff in (x for x in self.affected
                            if self.ped.fid_from_iid(x) in p.families):
                    if p.allele_counts[aff] == 1:
                        aff_hets[aff].append(pid)
                        carrier_fams.add(self.ped.fid_from_iid(aff))
                    elif p.allele_counts[aff] == 2:
                        if pid in assoc_alleles:
                            biallelics[aff].append(tuple([pid]))
            incompatibles = []  # create a list of sets of incompatible hets
            for hets in un_hets.values():
                if len(hets):
                    incompatibles.append(set(hets))
            for aff, hets in aff_hets.items():
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
                            # phase groups indicate alleles in cis
                            continue
                        if (hets[i] in assoc_alleles or hets[j] in
                                assoc_alleles):
                            biallelics[aff].append(tuple([hets[i], hets[j]]))
            if not biallelics:
                continue
            # see if all affecteds in the same family share the same biallelics
            for fid, affs in self._fam_to_aff.items():
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
                            segs, de_novo = self._check_parents(feat, bi, affs)
                            if not segs:
                                continue
                            if len(bi) == 1:
                                model = 'homozygous'
                                phase_known_biallelics.add(fid)
                            else:
                                model = 'compound_het'
                                if any(self.ped.individuals[x].parents for x in
                                       affs if x in self.samples):
                                    if all(any([y in x for x in
                                                de_novo.values()]) for y in
                                           affs):  # phase of de novos unknown
                                        phase_unknown_biallelics.add(fid)
                                    else:
                                        phase_known_biallelics.add(fid)
                                else:  # no parents, phase unknown
                                    phase_unknown_biallelics.add(fid)
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
            seg_fams = set([fam for tup in feat_segregating for fam in
                            tup[2]])
            fam_count = len(seg_fams)
            if fam_count >= self.min_families:
                for tp in feat_segregating:
                    if tp[0] in segregating:
                        segregating[tp[0]].add_samples(*tp[1:6])
                    else:
                        segregating[tp[0]] = SegregatingVariant(*tp)
            carrier_fams = carrier_fams.difference(seg_fams)
            bi_and_carrier_fams[feat] = {'phased': phase_known_biallelics,
                                         'unphased': phase_unknown_biallelics,
                                         'carrier': carrier_fams}
        var_to_segregants = OrderedDict()
        for sb in segregating.values():
            sb.annotate_record(self.report_file, self.annot_fields)
            self._annotate_phase_and_carrier_info(sb, bi_and_carrier_fams)
            if sb.segregant.var_id in var_to_segregants:
                var_to_segregants[sb.segregant.var_id].append(sb.segregant)
            else:
                var_to_segregants[sb.segregant.var_id] = [sb.segregant]
        # clear the cache except for the last entry which will be a new gene
        self._potential_recessives = OrderedDict(
            (k, v) for k, v in self._potential_recessives.items() if k in
            self._current_features)
        return var_to_segregants

    def _filter_incidentals(self, feat_segregating, comp_hets):
        filtered_segregating = []
        filtered = False
        alt2affs = defaultdict(list)
        alt2fams = defaultdict(list)
        for ft in feat_segregating:
            # each alt may occur more than once under different model or family
            alt2affs[ft[0].alt_id].extend(ft[1])
            alt2fams[ft[0].alt_id].extend(ft[2])
        for ft in feat_segregating:
            incidental_affs = (k for k, v in ft[0].allele_counts.items() if k
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
        # reassess compound hets after filtering
        still_segregating = []
        hets = defaultdict(list)
        for ft in filtered_segregating:
            if ft[3] == 'homozygous':
                still_segregating.append(ft)
            else:
                for fam in ft[2]:
                    hets[fam].append(ft)
        for fam, het_list in hets.items():
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

    def _annotate_phase_and_carrier_info(self, sb, bi_fams):
        info = defaultdict(list)
        for feat in sorted(sb.features):
            for k, v in bi_fams[feat].items():
                f = 'ReVERSe_' + k + '_families'
                info[f].append('|'.join(v) or '.')
        sb.segregant.record.add_info_fields(info)
