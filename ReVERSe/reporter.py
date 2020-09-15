import os
from collections import defaultdict
from tempfile import NamedTemporaryFile
import pandas as pd
from vase.vase_reporter import VaseReporter
from vase.vcf_reader import VcfReader

supported_formats = ['xlsx', 'csv', 'tsv']


class ReverseReporter(VaseReporter):
    ''' Read a ReVERSe_seg annotated VCF and output summary data.'''

    def __init__(self, vcf, out, output_format='csv', ogee_csv=None,
                 info_fields=[], **kwargs):
        tmpf = NamedTemporaryFile(delete=False)
        tmpf.close()
        self.table_out = out
        self.table_format = output_format
        if self.table_format not in supported_formats:
            raise ValueError("--output_type  must be one of the following: " +
                             ",".join(supported_formats))
        if not self.table_out.endswith('.' + self.table_format):
            self.table_out += '.' + self.table_format
        if not ogee_csv:
            ogee_csv = os.path.join(os.path.dirname(__file__),
                                    "data",
                                    "ogeev2_frac_ess.csv")
        self.ogee_df = pd.read_csv(ogee_csv)
        feat_annots = {
            'ReVERSe_biallelic_families': 'ReVERSe_biallelic_features'}
        vcfreader = VcfReader(vcf)
        info_fields.extend(self._get_count_fields(vcfreader))
        super().__init__(vcfreader, out=tmpf.name, output_type='json',
                         info_fields=info_fields,
                         custom_feat_annots=feat_annots, **kwargs)

    def _get_count_fields(self, vcf):
        f = [x for x in vcf.header.info if 'ReVERSe' in x and 'biallelic'
             not in x]
        if not f:
            raise ValueError("No recognised ReVERSe count fields identified " +
                             "in input VCF ({}). ".format(vcf.filename) +
                             "Was ReVERSe_count run for this VCF?")
        return f

    def _get_seg_fields(self):
        inheritance_fields = dict()
        if 'ReVERSe_biallelic_families' in self.vcf.header.info:
            inheritance_fields['ReVERSe_biallelic_families'] = 'recessive'
            self.logger.info("Found ReVERSe biallelic annotations.")
        else:
            raise RuntimeError("no ReVERSe recessive biallelic annotations " +
                               "found in vcf. Please run ReVERSe_seg first.")
        return inheritance_fields

    def convert_seg_data(self):
        df = defaultdict(list)
        rep_cols = self._get_header_columns(family=None)
        for fam in self.json_dict:
            affected = list()
            mother = list()
            father = list()
            for aff in self.ped.families[fam].get_affected():
                affected.append(aff)
                if self.ped.individuals[aff].father:
                    father.append(self.ped.individuals[aff].father)
                if self.ped.individuals[aff].mother:
                    mother.append(self.ped.individuals[aff].mother)
            for var in self.json_dict[fam]:
                df['Family'].append(fam)
                df['Affected'].append("|".join(affected))
                df['Father'].append("|".join(father))
                df['Mother'].append("|".join(mother))
                for c in rep_cols:
                    if c in var:
                        df[c].append(var[c])
                    else:
                        df[c].append('')
                aff_gt = []
                father_gt = []
                mother_gt = []
                for aff in affected:
                    if aff in var:
                        aff_gt.append(var[aff])
                for fat in father:
                    if fat in var:
                        father_gt.append(var[fat])
                for mot in mother:
                    if mot in var:
                        mother_gt.append(var[mot])
                df['Affected_GT'].append("|".join(aff_gt))
                df['Father_GT'].append("|".join(father_gt))
                df['Mother_GT'].append("|".join(mother_gt))
        return pd.DataFrame.from_dict(df)

    def row_pop_max_p(self, row):
        pvals = [getattr(row, x) for x in self.p_cols]
        i = pvals.index(row.Max_P)
        return self.p_cols[i]

    def post_process_df(self, df):
        second_hits = []
        n_second_hits = []
        self.p_cols = [x for x in df.columns if x.endswith('_P')]
#        df['Family'] = df.Family.astype(int)
        for row in df.itertuples():
            gt = row.Affected_GT.split(":")[0].split("/")
            hit_2 = []
            if len(set(gt)) > 1:
                mask = (df.Feature == row.Feature) & (df.Affected ==
                                                      row.Affected)
                other_hits = df[mask][["HGVSc", "HGVSp"]]
                for orow in other_hits.itertuples():
                    if orow.HGVSc == row.HGVSc:
                        continue
                    h = orow.HGVSc
                    if orow.HGVSp:
                        h += "|{}".format(orow.HGVSp)
                    hit_2.append(h)
            else:
                hit_2.append("homozygous")
            second_hits.append(",".join(hit_2))
            n_second_hits.append(len(hit_2))
        df['Max_P'] = df.apply(lambda x: max([getattr(x, a) for a in
                                              self.p_cols]), axis=1)
        df['Pop_Max_P'] = df.apply(lambda x: self.row_pop_max_p(x), axis=1)
        df['Second_hit'] = second_hits
        df['N_second_hit'] = n_second_hits
        df['Carrier_Families'] = df.ReVERSe_carrier_families.apply(
            lambda x: len(x.split(',')[0].split('|')) if x.split(',')[0] != '.'
            else 0)
        df['N_Families'] = df.SYMBOL.apply(
            lambda x: len(df[df.SYMBOL == x].Family.unique()))
        df['Frac_E'] = df.SYMBOL.apply(
            lambda x:
                self.ogee_df[self.ogee_df.symbols == x]['Frac_E'].values[0]
            if x in self.ogee_df.symbols.values else 0)
        df['Rank'] = ((1 + df.Carrier_Families) *
                      (1.1 - df.Frac_E) *
                      (1 + df.N_second_hit))
        df = df.sort_values(by=['Rank', 'Max_P'])
        return df

    def _finish_up(self):
        df = self.convert_seg_data()
        df = self.post_process_df(df)
        if self.table_format == 'csv':
            df.to_csv(self.table_out, index=False)
        elif self.table_format == 'tsv':
            df.to_csv(self.table_out, index=False, sep='\t')
        elif self.table_format == 'xlsx':
            df.to_excel(self.table_out, index=False)
        self.out_fh.close()
