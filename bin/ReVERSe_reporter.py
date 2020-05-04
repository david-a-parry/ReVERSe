#!/usr/bin/env python3
import argparse
from ReVERSe.reporter import ReverseReporter


def parse_args():
    parser = argparse.ArgumentParser(
                   description='''Write per family reports from a ReVERSe_seg
                                  annotated VCF.''')
    parser.add_argument("vcf", help='''ReVERSe annotated VCF file.''')
    parser.add_argument("out", help='''Name for output XLSX/JSON file.''')
    parser.add_argument("ped", help='''PED file (same as used with ReVERSe
                        segregation analysis).''')
    parser.add_argument("-f", "--families", nargs='+',
                        help='''One or more families to report variants
                        for.''')
    parser.add_argument("--blacklist", help='''A file containing a list of
                        Ensembl feature IDs to ignore, one ID per line.''')
    parser.add_argument("--choose_transcript", action='store_true', help=
                        '''Pick a single transcript for each variant. The
                        highest impact transcript will be chosen with canonical
                        and protein coding transcripts preferred where impacts
                        are the same.''')
    parser.add_argument("-g", "--g2p", help='''G2P CSV file - if provided
                        additional fields will be added for G2P genes.''')
    parser.add_argument("--filter_non_g2p", action='store_true', help=
                        '''If using a G2P file, only output records for genes
                        present in the G2P file.''')
    parser.add_argument("--allelic_requirement", action='store_true', help=
                        '''If using --filter_non_g2p, only output records if
                        the inheritance pattern is consistent with the
                        'allelic requirement' annotation from G2P.''')
    parser.add_argument("--mutation_requirement", action='store_true', help=
                        '''If using --filter_non_g2p, only output records if
                        the variant consequence is consistent with the
                        'mutation consequence' annotation from G2P.''')
    parser.add_argument("-c", "--gnomad_constraint", help='''Constraint txt
                        file from gnomAD. If provided, pLI, pRec, pNull, mis_z
                        and syn_z columns will be added from this file for
                        matching transcripts/genes.''')
    parser.add_argument("-l", "--mygene_lookups", action="store_true", help=
                        '''Add additional fields by connecting to MyGene.info
                        service. This provides the following extra columns:
                        ENTREZ_ID, Name, Summary, GO_BP, GO_CC, GO_MF, MIM and
                        GeneRIFs.''')
    parser.add_argument("-r", "--rest_lookups", action="store_true", help=
                        '''Add additional fields by connecting to Ensembl's
                        REST server. Not recommended for large numbers of
                        variants as this will cause significant slowdown. The
                        following extra columns will be added: ENTREZ,
                        Full_Name, GO, REACTOME, MOUSE_TRAITS, MIM_MORBID.''')
    parser.add_argument("--grch37", action="store_true",
                        help='''Use GRCh37 REST server
                        (http://grch37.rest.ensembl.org/) for REST lookups.''')
    parser.add_argument("-i", "--info_fields", nargs='+', default=[], help='''
                        Additional INFO field(s) to add as a column in your
                        output. Fields that are annotated per allele
                        (Number=A or Number=R) will be annotated for the
                        relevant allele otherwise the whole INFO field will be
                        added.''')
    parser.add_argument("-t", "--timeout", type=float, default=2.0, help=
                        '''Timeout (in seconds) for REST lookups.
                        Default=2.0''')
    parser.add_argument("-m", "--max_retries", type=int, default=2, help=
                        '''Number of reattempts for REST lookups that fail.
                        Default=2''')
    parser.add_argument("--prog_interval", type=int, metavar="N", help=
                        '''Report progress every N variants. Defaults to 1000
                        unless using --rest_lookups in which case it defaults
                        to 100.''')
    parser.add_argument('--hide_empty', action='store_true', help=
                        '''Hide worksheets that do not contain any variants in
                        the output.''')
    parser.add_argument('--quiet', action='store_true', help=
                        '''Do not output progress information to STDERR.''')
    parser.add_argument('--debug', action='store_true', help=
                        '''Output debugging information to STDERR.''')
    parser.add_argument('--force', action='store_true', help=
                        '''Force overwrite of existing output files.''')
    return parser


if __name__ == '__main__':
    parser = parse_args()
    args = parser.parse_args()
    runner = ReverseReporter(**vars(args))
    try:
        runner.write_report()
    finally:
        runner.out_fh.close()
