#!/usr/bin/env python3
import sys
import argparse
import pandas as pd


def get_options():
    parser = argparse.ArgumentParser(
        description='''Convert consolidated OGEE results into a table readable
        by ReVERSe reporter.''',
        epilog='''
Download consolidated results table from here:

        http://ogee.medgenius.info/browse/Homo%20sapiens

Use the downloaded file as input to this script.
        '''
    )
    parser.add_argument('csv',
                        help='''Input csv file. This should be the
                        consolidated results for your species of interest
                        from OGEEv2.''')
    parser.add_argument('-o', '--output', default=sys.stdout,
                        help='''Output file name. If not specified output will
                        be written to STDOUT''')
    return parser


def parse_table(csv, output=sys.stdout):
    out_cols = ['locus', 'symbols', 'NE', 'E', 'Frac_E']
    ogee_df = pd.read_csv(csv, sep='\t')
    ogee_df['NE'] = ogee_df['essentiality status'].apply(
        lambda x: x.count("NE"))
    ogee_df['E'] = ogee_df['essentiality status'].apply(
        lambda x: 1 + x.count(",") - x.count("NE"))
    ogee_df['Frac_E'] = ogee_df['E']/(ogee_df['E'] + ogee_df['NE'])
    ogee_df.sort_values(by='Frac_E', ascending=False)
    ogee_df[out_cols].to_csv(output, index=False)


if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    parse_table(**vars(args))
