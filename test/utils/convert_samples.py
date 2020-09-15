#!/usr/bin/env python3
import sys
import os
import gzip


def main(ped, vcf):
    fams = dict()
    f_count, i_count = 0, 0
    indvs = {'0': '0'}
    ped_out = os.path.splitext(ped)[0] + '.renamed.ped'
    with open(ped, 'rt') as infile, open(ped_out, 'wt') as outfile:
        for line in infile:
            fid, *rest = line.split()
            if rest[0] == '0':
                sys.exit("Individual ID can not be '0' - exiting\n")
            if fid not in fams:
                f_count += 1
                fams[fid] = "F{}".format(f_count)
            converted = [fams[fid]]
            for s in rest[:3]:
                if s not in indvs:
                    i_count += 1
                    indvs[s] = "I{}".format(i_count)
                converted.append(indvs[s])
            outfile.write("\t".join(converted + rest[3:]) + '\n')
    if vcf.endswith(".gz"):
        op = gzip.open
        vcf_out = vcf.replace(".vcf.gz", "") + ".renamed.vcf.gz"
    else:
        op = open
        vcf_out = os.path.splitext(vcf)[0] + '.renamed.vcf'
    with op(vcf, 'rt') as infile, op(vcf_out, 'wt') as outfile:
        for line in infile:
            if line.startswith('#CHROM'):
                cols = line.split()
                new_head = cols[:9]
                samples = cols[9:]
                for s in samples:
                    if s not in indvs:
                        i_count += 1
                        indvs[s] = "I{}".format(i_count)
                    new_head.append(indvs[s])
                outfile.write("\t".join(new_head) + '\n')
            else:
                outfile.write(line)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("Usage: {} input.ped input.vcf.gz".format(sys.argv[0]))
    main(sys.argv[1], sys.argv[2])
