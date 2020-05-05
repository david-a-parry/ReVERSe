# ReVERSe

**R**ar**e** **V**ariant **E**nrichment and **R**ecessive **Se**gregation

## Installation
 
### Method 1:

ReVERSe requires several modules from VASE (https://github.com/david-a-parry/vase).
First install VASE using pip (or for other methods look at the instructions in
the README at the [VASE github repo](https://github.com/david-a-parry/vase)) 
and then install ReVERSe from either the gitlab or github repo.

    python3 -m pip install git+git://github.com/david-a-parry/vase.git --user
    python3 -m pip install git+https://git.ecdf.ed.ac.uk/dparry/reverse.git  --user

The --user flag is optional, but circumvents permissions issues that may arise
depending on how your python libraries are configured.

### Method 2:

First clone the repo and install with dependancies.

    git clone https://git.ecdf.ed.ac.uk/dparry/reverse.git
    cd reverse
    python3 -m pip install . --process-dependency-links --user

If you get an error "no such option: --process-dependency-links" (with newer
versions of pip) with the third command see method 3.

### Method 3:

Note that the first two steps are not necssary if you've already run them while
attempting Method 2.

    git clone https://git.ecdf.ed.ac.uk/dparry/reverse.git
    cd reverse
    python3 -m pip install -r requirements.txt --user
    python3 -m pip install . --user

## Synopsis

    # Get counts from gnomAD data - make sure you are using the same genome 
    # build as the gnomAD VCF!

    ReVERSe_count -i cohort.vcf.gz \
        -ped cohort.ped \
        -e gnomad.exomes.r2.1.sites.vcf.gz \
        -g gnomad.genomes.r2.1.1.sites.vcf.bgz \
        --exome_coverage_files gnomad.exomes.coverage.summary.tsv.bgz \
        --genome_coverage_files gnomad.genomes.coverage.summary.tsv.bgz \
        --gnomad_version 2.1 \
        -v cohort_rev_counts.vcf.gz \
        -t cohort.rev_counts.txt.gz 
    
    # Perform recessive segregation analysis requiring either biallelic
    # enriched alleles or enriched alleles in trans with VEP HIGH impact
    # variants, looking at variants with a maximum frequency of 0.5 %. Alleles
    # with a Fisher's exact test P-value <= 1e-5 are considered enriched.
    
    ReVERSe_seg -i cohort_rev_counts.vcf.gz \
        --ped cohort.ped \
        --freq 0.005 \
        --canonical \
        --pops afr amr eas fin nfe sas \
        -v 1e-5 \
        --impact HIGH \
        -o cohort_rev_seg.vcf.gz

    # Output ranked variants in CSV format
    
    ReVERSe_reporter cohort_rev_seg.vcf.gz \
        cohort_rev_seg.report.csv \
        cohort.ped

Note that each of the programs have many different options that can be provided
to tweak parameters. For detailed help run the desired program with the --help
flag.

gnomAD VCFs and coverage data are available from https://gnomad.broadinstitute.org/downloads#v2-coverage
