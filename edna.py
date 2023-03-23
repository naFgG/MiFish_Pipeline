#!/usr/bin/env python
__author__ = 'FAN Gong'
__version__ = '25 October 2022, created'
__coding__ = 'utf-8'

import os
import pandas as pd
import numpy as np
import glob
from multiprocessing import Pool, Process
import argparse

parser = argparse.ArgumentParser(usage='python3 <PATH/to/this/script>/edna.py [OPTIONS]',
                                 description='eDNA pipeline')
parser.add_argument('--raw', help='Directory contains the zipped raw reads', required=False)
parser.add_argument('--primer', help='the primers, samples txt', required=False, default='primer.txt')
parser.add_argument('--qc', help='Quality control', required=False, action='store_true')
parser.add_argument('--trim', help='trimming', required=False, action='store_true')
parser.add_argument('--assemble', help='if need to assemble', required=False, action='store_true')
parser.add_argument('--mapping', help='if need to map', required=False, action='store_true')
parser.add_argument('--zotu', help='denoise, chimera filtration and zotu clustering', required=False, action='store_true')
parser.add_argument('--otu', help='chimera filtration and otu clustering', required=False, action='store_true')
parser.add_argument('--ref', help='Path of reference for mapping', required=False, default=0)
parser.add_argument('--zblast', help='blastn zotu to database', required=False, action='store_true')
parser.add_argument('--blast', help='blastn otu to database', required=False, action='store_true')
parser.add_argument('--zfeature', help='output final zotu result', required=False, action='store_true')
parser.add_argument('--feature', help='output final otu result', required=False, action='store_true')
args = parser.parse_args()


def unzip(gzs):
    """
    Unzipping the raw reads
    :param gzs: zipped raw reads
    :return: 0
    """
    os.system(f"gunzip -c {gzs} > ../unzipped/{gzs.split('.gz')[0]}")
    return 0


def qc(sample):
    """
    Quality control
    :param sample: unzipped raw reads
    :return: 0
    """
    os.system(f"{fastp} -i {sample}.R1.fq -I {sample}.R2.fq "
              f"-o ../cleaned/{sample}.R1.fq -O ../cleaned/{sample}.R2.fq "
              f"-j {sample}.json "
              f"-h {sample}.html "
              f"-A "
              f"-q 30 "
              f"-n 5 "
              f"-x 10 "
              f"-l 15 "
              f"-w 4")
    return 0


def assembling(fq):
    """
    Merging trimmed reads
    :param fq: qced
    :return: 0
    """
    os.system(f'{flash} -M 80 -t 4 -o {fq} '
              f'../cleaned/{fq}.R1.fq '
              f'../cleaned/{fq}.R2.fq')
    os.system(f'mv {fq}.extendedFrags.fastq {fq}.fq')
    return 0


def trimming():
    """
    Trimming primer
    :return: 0
    """
    os.chdir('merged')
    os.system('for i in `ls *.fq`;do ln -s %s/merged/${i} %s/trimmed/${i};done' % (submit, submit))
    os.chdir(f'{submit}/trimmed')

    # def trim1(sample, pf):
    #     os.system(f'perl {tagcleaner} -64 '
    #               f'-fastq {sample}.R1.fq --out_format 3 -minlen 30 '
    #               f'-tag5 {pf} --verbose -out {sample}_tagcleaned.R1')
    #     return 0
    #
    # def trim2(sample, pr):
    #     os.system(f'perl {tagcleaner} -64 '
    #               f'-fastq {sample}.R2.fq --out_format 3 -minlen 30 '
    #               f'-tag5 {pr} --verbose -out {sample}_tagcleaned.R2')
    #     return 0
    base_dict = {'A': 'T',
                 'T': 'A',
                 'C': 'G',
                 'G': 'C',
                 'N': 'N'}

    if primer['primer'].str.contains('mix').sum() >= 1:
        mix = primer[primer['primer'].str.contains('mix')]
        mix.reset_index(drop=True, inplace=True)
        for i in mix.index:
            samples = mix.iat[i, 1].split(',')
            primersF = mix.iat[i, 2].split(',')
            primersR = mix.iat[i, 3].split(',')
            for j in range(len(samples)):
                for k in range(len(primersF)):
                    pR = ''.join([base_dict[i] for i in primersR[k][::-1]])
                    if k == 0:
                        os.system(f'perl {tagcleaner} -64 '
                                  f'-fastq {samples[j]}.fq --out_format 3 -minlen 30 '
                                  f'-mm5 2 -mm3 2 '
                                  f'-tag5 {primersF[k]} '
                                  f'-tag3 {pR} '
                                  f'--verbose -out {samples[j]}_{k}')
                    else:
                        os.system(f'perl {tagcleaner} -64 '
                                  f'-fastq {samples[j]}_{k-1}.fastq --out_format 3 -minlen 30 '
                                  f'-mm5 2 -mm3 2 '
                                  f'-tag5 {primersF[k]} '
                                  f'-tag3 {pR} '
                                  f'--verbose -out {samples[j]}_{k}')
                        os.system(f'rm {samples[j]}_{k - 1}.fastq')
                        if k == (len(primersF) - 1):
                            os.system(f'mv {samples[j]}_{k}.fastq {samples[j]}_tagcleaned.fastq')
                    # if k == 0:
                    #     os.system(f'perl {tagcleaner} -64 '
                    #               f'-fastq {samples[j]}.R1.fq --out_format 3 -minlen 30 '
                    #               f'-tag5 {primersF[k]} --verbose -out {samples[j]}_{k}.R1')
                    #     os.system(f'perl {tagcleaner} -64 '
                    #               f'-fastq {samples[j]}.R2.fq --out_format 3 -minlen 30 '
                    #               f'-tag5 {primersR[k]} --verbose -out {samples[j]}_{k}.R2')
                    # else:
                    #     os.system(f'perl {tagcleaner} -64 '
                    #               f'-fastq {samples[j]}_{k-1}.R1.fastq --out_format 3 -minlen 30 '
                    #               f'-tag5 {primersF[k]} --verbose -out {samples[j]}_{k}.R1')
                    #     os.system(f'perl {tagcleaner} -64 '
                    #               f'-fastq {samples[j]}_{k-1}.R2.fastq --out_format 3 -minlen 30 '
                    #               f'-tag5 {primersR[k]} --verbose -out {samples[j]}_{k}.R2')
                    #     os.system(f'rm {samples[j]}_{k - 1}.R*.fastq')
                    #     if k == (len(primersF) - 1):
                    #         os.system(f'mv {samples[j]}_{k}.R1.fastq {samples[j]}_tagcleaned.R1.fastq')
                    #         os.system(f'mv {samples[j]}_{k}.R2.fastq {samples[j]}_tagcleaned.R2.fastq')
        unmix = primer[~primer['primer'].str.contains('mix')]
        unmix.index = range(len(unmix))
        for i in unmix.index:
            samples = unmix.iat[i, 1].split(',')
            primersF = unmix.iat[i, 2]
            primersR = unmix.iat[i, 3]
            for j in samples:
                pR = ''.join([base_dict[i] for i in primersR[::-1]])
                os.system(f'perl {tagcleaner} -64 '
                          f'-fastq {j}.fq --out_format 3 -minlen 30 '
                          f'-mm5 2 -mm3 2 '
                          f'-tag5 {primersF} '
                          f'-tag3 {pR} '
                          f'--verbose -out {j}_tagcleaned')
               # did not filter the reads not matching the tag at 5'-end
               # os.popen(f'perl {tagcleaner} -64'
               #        f'-fastq {j}_R1.fq --out_format 3 -minlen 30 '
               #        f'-tag5 {primersF} --verbose -out {j}_tagcleaned_R1')
               # filter the reads not matching the tag at 5'-end
               # os.popen(f'perl {tagcleaner} -64'
               #        f'-fastq {j}_R1.fq --out_format 3 -minlen 30 '
               #        f'-tag5 {primersF} -nomatch 1 --verbose -out {j}_tagcleaned_R1')
               #  t1 = Process(target=trim1, args=(j, primersF))
               #  t2 = Process(target=trim2, args=(j, primersR))
               #  t1.start()
               #  t2.start()
               #  t1.join()
               #  t2.join()
    else:
        for i in primer.index:
            samples = primer.iat[i, 1].split(',')
            primersF = primer.iat[i, 2]
            primersR = primer.iat[i, 3]
            for j in samples:
                pR = ''.join([base_dict[i] for i in primersR[::-1]])
                os.system(f'perl {tagcleaner} -64 '
                          f'-fastq {j}.fq --out_format 3 -minlen 30 '
                          f'-mm5 2 -mm3 2 '
                          f'-tag5 {primersF} '
                          f'-tag3 {pR} '
                          f'--verbose -out {j}_tagcleaned')
                # t1 = Process(target=trim1, args=(j, primersF))
                # t2 = Process(target=trim2, args=(j, primersR))
                # t1.start()
                # t2.start()
                # t1.join()
                # t2.join()
    os.system('for i in `find ./ -type l`;do rm ${i};done')
    os.chdir(submit)
    return 0


def mapping(fq):
    """
    Trimmed reads map to a reference
    :param fq: trimmed reads
    :return: 0
    """
    os.system(f'bwa mem {ref} {fq}_tagcleaned.fastq > {fq}.sam')
    # Extracting reads which mapped to the ref
    # os.system(f'samtools view -bF 4 {fq}.sam > {fq}.bam')
    # Converting all reads
    os.system(f'samtools view -bS {fq}.sam > {fq}.bam')
    os.system(f'samtools sort {fq}.bam > {fq}.sorted.bam')
    os.system(f'samtools index {fq}.sorted.bam')
    # Statistics
    os.system(f'samtools flagstat {fq}.sorted.bam > {fq}_stat.txt')
    return 0


def usearching(fq):
    """
    Denoise, chimera filtration and zotu clustering using unoise3 algorithm.
    :param fq: trimmed reads
    :return: 0
    """
    # without pooling samples
    # Quality filter
    os.system(f'{usearch} '
              f'-fastq_filter ../trimmed/{fq}_tagcleaned.fastq '
              f'-fastq_maxee 1.0 '
              f'-fastaout {fq}_filtered.fa '
              f'-threads 4')
    # Find unique read sequences and annotate abundances
    os.system(f'{usearch} '
              f'-fastx_uniques {fq}_filtered.fa '
              f'-sizeout '
              f'-relabel Uniq '
              f'-fastaout {fq}_uniques.fa')
    if args.zotu:
        # Denoise: predict biological sequences and filter chimeras
        os.system(f'{usearch} '
                  f'-unoise3 {fq}_uniques.fa '
                  f'-zotus {fq}_zotus.fa')
        os.system(f'{usearch} '
                  f'-otutab ../trimmed/{fq}_tagcleaned.fastq '
                  f'-zotus {fq}_zotus.fa '
                  f'-otutabout {fq}_zotu.xls '
                  f'-threads 4')
                 # -mapout zmap.txt
        os.system(f"sed -i '1d' {fq}_zotu.xls && sed -i '1i #OTU_ID\t{fq}' {fq}_zotu.xls")
    if args.otu:
        # Make 97% OTUs and filter chimeras
        os.system(f'{usearch} '
                  f'-cluster_otus {fq}_uniques.fa '
                  f'-otus {fq}_otus.fa '
                  f'-relabel Otu')
        os.system(f'{usearch} '
                  f'-otutab ../trimmed/{fq}_tagcleaned.fastq '
                  f'-otus {fq}_otus.fa '
                  f'-otutabout {fq}_otu.xls '
                  f'-threads 4')
        os.system(f"sed -i '1d' {fq}_otu.xls && sed -i '1i #OTU_ID\t{fq}' {fq}_otu.xls")
    return 0


def blastn(fq):
    """
    Blast
    :param fq: samples names
    :return: 0
    """
    print(f'###### start {fq} blasting ######')
    if args.zblast:
        os.system(f'blastn '
                  f'-query {fq}_zotus.fa '
                  f'-db {mitofish} '
                  f'-task blastn '
                  f'-outfmt 6 '
                  f'-num_alignments 1 '
                  f'-num_threads 4 '
                  f'-out {fq}_aln.xls')
    if args.blast:
        os.system(f'blastn '
                  f'-query {fq}_otus.fa '
                  f'-db {mitofish} '
                  f'-task blastn '
                  f'-outfmt 6 '
                  f'-num_alignments 1 '
                  f'-num_threads 4 '
                  f'-out {fq}_aln.xls')
    os.system(f"sed -i '1i query id\tsubject id\t% identity\talignment length\t"
              f"mismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue"
              f"\tbit score' {fq}_aln.xls")
    return 0


def featuring(fq):
    """
    featuring the result
    :param fq: samples names
    :return: 0
    """
    print(f'###### start {fq} featuring ######')
    try:
        if args.zfeature:
            otu_table = pd.read_table(f'{fq}_zotu.xls', sep='\t')
        if args.feature:
            otu_table = pd.read_table(f'{fq}_otu.xls', sep='\t')
        otu_table.rename(columns={'#OTU_ID': 'first'}, inplace=True)
        otu_table.reset_index(inplace=True)
        otu_table.rename(columns={'index': '#OTU_ID'}, inplace=True)
        if len(otu_table.columns) > 2:
            otu_table['sum'] = otu_table.iloc[:, 1:].sum(1)
            otu_table = otu_table[['#OTU_ID', 'sum']]
            otu_table.rename(columns={'sum': f'{fq}'}, inplace=True)
        otu_aln = pd.read_table(f'{fq}_aln.xls', sep='\t')
        otu_aln = otu_aln[otu_aln['evalue'] < 1e-10]
        otu_aln = otu_aln[otu_aln['alignment length'] > 150]
        otu_aln = otu_aln.iloc[:, : 3]
        otu_aln.rename(columns={'query id': '#OTU_ID', 'subject id': 'taxonomy'}, inplace=True)
        otu_aln.sort_values(by='% identity', ascending=False, inplace=True)
        otu_aln.drop_duplicates(subset='#OTU_ID', keep='first', ignore_index=True, inplace=True)
        otu_t_a = pd.merge(otu_table, otu_aln, on='#OTU_ID')
        otu_t_a.insert(3, 'Taxonomy', np.ones(len(otu_t_a.index)))
        otu_t_a.insert(4, 'Order', np.ones(len(otu_t_a.index)))
        otu_t_a.insert(5, 'Family', np.ones(len(otu_t_a.index)))
        for i in otu_t_a.index:
            otu_t_a.iat[i, 3] = tax[tax['taxonomy'].str.contains(otu_t_a.iat[i, 2].split('|')[1])].iat[0, 1]
            otu_t_a.iat[i, 4] = tax[tax['taxonomy'].str.contains(otu_t_a.iat[i, 2].split('|')[1])].iat[0, 2]
            otu_t_a.iat[i, 5] = tax[tax['taxonomy'].str.contains(otu_t_a.iat[i, 2].split('|')[1])].iat[0, 3]
        otu_t_a.drop(labels='taxonomy', axis=1, inplace=True)
        otu_t_a.to_csv(f'{fq}_taxonomy.xls', sep='\t', index=False)
    except FileNotFoundError:
        print(f'-!- {fq} have no zotu/otu -!-')
    return 0


if __name__ == '__main__':
    fastp = '/mnt/disk2/Lab_Users/fangong/edna_programs/fastp/fastp'
    tagcleaner = '/mnt/disk2/Lab_Users/fangong/edna_programs/tagcleaner-standalone-0.16/tagcleaner.pl'
    flash = '/mnt/disk2/Lab_Users/fangong/edna_programs/FLASH-1.2.11/flash'
    usearch = '/mnt/disk2/Lab_Users/fangong/edna_programs/usearch11/usearch11'
    # makeblastdb -in complete_partial_mitogenomes.fa -dbtype nucl -parse_seqids -out mitofish
    mitofish = '/mnt/disk2/Lab_Users/fangong/edna_programs/mitofish_v3.87/mitofish_v3_87'
    submit = os.getcwd()
    primer = pd.read_table(args.primer, sep='\t', header=0, comment='#', converters={'primer': str})
    samplefq = [j for i in primer['sample'] for j in i.split(',')]
    tax = pd.read_table('/mnt/disk2/Lab_Users/fangong/edna_programs/mitofish/taxonomy.txt')

    if args.raw:
        os.popen('mkdir unzipped')
        os.chdir(args.raw)
        gz = glob.glob('*.gz')

        gunzip = Pool(20)
        gunzip.map(unzip, gz)
        gunzip.close()
        gunzip.join()

        os.chdir(submit)

    if args.qc:
        os.popen('mkdir cleaned')
        os.chdir('unzipped')
        fqs = [i.strip() for i in os.popen("ls *.fq | sed 's/.R..fq//' | uniq").readlines()]

        fp = Pool(10)
        fp.map(qc, fqs)
        fp.close()
        fp.join()

        os.chdir(submit)
        os.system('mkdir -p intermediate/fastp/html_jason')
        os.system('mv unzipped/*.html unzipped/*.json intermediate/fastp/html_jason')

    if args.assemble:
        os.system('mkdir merged')
        os.chdir('merged')
        assembly = Pool(10)
        assembly.map(assembling, samplefq)
        assembly.close()
        assembly.join()
        os.chdir(submit)

        os.system('mkdir intermediate/flash')
        os.system('mv merged/*.fastq merged/*.his* intermediate/flash')

    if args.trim:
        os.system('mkdir trimmed')
        trimming()

    if args.mapping:
        if args.ref == 0:
            raise Exception('Reference did not specified')
        else:
            os.chdir('trimmed')
            ref = os.path.basename(args.ref)
            os.system(f'ln -s {args.ref} {ref}')
            os.system(f'bwa index {ref}')
            mapp = Pool(10)
            mapp.map(mapping, samplefq)
            mapp.close()
            mapp.join()
            if os.path.exists('mapping_result'):
                os.remove('mapping_result')
            os.mkdir('mapping_result')
            os.system('mv *.txt *.sam *bam *fasta *fasta.* *.bai mapping_result/')
            os.system(f'mv mapping_result/ {submit}')
            os.chdir(submit)

    if args.zotu or args.otu:
        if args.zotu:
            os.system('mkdir zotu')
            os.chdir('zotu')
        if args.otu:
            os.system('mkdir otu')
            os.chdir('otu')

        otu = Pool(10)
        otu.map(usearching, samplefq)
        otu.close()
        otu.join()

        os.chdir(submit)
        os.system('mkdir intermediate/usearch')
        if args.zotu:
            os.system('mv zotu/*filtered.fa zotu/*uniques.fa intermediate/usearch')
        if args.otu:
            os.system('mv otu/*filtered.fa otu/*uniques.fa intermediate/usearch')

    if args.zblast or args.blast:
        if args.zblast:
            os.chdir('zotu')
        if args.blast:
            os.chdir('otu')

        bl = Pool(10)
        bl.map(blastn, samplefq)
        bl.close()
        bl.join()

        os.chdir(submit)

    if args.zfeature or args.feature:
        if args.zfeature:
            os.chdir('zotu')
        if args.feature:
            os.chdir('otu')

        feat = Pool(10)
        feat.map(featuring, samplefq)
        feat.close()
        feat.join()

        os.system('mkdir feature')
        os.system('mv *taxonomy.xls feature/')
        os.system(f'mv feature/ {submit}')
        os.chdir(submit)
