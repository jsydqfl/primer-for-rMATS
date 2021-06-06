#!/usr/bin/env python
import os
import sys
import io
import primer3
from optparse import OptionParser as OP
import pandas as pd


# 前期准备
# 1. samtools放在环境变量里
# 2. 安装primer3：c

# options
def cmdparameter(argv):
    if len(argv) <= 0:
        global desc
        print(desc)
        cmd = 'python' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -f .txt -g .fa -o ./"
    parser = OP(usage=usages)
    parser.add_option("-f", "--filename", dest="filename",
        default="SE.txt", help="SE file name")
    parser.add_option("-g", "--genome_filename", dest="genome_filename",
        default="All-GRCh38.fa", help="Genome fasta file name")
    parser.add_option("-o", "--output", dest="output",
        default="./", help="Output prefix")
    (options, args) = parser.parse_args(argv[1:])
    return (options, args)
#------------END options-----------------------------------------

# primerdesign
def primer3Code(seq_id, seq_seq, seq_start, seq_end):
    # Primer3 全局参数
    # 前提条件 PRIMER_MAX_SIZE<= min PRIMER_PRODUCT_SIZE_RANGE; len(PRIMER_PRODUCT_SIZE_RANGE) >=2
    sl = len(seq_seq)
    global_args = {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': [[sl,sl+25],[sl+25,sl+50],[sl+50,sl+75],[sl+75,sl+100]],
    }
        # design primer
    seq_args = {
            'SEQUENCE_ID': seq_id,
            'SEQUENCE_TEMPLATE': seq_seq,}
    primer3_result = primer3.bindings.designPrimers(seq_args, global_args)

    # result format
    primer3_result_table_dict = {}
    for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
        primer_id = str(i)
        for key in primer3_result:
            if primer_id in key:
                info_tag = key.replace("_" + primer_id, "")

                try:
                    primer3_result_table_dict[info_tag]
                except:
                    primer3_result_table_dict[info_tag] = []
                finally:
                    primer3_result_table_dict[info_tag].append(primer3_result[key])

    index = []
    for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):

        index.append("PRIMER_PAIR_" + str(i))

    primer3_result_df = pd.DataFrame(primer3_result_table_dict, index=index)

    primer3_result_df = primer3_result_df.T
    return(primer3_result_df)
#------------END Primer Code------

if __name__ == '__main__':
    options, args = cmdparameter(sys.argv)
    n = 0
    with open(options.filename) as fn:
        for line in fn:
            lineL = line.strip("\n").split("\t")

            if n == 0:
                ID_loc = lineL.index('ID')
                GeneID_loc = lineL.index('GeneID')
                chr_loc = lineL.index('chr')
                exonStart_0base_loc = lineL.index('exonStart_0base')
                exonEnd_loc = lineL.index('exonEnd')
                upstreamEE_loc = lineL.index('upstreamEE')
                downstreamES_loc = lineL.index('downstreamES')
                n += 1
            else:
                seq_id = lineL[GeneID_loc] + "_" + lineL[ID_loc]
                print(seq_id)
                exon_20 = lineL[chr_loc] + ":" + str(int(lineL[exonEnd_loc])-150) + "-" + str(int(lineL[downstreamES_loc])+150)
                seq_seq = "".join(os.popen("samtools faidx " + options.genome_filename + " " + exon_20).read().split("\n")[1:])

                if len(seq_seq) != 0:
                    print(exon_20, seq_seq)
                    primer3_result_df1 = primer3Code(seq_id + "_exon", seq_seq, 0, len(seq_seq))
                else:
                    print("%s is zero length." % exon_20)

                upEdownS_20 = lineL[chr_loc] + ":" + str(int(lineL[upstreamEE_loc])-150) + "-" + str(int(lineL[downstreamES_loc])+150)
                seq_seq = "".join(os.popen("samtools faidx " + options.genome_filename + " " + upEdownS_20).read().split("\n")[1:])

                if len(seq_seq) != 0:
                    print(upEdownS_20, seq_seq)
                    primer3_result_df = primer3Code(seq_id + "_ES", seq_seq, 0, len(seq_seq))
                else:
                    print("%s is zero length." % exon_20)

                primer3_result_df2 = pd.concat([primer3_result_df1,primer3_result_df],axis=0)
                if primer3_result_df2.shape[0] != 0:
                    primer3_result_df2.to_csv(seq_id.replace('\"','') + "_primer3_result.csv")