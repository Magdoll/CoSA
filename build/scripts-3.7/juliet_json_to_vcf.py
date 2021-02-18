#!/home/UNIXHOME/etseng/anacondaPy37/envs/anaCogentPy37/bin/python

__VCF_EXAMPLE__ = \
"""
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
20      1       .       G       A,T     .       PASS    AF=0.5       GT
"""
import pdb
import os, sys
import json
import vcf



def juliet_json_to_vcf(json_filename, vcf_filename, gene_pos_info, ref_name='NC_045512v2', sample_name='UnknownSample'):

    with open('template.vcf', 'w') as f:
        f.write(__VCF_EXAMPLE__ + '\n')
    reader = vcf.VCFReader(open('template.vcf'))
    reader.samples = [sample_name]

    f_vcf = vcf.Writer(open(vcf_filename, 'w'), reader)

    h = open(json_filename)
    data = json.load(h)

    for g in data['genes']:
        for v in g['variant_positions']:
            cov = v['coverage']
            ref_codon = v['ref_codon']
            abs_pos = gene_pos_info[g['name']] + 3 * (v['ref_position'] - 1)
            #flag_is_primer = has_pos_overlap(abs_pos, primer_regions)
            ind = 0
            for vac in v['variant_amino_acids']:
                for cur in vac['variant_codons']:
                    ind += 1
                    _id = "{g}.{r}.{ind}".format(g=g['name'], r=v['ref_position'], ind=ind)
                    codon = cur['codon']
                    codon_offset = 0

                    for codon_offset in range(3):
                        if ref_codon[codon_offset] != codon[codon_offset]:
                            freq = "{0:.6f}".format(cur['frequency'])
                            #pdb.set_trace()
                            rec = vcf.model._Record(CHROM=ref_name,
                                                    POS=abs_pos+codon_offset,
                                                    ID=_id,
                                                    REF=ref_codon[codon_offset],
                                                    ALT=[vcf.model._Substitution(codon[codon_offset])],
                                                    QUAL='.', FILTER='PASS',
                                                    INFO={'AF': freq, 'DP': cov},
                                                    FORMAT="GT",
                                                    sample_indexes=None)
                            samp_ft = vcf.model.make_calldata_tuple(['GT'])
                            rec.samples.append(vcf.model._Call(rec, sample_name, samp_ft(*["0|1"])))
                            f_vcf.write_record(rec)
    f_vcf.close()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("juliet_json", help=".json output from running juliet")
    parser.add_argument("config_json", help=".json config file used in running juliet")
    parser.add_argument("output_vcf", help="output vcf filename")
    #parser.add_argument("primer_bed", help="primer BED file")

    args = parser.parse_args()

    gene_pos_info = {}
    h = open(args.config_json)
    d = json.load(h)
    for g in d['genes']: gene_pos_info[g['name']] = g['begin']

    #primer_regions = read_bed(args.primer_bed)

    juliet_json_to_vcf(args.juliet_json, args.output_vcf, gene_pos_info)