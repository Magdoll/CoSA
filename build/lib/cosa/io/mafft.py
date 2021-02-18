"""
Process mafft which must contain the reference "NC_045512v2"
"""
import os, sys
from Bio import SeqIO

REFNAME = 'NC_045512v2'

def get_coord_mapping_for_ref(aln):
    """
    :param aln: sequence alignnment for the reference
    :return: dict of (0-based aln pos, 0-based ref pos) and dict of (0-based ref pos, 0-based aln pos)
    """
    aln_i = 0
    ref_i = 0

    map_aln_to_ref = {}
    map_ref_to_aln = {}
    for i,s in enumerate(aln):
        if s != '-':
            map_ref_to_aln[ref_i] = aln_i
            map_aln_to_ref[aln_i] = ('M', ref_i)
            ref_i += 1
        else:
            # insertion w.r.t reference
            map_aln_to_ref[aln_i] = ('I', ref_i)
        aln_i += 1
    return map_aln_to_ref, map_ref_to_aln



