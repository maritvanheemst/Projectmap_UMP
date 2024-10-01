#!/usr/bin/env python

"""
Sequence-based structural alignment of two proteins.
"""

from __future__ import print_function, division

import argparse
import os
from argparse import Namespace

from Bio.PDB import PDBParser, FastMMCIFParser, Superimposer, PDBIO
from Bio.PDB.Polypeptide import is_aa

from Bio import pairwise2
#from Bio.SubsMat import MatrixInfo as matlist
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

def align_sequences(structA, structB, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    """

    def _calculate_identity(sequenceA, sequenceB):
        """
        Returns the percentage of identical characters between two sequences.
        Assumes the sequences are aligned.
        """

        sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
        matches = [sa[i] == sb[i] for i in range(sl)]
        seq_id = (100 * sum(matches)) / sl

        gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
        gap_id = (100 * sum(matches)) / gapless_sl
        return (seq_id, gap_id)

    def _get_pdb_sequence(structure):
        """
        Retrieves the AA sequence from a PDB structure.
        """

        _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
        seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
        return seq

    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)
    trim_ends = kwargs.get('trim_ends', True)

    resseq_A = _get_pdb_sequence(structA)
    resseq_B = _get_pdb_sequence(structB)

    sequence_A = ''.join([i[1] for i in resseq_A])
    sequence_B = ''.join([i[1] for i in resseq_B])
    alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                    matrix, gap_open, gap_extend,
                                    penalize_end_gaps=(False, False) )

    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln

    # Equivalent residue numbering
    # Relative to reference
    mapping = {}
    aa_i_A, aa_i_B = 0, 0
    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
        if aa_aln_A == '-':
            if aa_aln_B != '-':
                aa_i_B += 1
        elif aa_aln_B == '-':
            if aa_aln_A != '-':
                aa_i_A += 1
        else:
            assert resseq_A[aa_i_A][1] == aa_aln_A
            assert resseq_B[aa_i_B][1] == aa_aln_B
            mapping[resseq_A[aa_i_A][0]] = resseq_B[aa_i_B][0]
            aa_i_A += 1
            aa_i_B += 1

    # Gapless alignment
    # def _trimmer(sequence):
    #     """Returns indices of first and last ungapped position"""

    #     leading = [i for (i, aa) in enumerate(sequence) if aa != '-'][0]
    #     trailing = [i for (i, aa) in enumerate(sequence[::-1]) if aa != '-'][0]

    #     trailing = len(sequence) - trailing
    #     return (leading, trailing)

    # lead_A, trail_A = _trimmer(aligned_A)
    # lead_B, trail_B = _trimmer(aligned_B)

    # lead = max(lead_A, lead_B)
    # trail = min(trail_A, trail_B)
    # trim_aln_A = aligned_A[lead:trail]
    # trim_aln_B = aligned_B[lead:trail]
    # mismatch = ''.join(['+' if a!=b else ' ' for (a,b) in zip(trim_aln_A, trim_aln_B)])

    # Calculate (gapless) sequence identity
    seq_id, g_seq_id = _calculate_identity(aligned_A, aligned_B)
    return ((aligned_A, aligned_B), seq_id, g_seq_id, mapping)
    # return ((trim_aln_A, trim_aln_B, mismatch), seq_id, g_seq_id, mapping)

def parse_structure(spath):
    """Parses a PDB/cif structure"""

    if not os.path.isfile(spath):
        return IOError('File not found: {0}'.format(spath))

    if spath.endswith(('pdb', 'ent')):
        parser = PDBParser()
    elif spath.endswith('cif'):
        parser = FastMMCIFParser()
    else:
        raise Exception('Format not supported ({0}). Must be .pdb/.ent or .cif'.format(spath))

    sname = os.path.basename(spath.split('.')[0])
    return parser.get_structure(sname, spath)
##

def run(refe, mobi, r_chain='A',m_chain='A'):

    #p= argparse.ArgumentParser(description=__doc__)
    #p.add_argument('refe', help='Reference Structure')
    #p.add_argument('--r_chain', default='A', help='Reference Structure Chain')
    #p.add_argument('mobi', help='Mobile Structure')
    #p.add_argument('--m_chain', default='A', help='Reference Structure Chain')
    #line = ap.parse_args()
    
    args= {'refe':refe,'mobi':mobi,'r_chain':r_chain,'m_chain':m_chain}
    
    cline = Namespace(**args)

    # Parse structures & take only the necessary chain
    s_reference = parse_structure(cline.refe)
    try:
        reference = s_reference[0][cline.r_chain]
    except KeyError:
        raise Exception('Chain {0} not found in reference structure'.format(cline.r_chain))

    s_mobile = parse_structure(cline.mobi)
    try:
        mobile = s_mobile[0][cline.m_chain]
    except KeyError:
        raise Exception('Chain {0} not found in mobile structure'.format(cline.m_chain))

    # Align sequences to get mapping between residues
    aln, seq_id, gapless_id, res_map = align_sequences(reference, mobile)

    refe_ca_list, mobi_ca_list = [], []
    for refe_res in res_map:
        refe_ca_list.append(reference[refe_res]['CA'])
        mobi_ca_list.append(mobile[res_map[refe_res]]['CA'])

    # Superimpose matching residues
    si = Superimposer()
    si.set_atoms(refe_ca_list, mobi_ca_list)

    # Transform & Write Mobile
    si.apply(mobile.get_atoms())

    io = PDBIO()
    io.set_structure(mobile)
    m_transformed_name = '{0}_transformed.pdb'.format(s_mobile.id)
    io.save(m_transformed_name)