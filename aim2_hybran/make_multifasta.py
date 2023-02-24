'''
author: Afif Elghraoui
modified by: Emma Bishop
'''

import tempfile
from Bio import SeqIO


def main():
    proj_dir = '/home/ebishop/thesis_projects/aim2_hybran/homology/'
    fofn = proj_dir + 'gbk_fofn'  
    orfs_f = proj_dir + 'orfs_list'
    out_fasta = proj_dir + 'novel_orfs.fasta'

    # Load GenBank paths
    with open(fofn, 'r') as f:
        gbks = [line.strip() for line in f]

    # Load non-Rv/Ra ORFs
    with open(orfs_f, 'r') as f:
        orfs = [line.strip() for line in f]

    # Loop over the GBKs, getting fasta records for all the novel genes
    all_seqs = []
    for g in gbks:
        found_seqs, found_orfs = get_seq(g, orfs)
        # Update the list of orfs to search, excluding ones we found already
        orfs = [i for i in orfs if i not in found_orfs]
        # Update list of FASTA records
        all_seqs = all_seqs + found_seqs

    # Write the records to a multifasta
    with open(out_fasta, 'a') as f:
        for seq in all_seqs:
            print(SeqIO.FastaIO.as_fasta(seq), file=f)


def id(feature):
    id = feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers.keys():
        id += '|' + feature.qualifiers['gene'][0]
    return id


def extractn(feature, record):
    seq = feature.extract(record)
    seq.id = id(feature)
    return seq


def get_seq(in_gbk, in_orfs):
    '''
    Takes a GBK and list of gene names and searches for those genes in the GBK.
    Returns a list of FASTA records and a list of gene names for the found
    genes.

    Parameters:
    in_gbk (str): GBK filepath
    in_orfs (list): gene names to search
    
    Returns:
    query_list, orfs_list (tuple of lists of str): FASTA records and gene names
    '''
    # Parse GBK to get gene info
    gbk = list(SeqIO.parse(in_gbk, "genbank"))[0]
    gene_seq = dict()
    for rec in gbk.features:
        # Some sanity checking
        if rec.type not in ['CDS','rRNA','tRNA']:
            continue
        if False and 'translation' not in rec.qualifiers.keys():
            continue
        # Add to query_seq    
        if 'gene' in rec.qualifiers.keys():
            gene_seq[rec.qualifiers['gene'][0]] = rec
    # Get FASTA records
    query_list = []
    orfs_list = []
    for orf in in_orfs:
        try:
            query = extractn(gene_seq[orf], gbk)
        except KeyError:
            continue
        query_list.append(query)
        orfs_list.append(orf)
    # Return the FASTA records and gene names
    return query_list, orfs_list


if __name__ == "__main__":
    main()

