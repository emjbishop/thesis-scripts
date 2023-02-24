#!/usr/bin/python3.7
import gzip
import argparse


def main():
    parser = argparse.ArgumentParser(description='Converts merged vcf to tsv')
    parser.add_argument('-i', '--infile', help='Name of gz zipped merged vcf file', default='merged.vcf.gz')
    parser.add_argument('-o', '--outfile', help='Name of output file', default='merged-vcf.csv')
    args = parser.parse_args()

    inlines = []
    with gzip.open(args.infile, 'r') as f:
        for line in f:
            inlines.append(line.decode('utf-8'))

    newlines = []
    for linestring in inlines:
        if linestring.startswith('#CHROM'):  # header line, includes sample ids after vcf header
            tablist = linestring.split('\t')
            samples = tablist[9:]  # sample ids
            varinfo = ['POS', 'REF', 'ALT', 'Consequence', 'Symbol', 'Gene', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'MutID', 'GeneMutID']
            linelist = varinfo + samples
            newline = '\t'.join(linelist)
            newlines.append(newline)
        if linestring[0] != '#':  # non-header lines, includes vcf variant data and 0/1 for presence/absence in each sample
            tablist = linestring.split('\t')
            samples = tablist[9:]  # 0/1 presence or absence of variant for each sample
            pos = tablist[1]
            ref = tablist[3]
            alt = tablist[4] 
            infolist = sublist = tablist[7].split('|')  # INFO column contains rest of desired variant data. it is nested, deliminated by |
            if len(infolist) < 2:   # some variants (just one in test set) dont have all the info fields. not even empty. was VEP not run on the isolates with them?
                consequence = 'none'
                genesymbol = 'none'
                genetag = 'none'
                cpos = 'none'
                ppos = 'none'
                aminoacids = 'none'
                codons = 'none'
                mutid = ref + pos + alt
            else:
                consequence = infolist[1]
                genesymbol = infolist[3]  # Symbol  aka gene name
                genetag = infolist[4]  # Gene aka locus tag in reference
                if not genesymbol:  # some variants are missing value for SYMBOL but not Gene
                    genesymbol = genetag
                if consequence == 'upstream_gene_variant':
                   cpos = '-' + infolist[18] # Distance field in INFO column gives distance from variant to downstream CDS start
                   ppos = 'none'
                   aminoacids = 'none'
                   codons = 'none'
                   mutid = ref + cpos + alt
                else:
                    cpos = infolist[13]
                    ppos = infolist[14]
                    aminoacids = infolist[15]
                    codons = infolist[16]
                    if not aminoacids:  # if amino acids field is empty
                       mutid = ref + cpos + alt
                    elif '/' in aminoacids:  # if amino acids changed with variant, ref and alt acids on either side of /
                        aminoacidlist = aminoacids.split('/')
                        mutid = aminoacidlist[0] + ppos + aminoacidlist[1]
                    else:  # if amino acid field is not empty, but didn't change with variant. ususally a synonomous variant
                        mutid = aminoacids + ppos + aminoacids
            genemutid = genesymbol + '_' + mutid
            varinfo = [pos, ref, alt, consequence, genesymbol, genetag, cpos, ppos, aminoacids, codons, mutid, genemutid]
            linelist = varinfo + samples
            newline = '\t'.join(linelist)
            newlines.append(newline)

    with open(args.outfile, 'w') as f:
        for line in newlines:
            f.write(line)


if __name__ == '__main__':
    main()
