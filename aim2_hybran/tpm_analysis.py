import os
import pandas as pd

tpm_file = '/home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/data_intermed/6_sorted_genes.out'
gtf_file = '/home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/data_intermed/5_H37Rv_chr1.gtf'
outdir = '/home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/results/'

# Load TPMCalculator output
tpm_raw = pd.read_csv(tpm_file, sep='\t', usecols=['Gene_Id', 'TPM'])
print('TPM scores BEFORE removing rRNA:\n', tpm_raw.describe())

# Remove rRNA: 16S (H37Rv-NCBI_01403), 23S (H37Rv-NCBI_01403), 5S (H37Rv-NCBI_01403)
rrna = ['H37Rv-NCBI_01403', 'H37Rv-NCBI_01404', 'H37Rv-NCBI_01405']
tpm = tpm_raw.loc[~tpm_raw['Gene_Id'].isin(rrna)]
print('TPM scores AFTER removing rRNA:\n', tpm.describe())

with open(outdir + 'tpm_stats.txt', 'w') as f:
    print('TPM scores BEFORE removing rRNA:\n', tpm_raw.describe(), file=f)
    print('\nTPM scores AFTER removing rRNA:\n', tpm.describe(), file=f)

# Get novel gene info (need locus tag to match with TPMCalc output)
novel = []
with open(gtf_file, 'r') as g:
    for line in g:
        if line.startswith('#'):
            continue
        column = line.rstrip('').split('\t')
        if len(column) >= 8 and column[2] == 'CDS':
            start = column[3]
            end = column[4]
            info = column[8].split(';')
            locus_tag = ''.join([i.split(' ')[2] for i in info if i.startswith(' locus_tag')]).strip('"')
            translation = ''.join([i.split(' ')[2] for i in info if i.startswith(' translation')]).strip('"')
            gene = ''.join([i.split(' ')[2] for i in info if i.startswith(' gene')]).strip('"')
            if gene.startswith('ORF') and 'Rv' not in gene and 's' not in gene and len(translation) > 0:
                novel.append([locus_tag, gene, start, end, translation])

# Make novel genes dataframe
novel_df = pd.DataFrame(novel, columns=['Gene_Id', 'gene', 'start', 'end', 'translation'])

# Get novel ORFs with TPM >= third quartile
third = tpm.TPM.quantile(0.75)
novel_expressed = tpm.loc[(tpm['Gene_Id'].isin(novel_df['Gene_Id'])) & (tpm['TPM'] >= third)]
novel_tags = list(novel_expressed['Gene_Id'])
print('{} novel genes are being expressed (TPM score >= {})'.format(novel_expressed.shape[0], third))

# Include other GTF info and write to file
out_df = pd.merge(novel_expressed, novel_df, on='Gene_Id')
out_df.to_csv(outdir + 'expressed_novel_proteins.tsv', sep='\t', index=False)
print('Expressed genes written to', outdir + 'expressed_novel_proteins.tsv')

# Write fasta files
fasta_outdir = outdir + 'fastas/'
if not os.path.exists(fasta_outdir):
    os.mkdir(fasta_outdir)

for i in out_df.index:
    orf =  out_df['gene'][i]
    aa = out_df['translation'][i]
    with open(fasta_outdir + orf + '.fa', 'w') as f:
        print('>' + orf + '\n' + aa, file=f)
print('Expressed gene fastas written to', fasta_outdir)
