import pandas as pd

# Load data
proj_dir = '/home/ebishop/thesis_projects/aim2_hybran/homology/'
x_in = proj_dir + 'blastx.out'
rx_in = proj_dir + 'ra_blastx.out'
n_in = proj_dir + 'blastn.out'
orfs_in = proj_dir + 'orfs_list'

all_orfs = pd.read_csv(orfs_in, sep=' ', header=None)
all_orfs.rename({0:'orf'}, axis=1, inplace=True)
blastn = pd.read_csv(n_in, sep='\t', header=None) 
blstx_notra = pd.read_csv(x_in, sep='\t', header=None)
blstx_ra = pd.read_csv(rx_in, sep='\t', header=None)
blastx = pd.concat([blstx_notra, blstx_ra])


def analyze_blast(df, blast_type, all_orfs=all_orfs):
    '''Count up each kind of result (hit, incomplete, no hit)'''
    # Rename and split columns
    df.rename({0:'qry', 1:'hit', 2:'evalue', 3:'bitscore'}, axis=1, inplace=True)
    df[['isolate', 'orf']] = df['qry'].str.split('|', expand=True)
    df.drop(columns='qry', inplace=True)
    # Get duplicated orfs (had incomplete alignment)
    mask = df.orf.duplicated(keep=False)
    dup = df[mask][['orf']].drop_duplicates(ignore_index=True)
    dup[blast_type] = 'incomplete'
    # Get missing orfs (had no hits and don't appear in blast output)
    missing = all_orfs.copy()
    missing = missing[~missing.orf.isin(df.orf)]
    missing[blast_type] = 'no_hit'
    # Get orfs with blast hits with complete alignments
    exclude = pd.concat([dup.orf, missing.orf])
    for_plt = df[~df.orf.isin(exclude)].reset_index()
    # Create df with result for each orf
    hit_orfs = for_plt.copy()
    hit_orfs = hit_orfs[['orf']]
    hit_orfs[blast_type] = 'hit'
    orfs_dat = pd.concat([hit_orfs, dup, missing]).sort_values('orf', ignore_index=True)

    return for_plt, orfs_dat


x_to_plot, x_dat = analyze_blast(blastx, 'blastx')
n_to_plot, n_dat = analyze_blast(blastn, 'blastn')
orf_dat = x_dat.merge(n_dat)
# Save
orf_dat.to_csv(proj_dir + 'orf_blast_status.csv', index=False)

# Print total hit, incomplete, no_hit counts for blastx/n
print('Number of each BLASTX/N result')
print(orf_dat.value_counts('blastx').sort_index(), '\n\n', orf_dat.value_counts('blastn').sort_index())

# Make dataframe of BLASTX homology results
known = len(x_to_plot[x_to_plot['bitscore'] >= 90])
homol = len(x_to_plot[x_to_plot['bitscore'].between(40, 90)])
unk = len(x_to_plot[x_to_plot['bitscore'] <= 40])
# Total ORFs - ORFs with complete, aligned BLASTX hits
no_hit_incomplete = orf_dat.shape[0] - x_to_plot.shape[0]

homol = pd.DataFrame({'result_type': ['known', 'homologous', 'unknown', 'no_hit_incomplete'],
                     'count': [known, homol, unk, no_hit_incomplete]})
# Save
homol.to_csv(proj_dir + 'blastx_homology_status.csv', index=False)
