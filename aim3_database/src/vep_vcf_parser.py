"""VEP VCF Parser

Parses a VCF file producd by VEP and inserts data into the database.
"""

import sys
import os
import vcf
import sqlite3


vcf_file = sys.argv[1]

def parse_vcf(vcf_file):
    # Read in VCF file and get isolate ID
    reader = vcf.Reader(filename=vcf_file)
    isolate = reader.filename.split('/')[-1].split('.')[0]

    # Connect to the database
    home = os.environ['GROUPHOME']
    conn = sqlite3.connect(home + '/data/depot/database/database')
    cur = conn.cursor()
    execute = cur.execute

    # Prepare sqlite commands
    # variants table
    var_ins = "INSERT OR IGNORE INTO variant(position, ref, alt) VALUES"
    var_sel = "SELECT id FROM variant WHERE (position, ref, alt) = "
    # isolates table
    iso_cmd =  "INSERT OR IGNORE INTO iso_var(iso, var) VALUES"
    # effects table
    effect_cmd = "INSERT OR IGNORE INTO effect(var, gene, locus_tag, consequence, codon, aa, distance, strand) VALUES"

    var_ids = []
    lstr = str
    # Loop over records
    for record in reader:
        info = record.INFO['CSQ']
        pos = record.POS
        ref = record.REF
        alt = lstr(record.ALT[0])
        # Insert or ignore into variants table
        cur.execute(var_ins + lstr((pos, ref, alt)))

        # Now use the select query to see if this variant's pos/ref/alt is in the table
        cur.execute(var_sel + lstr((pos, ref, alt)))
        var_id = cur.fetchone()
        var_ids.append(var_id)

        # Some variants don't have CSQ data in the INFO column
        if len(info) <= 1:
            continue
        else:
            # Pull out all the nested INFO rows
            info_list = [i.split('|') for i in info]
            # Pull out the releveant information from each INFO row
            info_row = [(var_id[0],l[3],l[4],l[1],l[14],l[15],l[18],l[19]) for l in info_list]
            # Alter the command to insert or ignore into the effects table
            effect_cmds = [effect_cmd + lstr(i) for i in info_row]
            # Insert or ignore into effects table
            result = list(map(execute, effect_cmds))

    # Insert isolate and var_id into isolates table
    values = [(isolate, i[0]) for i in var_ids]
    iso_cmds = [iso_cmd + lstr(i) for i in values]
    result = list(map(execute, iso_cmds))

    # Save changes and close connection
    conn.commit()
    conn.close()


parse_vcf(vcf_file)

