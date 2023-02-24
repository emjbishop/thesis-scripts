# Database

A SQLite database to store DST and variant data for _M. tuberculosis_ clinical isolates.

### Dependencies

Everything you need is already installed on the lab's computers or base conda environment. You only need to activate the conda environment to import variant data.

## Building the database

These steps describe how the database was actually built (somewhat inefficiently).

### 1. Create an empty database

Create the database file (if it doesn't exist) and launch the SQLite console:
```
sqlite3 $GROUPHOME/data/depot/database/database
```

Paste this schema at the prompt and hit `Enter` to build the tables:
```
CREATE TABLE isolate (
    id TEXT NOT NULL UNIQUE,
    bioproject TEXT
);
CREATE TABLE dst (
    isolate TEXT NOT NULL UNIQUE,
    inh TEXT,
    rif TEXT,
    amk TEXT,
    cap TEXT,
    kan TEXT,
    mox TEXT,
    ofx TEXT,
    pza TEXT,
    FOREIGN KEY (isolate)
    REFERENCES "isolate" (id)
        ON UPDATE CASCADE
        ON DELETE CASCADE
    UNIQUE (isolate, inh, rif, amk, cap, kan, mox, ofx, pza)
);
CREATE TABLE effect (
    var INTEGER,
    gene TEXT,
    locus_tag TEXT,
    consequence TEXT,
    codon TEXT,
    aa TEXT,
    distance TEXT,
    strand TEXT,
    UNIQUE (var, gene, locus_tag, consequence, codon, aa, distance, strand)
);
CREATE TABLE variant (
    id INTEGER PRIMARY KEY,
    position TEXT NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    UNIQUE (position, ref, alt)
);
CREATE TABLE iso_var (
    iso TEXT NOT NULL,
    var INT NOT NULL,
    UNIQUE (iso, var)
);
CREATE INDEX idx_pos_ref_alt ON variant (position, ref, alt);
```

### 2. Import data

##### For the `isolates` table:
1. Manually extract all the isolate IDs we care about:
     - All the IDs from `$GROUPHOME/metadata/gcdd.lineage.txt`
     - Just the SEA* IDs from `$GROUPHOME/metadata/dst.txt`
2. Store these in `$GROUPHOME/data/depot/database/intermediate_files/isolates.csv`
3. Import and update the `bioproject` column (these are all GCDD isolates):
    ```
    sqlite3 $GROUPHOME/data/depot/database/database
    .import $GROUPHOME/data/depot/database/intermediate_files/isolates.csv isolate
    UPDATE isolates SET bioproject = 'GCDD';
    ```

##### For the `dst` table:
```
# From the top level of the repository:
./src/isolate_dst.sh
sqlite3 $GROUPHOME/data/depot/database/database
.import $GROUPHOME/data/depot/database/intermediate_files/isolate_dst.csv dst
```

##### For the `variant`, `iso_var`, and `effect` tables:
```
# From the top level of the repository:
conda activate
./src/import_variants.sh
```

## Example queries

You can paste these directly into the console or incorporate them into a Python script.

##### Tips:
- Exit console with `.exit`
- You can paste multiple lines at once (meaning you can paste an entire query from here to the prompt)
- ";" indicates the end of a command
- You can use bash hotkeys in the SQLite console (like `ctrl+c` to kill a running query)
- SQLite commands are not case sensitive but keywords are typically written in CAPS
- All the data is type TEXT except for the `var` and `variant.id` columns
- You can output results to a file instead of the default `stdout`:
  ```
  # Write comma separated output to a .csv file
  .mode csv
  .output /path/to/output/file.csv
  ```
  Just be sure to do `.output stdout` when you're done so your file isn't
  overwritten by the next query.

##### Examples:

Get the isolate IDs and inh dst results for isolates that are inhR:
```
SELECT isolate, inh FROM dst WHERE inh = 'R';
```

Get isolate and variant information for inhR isolates, looking only at
variants affecting katG and/or inhA and that are within 500bp of those genes:
```
SELECT isolate, position, ref, alt, codon, aa, consequence, distance, gene
FROM dst
INNER JOIN iso_var ON iso_var.iso = dst.isolate
INNER JOIN variant ON iso_var.var = variant.id
INNER JOIN effect ON effect.var = variant.id
WHERE inh = 'R' AND (gene = 'katG' OR gene ='inhA') AND CAST(distance as INTEGER) <= '500'
ORDER BY isolate;
```

For variants that appear in amkR isolates, get the variant ID and how many amkR isolates have that variant:
```
SELECT var, COUNT(iso) as isoCount
FROM iso_var
INNER JOIN dst ON dst.isolate = iso_var.iso
WHERE amk = 'R'
GROUP BY var;
```

Count how many results are returned from the sub (inner) query:
```
SELECT COUNT(*) as amk_variants from
(
SELECT id
FROM variant
INNER JOIN iso_var ON variant.id = iso_var.var
INNER JOIN dst ON dst.isolate = iso_var.iso
WHERE amk = 'R'
) subquery;
```

## Further reading

The schema, scripts, and queries here were primarily informed by Microsoft's
["Database Design Basics" page](https://support.microsoft.com/en-us/office/database-design-basics-eb2159cf-1e30-401a-8084-bd4f9c9ca1f5#bmtablerelationships)
and [SQLite Tutorial](https://www.sqlitetutorial.net/).
