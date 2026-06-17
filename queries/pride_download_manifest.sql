-- Files to download for curated PRIDE picks (result tables only, no RAW).
WITH picks AS (
  SELECT unnest([
    'PXD035902', 'PXD039536', 'PXD058744', 'PXD064117', 'PXD064782',
    'PXD033062', 'PXD042838', 'PXD033043', 'PXD036527', 'PXD063995',
    'PXD033026', 'PXD014785'
  ]) AS accession
),
files AS (
  SELECT
    f.project_accession,
    f.file_name,
    f.file_size_bytes,
    f.primary_public_location_value,
    CASE
      WHEN lower(f.file_name) IN ('psm.tsv', 'combined_modified_peptide.tsv') THEN 'fragpipe'
      WHEN lower(f.file_name) LIKE '%report.parquet' OR lower(f.file_name) = 'report.tsv' THEN 'diann'
      WHEN lower(f.file_name) LIKE '%diann%' AND lower(f.file_name) LIKE '%.zip' THEN 'diann_zip'
      WHEN lower(f.file_name) LIKE '%hexnac%site%' AND lower(f.file_name) LIKE '%.txt' THEN 'maxquant_sites'
      WHEN lower(f.file_name) LIKE '%.zip'
           AND (lower(f.file_name) LIKE '%hexnac%' OR lower(f.file_name) LIKE '%oglcnac%'
                OR lower(f.file_name) LIKE '%o-glcnac%' OR lower(f.file_name) LIKE '%og_site%'
                OR lower(f.file_name) LIKE '%st_site%') THEN 'site_zip'
      WHEN (lower(f.file_name) LIKE '%o-glcnac%site%' OR lower(f.file_name) LIKE '%oglcnac%site%'
            OR lower(f.file_name) LIKE '%og_site%') AND lower(f.file_name) LIKE '%.txt' THEN 'site_table'
      WHEN lower(f.file_name) LIKE '%.mztab' THEN 'mztab_result'
      WHEN lower(f.file_name) LIKE '%psm%' AND lower(f.file_name) LIKE '%.txt' THEN 'pd_psm'
      WHEN lower(f.file_name) LIKE '%reporterions%' AND lower(f.file_name) LIKE '%oglcnac%' THEN 'skyline_reporter'
      WHEN lower(f.file_name) LIKE '%.csv' AND (lower(f.file_name) LIKE '%oglcnac%' OR lower(f.file_name) LIKE '%hexnac%') THEN 'csv_table'
      WHEN lower(f.file_name) LIKE '%.xlsx' AND lower(f.file_name) LIKE '%glcnac%' THEN 'xlsx_table'
      ELSE NULL
    END AS engine_kind
  FROM read_parquet('data/pride/silver/project_files/snapshot_date=*/*.parquet', hive_partitioning := true) f
  INNER JOIN picks p ON f.project_accession = p.accession
  WHERE lower(f.file_name) NOT LIKE '%.raw'
    AND lower(f.file_name) NOT LIKE '%.mgf'
    AND lower(f.file_name) NOT LIKE '%.mzml'
    AND lower(f.file_name) NOT LIKE '%.wiff'
    AND f.primary_public_location_value IS NOT NULL
)
SELECT project_accession, file_name, file_size_bytes, primary_public_location_value, engine_kind
FROM files
WHERE engine_kind IS NOT NULL
ORDER BY project_accession, file_name;
