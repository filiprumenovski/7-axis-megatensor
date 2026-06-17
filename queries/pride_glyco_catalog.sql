-- File-first O-GlcNAc / HexNAc(ST) PRIDE catalog with deposited result tables.
-- Run against unpacked pride-ingest snapshot (bronze projects + silver project_files).
WITH projects AS (
  SELECT
    accession,
    title,
    countries,
    list_transform(instruments, x -> x.name) AS instrument_names,
    list_transform(softwares, x -> x.name) AS software_names,
    list_transform(labPIs, x -> x.country) AS pi_countries,
    list_contains(list_transform(keywords, k -> lower(k)), 'o-glcnac') AS kw_oglcnac,
    lower(title) LIKE '%o-glcnac%'
      OR lower(title) LIKE '%oglcnac%'
      OR (lower(title) LIKE '%glcnac%' AND lower(title) NOT LIKE '%n-glycan%') AS title_glyco,
    EXISTS (
      SELECT 1 FROM UNNEST(identifiedPTMStrings) t(ptm)
      WHERE lower(ptm.name) LIKE '%glcnac%' OR lower(ptm.name) LIKE '%hexnac%'
    ) AS ptm_glyco,
    array_to_string(list_transform(instruments, x -> x.name), ',') ILIKE '%astral%' AS has_astral,
    array_to_string(list_transform(instruments, x -> x.name), ',') ILIKE '%exploris%' AS has_exploris,
    array_to_string(list_transform(instruments, x -> x.name), ',') ILIKE '%tims%' AS has_timstof
  FROM read_parquet('data/pride/bronze/projects/snapshot_date=*/**/*.parquet', hive_partitioning := true)
),
files AS (
  SELECT
    project_accession,
    file_name,
    file_category_name,
    primary_public_location_value,
    CASE
      WHEN lower(file_name) IN ('psm.tsv', 'combined_modified_peptide.tsv') THEN 'fragpipe'
      WHEN lower(file_name) LIKE '%fragpipe%' AND lower(file_name) LIKE '%.zip' THEN 'fragpipe_zip'
      WHEN lower(file_name) LIKE '%msfragger%' AND lower(file_name) LIKE '%.zip' THEN 'fragpipe_zip'
      WHEN lower(file_name) = 'report.tsv' OR lower(file_name) LIKE '%report.parquet' THEN 'diann'
      WHEN lower(file_name) LIKE '%diann%' AND lower(file_name) LIKE '%.zip' THEN 'diann_zip'
      WHEN lower(file_name) LIKE '%hexnac%site%' AND lower(file_name) LIKE '%.txt' THEN 'maxquant_sites'
      WHEN lower(file_name) LIKE '%.zip'
           AND (lower(file_name) LIKE '%hexnac%' OR lower(file_name) LIKE '%oglcnac%'
                OR lower(file_name) LIKE '%o-glcnac%' OR lower(file_name) LIKE '%og_site%'
                OR lower(file_name) LIKE '%st_site%') THEN 'site_zip'
      WHEN (lower(file_name) LIKE '%o-glcnac%site%' OR lower(file_name) LIKE '%oglcnac%site%'
            OR lower(file_name) LIKE '%og_site%') AND lower(file_name) LIKE '%.txt' THEN 'site_table'
      WHEN lower(file_name) LIKE '%.mztab' THEN 'mztab_result'
      WHEN lower(file_name) LIKE '%psm%' AND lower(file_name) LIKE '%.txt' THEN 'pd_psm'
      WHEN lower(file_name) LIKE '%reporterions%' AND lower(file_name) LIKE '%oglcnac%' THEN 'skyline_reporter'
      WHEN lower(file_name) LIKE '%.csv' AND (lower(file_name) LIKE '%oglcnac%' OR lower(file_name) LIKE '%hexnac%') THEN 'csv_table'
      WHEN lower(file_name) LIKE '%.xlsx' AND lower(file_name) LIKE '%glcnac%' THEN 'xlsx_table'
      ELSE NULL
    END AS engine_kind,
    (lower(file_name) LIKE '%glcnac%' OR lower(file_name) LIKE '%hexnac%'
      OR lower(file_name) LIKE '%oglcnac%' OR lower(file_name) LIKE '%o-glcnac%') AS file_glyco_hint
  FROM read_parquet('data/pride/silver/project_files/snapshot_date=*/*.parquet', hive_partitioning := true)
  WHERE lower(file_name) NOT LIKE '%.raw'
    AND lower(file_name) NOT LIKE '%.mgf'
    AND lower(file_name) NOT LIKE '%.mzml'
    AND lower(file_name) NOT LIKE '%.wiff'
),
agg AS (
  SELECT
    project_accession,
    engine_kind,
    count(*) AS result_files,
    sum(file_glyco_hint::int) AS glyco_named_files,
    list_slice(list(DISTINCT file_name), 1, 6) AS sample_files,
    list_slice(list(DISTINCT primary_public_location_value), 1, 3) AS download_urls
  FROM files
  WHERE engine_kind IS NOT NULL
  GROUP BY 1, 2
)
SELECT
  p.accession,
  p.title,
  p.countries,
  p.instrument_names,
  p.software_names,
  a.engine_kind,
  a.result_files,
  a.glyco_named_files,
  a.sample_files,
  a.download_urls,
  (p.kw_oglcnac OR p.title_glyco OR p.ptm_glyco) AS metadata_glyco,
  p.has_astral,
  p.has_exploris,
  p.has_timstof,
  (
    lower(p.title) LIKE '%treatment%' OR lower(p.title) LIKE '% vs %' OR lower(p.title) LIKE '%knock%'
    OR lower(p.title) LIKE '%overexpression%' OR lower(p.title) LIKE '%resistance%'
    OR lower(p.title) LIKE '%insulin%' OR lower(p.title) LIKE '%light%' OR lower(p.title) LIKE '%heavy%'
    OR lower(p.title) LIKE '%diabetic%' OR lower(p.title) LIKE '%alzheimer%'
    OR array_to_string(a.sample_files, ',') ILIKE '%insulin%'
    OR array_to_string(a.sample_files, ',') ILIKE '%light%'
    OR array_to_string(a.sample_files, ',') ILIKE '%heavy%'
  ) AS likely_conditioned
FROM agg a
JOIN projects p ON a.project_accession = p.accession
WHERE p.kw_oglcnac OR p.title_glyco OR p.ptm_glyco OR a.glyco_named_files > 0
ORDER BY metadata_glyco DESC, likely_conditioned DESC, a.result_files DESC, p.accession;
