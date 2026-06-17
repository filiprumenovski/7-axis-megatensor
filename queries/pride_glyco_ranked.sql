-- Rank glyco PRIDE projects with deposited site/PSM tables for adapter feasibility.
WITH glyco AS (
  SELECT
    accession,
    title,
    countries,
    list_transform(instruments, x -> x.name) AS instrument_names,
    list_transform(softwares, x -> x.name) AS software_names
  FROM read_parquet('data/pride/bronze/projects/snapshot_date=*/**/*.parquet', hive_partitioning := true)
  WHERE list_contains(list_transform(keywords, k -> lower(k)), 'o-glcnac')
     OR lower(title) LIKE '%glcnac%'
     OR lower(title) LIKE '%o-glcnac%'
     OR EXISTS (
       SELECT 1 FROM UNNEST(identifiedPTMStrings) t(ptm)
       WHERE lower(ptm.name) LIKE '%glcnac%' OR lower(ptm.name) LIKE '%hexnac%'
     )
),
files AS (
  SELECT
    project_accession,
    file_name,
    primary_public_location_value,
    CASE
      WHEN lower(file_name) LIKE '%hexnac%site%' AND lower(file_name) LIKE '%.txt' THEN 'maxquant_sites'
      WHEN lower(file_name) LIKE '%glcnac%site%' AND lower(file_name) LIKE '%.txt' THEN 'site_table'
      WHEN lower(file_name) LIKE '%hexnac%site%' AND lower(file_name) LIKE '%.zip' THEN 'site_zip'
      WHEN lower(file_name) LIKE '%psm%' AND lower(file_name) LIKE '%.txt' THEN 'pd_psm'
      WHEN lower(file_name) LIKE '%report.tsv' OR lower(file_name) LIKE '%report.parquet' THEN 'diann_report'
      WHEN lower(file_name) IN ('psm.tsv', 'combined_modified_peptide.tsv') THEN 'fragpipe'
      ELSE NULL
    END AS result_kind
  FROM read_parquet('data/pride/silver/project_files/snapshot_date=*/*.parquet', hive_partitioning := true)
),
scored AS (
  SELECT
    g.accession,
    g.title,
    g.countries,
    g.instrument_names,
    g.software_names,
    list(DISTINCT f.result_kind) FILTER (WHERE f.result_kind IS NOT NULL) AS result_kinds,
    count(f.file_name) FILTER (WHERE f.result_kind IS NOT NULL) AS result_files,
    CASE
      WHEN list_contains(g.software_names, 'MaxQuant')
        OR array_to_string(g.software_names, ',') LIKE '%MaxQuant%' THEN 'maxquant'
      WHEN list_contains(g.software_names, 'Proteome Discoverer')
        OR array_to_string(g.software_names, ',') LIKE '%Proteome Discoverer%' THEN 'pd'
      WHEN array_to_string(g.software_names, ',') ILIKE '%fragpipe%'
        OR array_to_string(g.software_names, ',') ILIKE '%msfragger%' THEN 'fragpipe'
      WHEN array_to_string(g.software_names, ',') ILIKE '%diann%' THEN 'diann'
      WHEN array_to_string(g.software_names, ',') ILIKE '%skyline%' THEN 'skyline'
      ELSE 'other'
    END AS engine_bucket,
    (
      lower(g.title) LIKE '%treatment%'
      OR lower(g.title) LIKE '%vs%'
      OR lower(g.title) LIKE '%knock%'
      OR lower(g.title) LIKE '%overexpression%'
      OR lower(g.title) LIKE '%light%'
      OR lower(g.title) LIKE '%heavy%'
      OR lower(g.title) LIKE '%insulin%'
      OR lower(g.title) LIKE '%resistance%'
    ) AS likely_conditioned
  FROM glyco g
  LEFT JOIN files f ON g.accession = f.project_accession
  GROUP BY ALL
)
SELECT *
FROM scored
WHERE result_files > 0
ORDER BY likely_conditioned DESC, result_files DESC, accession;
