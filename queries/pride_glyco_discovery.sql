-- Phase 1: discover O-GlcNAc PRIDE projects (run against pride-ingest bronze)
WITH p AS (
  SELECT accession, title, countries, instruments, softwares, keywords, identifiedPTMStrings
  FROM read_parquet('data/pride/bronze/projects/snapshot_date=*/**/*.parquet', hive_partitioning := true)
)
SELECT
  accession,
  title,
  countries,
  list_transform(instruments, x -> x.name) AS instruments,
  list_transform(softwares, x -> x.name) AS softwares
FROM p
WHERE list_contains(list_transform(keywords, k -> lower(k)), 'o-glcnac')
   OR lower(title) LIKE '%glcnac%'
   OR EXISTS (
     SELECT 1 FROM UNNEST(identifiedPTMStrings) t(ptm)
     WHERE lower(ptm.name) LIKE '%glcnac%' OR lower(ptm.name) LIKE '%hexnac%'
   )
ORDER BY accession;
