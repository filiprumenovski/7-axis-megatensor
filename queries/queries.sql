-- §11.4 Live query demo (DuckDB). Run from repo root:
--   duckdb < queries/queries.sql
-- Or: duckdb -c ".read queries/queries.sql"

-- Headline: one O-GlcNAc site observed across canon + PRIDE sources
-- (Q9UHR5:233:T — Figure C candidate, Raptor S233)
WITH site AS (
  SELECT * FROM read_parquet('megatensor/union/staging/site_index.parquet')
  WHERE site_key = 'Q9UHR5:233:T'
)
SELECT * FROM site;

-- Same site: all PRIDE observations with provenance axes
SELECT
  o.dataset_id,
  o.protein_id_raw,
  o.residue_pos_raw,
  o.residue_aa,
  o.cond_treatment,
  o.metric_name,
  o.metric_value,
  o.loc_score,
  o.loc_method,
  o.inst_model,
  o.prov_country,
  o.source_engine
FROM read_parquet('megatensor/pride/staging/observations.parquet') o
WHERE o.protein_id_raw = 'Q9UHR5'
  AND o.residue_pos_raw = 233
  AND o.residue_aa = 'T'
ORDER BY o.dataset_id, o.cond_treatment, o.metric_name;

-- Five labs / five instruments / one site (sites in both layers)
SELECT
  site_key,
  layers,
  datasets,
  set_hits
FROM read_parquet('megatensor/union/staging/site_index.parquet')
WHERE list_contains(layers, 'canon') AND list_contains(layers, 'pride')
ORDER BY set_hits DESC
LIMIT 10;

-- Per-instrument slice (PRIDE)
SELECT
  inst_model,
  count(DISTINCT concat(protein_id_raw, ':', residue_pos_raw::VARCHAR, ':', residue_aa)) AS sites,
  count(*) AS obs_rows
FROM read_parquet('megatensor/pride/staging/observations.parquet')
WHERE inst_model IS NOT NULL
GROUP BY 1
ORDER BY obs_rows DESC;

-- Enriched site features for shared sites
SELECT e.*
FROM read_parquet('megatensor/union/enrichment/site_features.parquet') e
WHERE concat(e.protein_acc, ':', e.residue_pos::VARCHAR, ':', e.residue_aa) IN (
  SELECT site_key FROM read_parquet('megatensor/union/staging/site_index.parquet')
  WHERE n_layers >= 2
)
LIMIT 20;
