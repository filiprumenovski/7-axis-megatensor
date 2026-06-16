-- Megatensor union view (DuckDB). Point paths at megatensor/ under repo root.
CREATE OR REPLACE VIEW megatensor AS
SELECT
  c.set_uid,
  c.dataset_id,
  i.protein_acc,
  i.isoform,
  i.residue_pos,
  i.residue_aa,
  i.protein_level_only,
  p.ptm_label,
  p.ptm_unimod,
  c.loc_score,
  c.loc_method,
  c.loc_is_ambiguous,
  cd.condition_token,
  ad.acquisition_token,
  id.instrument_token,
  pd.provenance_token,
  m.metric_name,
  m.metric_value,
  m.metric_norm_state,
  m.metric_unit,
  m.qc_flags
FROM read_parquet('megatensor/sets/set_coordinates/dataset_id=*/*.parquet', hive_partitioning := true) c
JOIN read_parquet('megatensor/registry/identity_dim.parquet') i ON c.identity_id = i.identity_id
JOIN read_parquet('megatensor/registry/ptm_dim.parquet') p ON c.ptm_id = p.ptm_id
JOIN read_parquet('megatensor/registry/condition_dim.parquet') cd ON c.condition_id = cd.condition_id
JOIN read_parquet('megatensor/registry/acquisition_dim.parquet') ad ON c.acquisition_id = ad.acquisition_id
JOIN read_parquet('megatensor/registry/instrument_dim.parquet') id ON c.instrument_id = id.instrument_id
JOIN read_parquet('megatensor/registry/provenance_dim.parquet') pd ON c.provenance_id = pd.provenance_id
JOIN read_parquet('megatensor/metrics/set_metrics/dataset_id=*/*.parquet', hive_partitioning := true) m
  ON c.set_uid = m.set_uid AND c.dataset_id = m.dataset_id;
