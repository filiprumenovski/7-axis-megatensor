-- PRIDE experimental tensor (rich context axes)
SELECT c.*, i.protein_acc, i.residue_pos, i.residue_aa,
       id.inst_model_token, pd.provenance_token, m.metric_name, m.metric_value
FROM read_parquet('/run/media/filip/Data/7-axis-ptm-tensor/megatensor/pride/sets/set_coordinates/dataset_id=*/*.parquet', hive_partitioning := true) c
JOIN read_parquet('/run/media/filip/Data/7-axis-ptm-tensor/megatensor/pride/registry/identity_dim.parquet') i USING (identity_id)
JOIN read_parquet('/run/media/filip/Data/7-axis-ptm-tensor/megatensor/pride/registry/instrument_dim.parquet') id ON c.instrument_id = id.instrument_id
JOIN read_parquet('/run/media/filip/Data/7-axis-ptm-tensor/megatensor/pride/registry/provenance_dim.parquet') pd ON c.provenance_id = pd.provenance_id
JOIN read_parquet('/run/media/filip/Data/7-axis-ptm-tensor/megatensor/pride/metrics/set_metrics/dataset_id=*/*.parquet', hive_partitioning := true) m USING (set_uid, dataset_id);
