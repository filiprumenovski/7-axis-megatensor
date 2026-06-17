-- Canon reference tensor (identity backbone demo)
SELECT c.*, i.protein_acc, i.residue_pos, i.residue_aa, m.metric_name, m.metric_value
FROM read_parquet('/run/media/filip/Data/7-axis-ptm-tensor/megatensor/canon/sets/set_coordinates/dataset_id=*/*.parquet', hive_partitioning := true) c
JOIN read_parquet('/run/media/filip/Data/7-axis-ptm-tensor/megatensor/canon/registry/identity_dim.parquet') i USING (identity_id)
JOIN read_parquet('/run/media/filip/Data/7-axis-ptm-tensor/megatensor/canon/metrics/set_metrics/dataset_id=*/*.parquet', hive_partitioning := true) m USING (set_uid, dataset_id);
