import pandas as pd
import point_grid as pg


internal_count = 3
edge_count = 6


df = pd.read_csv("data/6-gon-mesh-with-internal-nodes.csv")
df = df[df.target_edge_length == 0.4]
df = df[df.internal_node_count != 0]
dataset = df.drop(["target_edge_length", "internal_node_count"], axis=1)

# Csv-file to write to.
new_csv_path = "data/8-gon-patch-test.csv"
pc = pg.generate_patch_collection(dataset, edge_count=edge_count)
pg.write_patch_collection_to_csv(pc, new_csv_path, edge_count)
