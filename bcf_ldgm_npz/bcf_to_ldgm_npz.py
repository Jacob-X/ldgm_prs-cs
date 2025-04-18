import pysam
import numpy as np
import os

def parse_bcf_to_ldgm_graph(bcf_path):
    vcf = pysam.VariantFile(bcf_path)
    positions = []
    diagonals = []
    neighbors = []
    weights = []

    for rec in vcf.fetch():
        info = rec.info
        positions.append(rec.pos)
        diagonals.append(float(info.get("LD_diagonal", 1.0)))

        neighbor_str = info.get("LD_neighbors")
        weight_str = info.get("LD_weights")

        if neighbor_str is not None:
            ns = list(map(int, neighbor_str.split(","))) if isinstance(neighbor_str, str) else list(neighbor_str)
        else:
            ns = []

        if weight_str is not None:
            ws = list(map(float, weight_str.split(","))) if isinstance(weight_str, str) else list(weight_str)
        else:
            ws = []

        neighbors.append(np.array(ns, dtype=np.int32))
        weights.append(np.array(ws, dtype=np.float32))

    return {
        "positions": np.array(positions, dtype=np.int32),
        "diagonals": np.array(diagonals, dtype=np.float32),
        "neighbors": neighbors,
        "weights": weights
    }

def save_ldgm_npz(graph, output_path):
    np.savez_compressed(
        output_path,
        positions=graph["positions"],
        diagonals=graph["diagonals"],
        neighbors=np.array(graph["neighbors"], dtype=object),
        weights=np.array(graph["weights"], dtype=object)
    )

def convert_all_bcf_to_npz(input_dir, output_dir, populations):
    os.makedirs(output_dir, exist_ok=True)
    for pop in populations:
        bcf_path = os.path.join(input_dir, f"1kg_ldgm.{pop}.bcf")
        output_path = os.path.join(output_dir, f"1kg_ldgm.{pop}.npz")
        print(f"Processing {pop}...")

        graph = parse_bcf_to_ldgm_graph(bcf_path)
        save_ldgm_npz(graph, output_path)

        print(f"Saved to {output_path}")

if __name__ == '__main__':
    input_dir = "/Volumes/data_files/LDGM/example/height/ldgms.GRCh38"
    output_dir = "/Volumes/data_files/LDGM/example/height/ldgm_npz"
    populations = ["AFR", "AMR", "EAS", "EUR", "SAS"]

    convert_all_bcf_to_npz(input_dir, output_dir, populations)