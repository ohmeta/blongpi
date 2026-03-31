#!/usr/bin/env python3
import networkx as nx
import pandas as pd
import os

# blongpi - Step 5: Operon Integrity Analysis
# Principle: Check if candidate genes are physically adjacent in the pangenome graph.

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRAPH_PATH = os.path.join(BASE_DIR, "results/pangenome/final_graph.gml")
OUT_PATH = os.path.join(BASE_DIR, "results/operon_integrity.csv")

# Targeted operons to check
# Format: { 'OperonName': ['Gene1', 'Gene2', 'Gene3'] }
TARGET_OPERONS = {
    'MtrAB-LpqB': ['mtrA', 'mtrB', 'lpqB'],
    'FructosePTS': ['fruG', 'fruF', 'fruK', 'fruE'],
    'TadPili': ['tadA', 'tadB', 'tadC', 'flp']
}

def analyze_operons():
    if not os.path.exists(GRAPH_PATH):
        print(f"Graph not found at {GRAPH_PATH}. Run Step 1 first.")
        return

    print(">>> Loading Pangenome Graph...")
    G = nx.read_gml(GRAPH_PATH)
    
    results = []
    
    for op_name, genes in TARGET_OPERONS.items():
        print(f"Checking integrity of {op_name}...")
        # Find node IDs for these genes
        node_map = {}
        for node, data in G.nodes(data=True):
            gene_name = data.get('label', '')
            if any(g in gene_name for g in genes):
                # Map standard name to Panaroo node ID
                matched_gene = [g for g in genes if g in gene_name][0]
                node_map[matched_gene] = node
        
        # Check adjacency
        found_genes = list(node_map.keys())
        is_intact = len(found_genes) == len(genes)
        
        # Check if they are connected in the graph
        connections = 0
        if is_intact:
            for i in range(len(genes)-1):
                if G.has_edge(node_map[genes[i]], node_map[genes[i+1]]):
                    connections += 1
        
        integrity_score = connections / (len(genes)-1) if len(genes) > 1 else 1.0
        results.append({
            'Operon': op_name,
            'Genes_Found': ",".join(found_genes),
            'Integrity_Score': integrity_score
        })

    df = pd.DataFrame(results)
    df.to_csv(OUT_PATH, index=False)
    print(f">>> Operon analysis saved to {OUT_PATH}")

if __name__ == "__main__":
    analyze_operons()
