**Table 2.** Synthetic-sweep benchmark (8 configurations) 
underlying figure 1. Each row is a single synthetic study generated
by `make_synthetic_study(n_cells, n_genes, n_clusters = 4, density = 0.10)`
and benchmarked end-to-end via `bench_bundle()`: bundle + shard write,
`load_study()` cold-start, and 25 random `gene_values()` fetches.
Bundles stay near-constant in size while the shard tree grows linearly
with `n_cells × n_genes`.

| n_cells | n_genes | bundle MB | shard MB | prepare s | peak Mb | load s | fetch median ms | fetch max ms |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 500 | 2,000 | 2.58 | 1.1 | 32.8 | 256 | 0.117 | 15.2 | 17.5 |
| 500 | 5,000 | 6.21 | 2.6 | 90.1 | 289 | 0.141 | 14.6 | 21.6 |
| 1,000 | 2,000 | 2.62 | 2.1 | 35.3 | 281 | 0.116 | 16.1 | 20.6 |
| 1,000 | 5,000 | 6.27 | 5.2 | 93.5 | 282 | 0.133 | 16.4 | 20.3 |
| 2,000 | 2,000 | 2.71 | 4.1 | 38.6 | 300 | 0.120 | 19.3 | 52.9 |
| 2,000 | 5,000 | 6.36 | 10.3 | 104.0 | 310 | 0.145 | 19.1 | 25.6 |
| 4,000 | 2,000 | 2.91 | 8.2 | 44.2 | 296 | 0.119 | 23.0 | 26.2 |
| 4,000 | 5,000 | 6.55 | 20.5 | 121.4 | 317 | 0.137 | 22.7 | 38.5 |
