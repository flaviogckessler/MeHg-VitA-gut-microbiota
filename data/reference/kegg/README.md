# KEGG REST reference results

[`kegg_api_find_ko_mercuric_orthologs.tsv`](kegg_api_find_ko_mercuric_orthologs.tsv)
is the unchanged tabular result of the manual KEGG REST `find/ko/mercuric`
query used to identify the orthologs considered for Figure S18. The supplied
file was renamed from `KEGG_rest.kegg.jp_find_ko_mercuric.txt` and moved out of
the figure directory so that small external references are stored under
`data/reference/`.

The query result lists seven KOs. Five occur in the batch 01 PICRUSt2 KO
abundance matrix, and the four used in Figure S18 (`K00520`, `K08363`, `K08364`,
and `K08365`) are also present in the batch 02 matrix. The official generation
script validates the query schema, both canonical KO matrices, and this explicit
four-KO cross-batch subset before plotting.
