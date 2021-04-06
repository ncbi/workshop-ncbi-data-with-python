# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from collections import Counter

import ncbi.datasets as datasets
import pandas as pd

help(datasets)

api_instance = datasets.GenomeApi(datasets.ApiClient())
assembly_accessions = ["GCF_000001405.39"]
genome_summary = api_instance.assembly_descriptors_by_accessions(
    assembly_accessions, limit="all"
)
print(genome_summary)

print(f"Number of assemblies: {genome_summary.total_count}")

for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
    print(
        assembly.assembly_accession,
        assembly.assembly_level,
        len(assembly.chromosomes),
        assembly.submission_date,
        sep="\t",
    )

tax_name = "primates"
genome_summary = api_instance.assembly_descriptors_by_taxon(taxon=tax_name, limit="all")
assm_counter = Counter()
for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
    if assembly.assembly_accession[:3] == "GCA":
        assm_counter["GenBank"] += 1
    elif assembly.assembly_accession[:3] == "GCF":
        assm_counter["RefSeq"] += 1
df = pd.DataFrame.from_dict(assm_counter, orient="index", columns=["count"])
df.plot(kind="pie", y="count", figsize=(6, 6), title="Assemblies by type")
