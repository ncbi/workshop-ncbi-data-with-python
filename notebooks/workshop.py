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

# # Getting started with NCBI data in Python
# In this workshop you will learn how to:
# - Make computational "lab notebooks" with Jupyter
# - Manipulate and analyze data with Python programming
# -

# Motivation:
# 1. Do more stuff faster.
# 2. Have a reproducible record of what you did.
# 3. Reuse, revise, extend your work.
#
# Discussion: problems with Excel.

# ## This is a Jupyter notebook

# ### What is a Jupyter notebook and what is it for?
# 1. A Jupyter notebook is a document that allows you to combine code, formatted text, and images
# 2. They are displayed and edited in a web browser
# 3. You can run the edit and run the code in place and display the output
# 4. They are useful for:
#     - Exploration: you can quickly test out ideas and see the results
#     - Documentation: a Jupyter notebook constitutes a record of precisely what you did. (Think of it as a "lab notebook" for your computational "experiments.")
#     - Communication: Jupyter notebooks make it easy to share what you did with colleagues (e.g. reports for your PI, interactive examples to accompany publications)
#
#
# ### How to use a Jupyter notebook
# To run a Jupyter notebook, you need a Jupyter server, either:
# 1. Remotely: i.e., on the internet; what we're doing today LINKS TO ONLINE SERVICES
# 2. Locally: i.e., on your own computer LINK TO INSTALLATION INSTRUCTIONS (today we won't cover installation because it may vary based on your system, but it's not difficult.)
# A Jupyter notebook is made up of discrete blocks called "cells"
# - Basic cell operation: edit, move, make new ones, run, changing cell type
# - Markdown cells are for text describing what you're doing
# - They're useful for explanation, notes-to-self and to dividing your document into sections
# - You can use Markdown syntax to add formatting to your text LINK TO MARKDOWN CHEATSHEET
# EXAMPLES: heading, numbered list
#
# Code cells are for code that you'd like to run
# - Running an existing block
# - Modifying code and re-running EXAMPLE: e.g., "Hello world"
# - Warning: you can run the cells in any order, but you should try not to. When you're done working, you should restart the notebook and run the cells in order to make sure that it works the way you expect. (this may be out of place)
#
# Raw cells are used when converting the notebook to other file types with the `nbconvert` utility. You don't need to worry about these today (or likely ever).
#
# The code is run by a program called the kernel
# - Restarting the kernel. Sometimes, you will want the program to "forget" the results of the cells you've run and start fresh. EXAMPLE of restarting the kernel
# - Interrupting the kernel. Sometimes a cell is taking too long to run, either because you made a mistake or because the task is bigger than you expected. EXAMPLE of interrupting an infinite loop

# ## Getting started with Python

# Variables: giving names to values

# Functions: transforming input to output

# Types: how data is represented
# - Int
# - Float
# - String
# - Boolean

# Data structures: storing lots of values
# - Lists
# - Dicts
# - pandas dataframes

# Modules: using other peoples' code

# ## Codon usage

import re
import subprocess
import tempfile
from io import BytesIO, StringIO
from pathlib import Path
from typing import Dict, List, Tuple
from zipfile import ZipFile

from Bio import AlignIO, Phylo, SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import (NNITreeSearcher, ParsimonyScorer,
                                        ParsimonyTreeConstructor)
from Bio.SeqRecord import SeqRecord
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets.openapi import GeneApi as DatasetsGeneApi
from ncbi.datasets.openapi import GenomeApi as DatasetsGenomeApi
from ncbi.datasets.openapi.models import AssemblyDatasetRequestResolution
from ncbi.datasets.package import dataset

# +
# ### Manipulating sequences
# - Make a string that could be a DNA sequence
# - Make Seq object out of string
# - Various manipulation exercises
#     - Complement/Reverse complement: Predict by hand, apply function, test equality
#     - Transcribe: predict by hand, apply function, test equality

# ### Project
# Motivating question: how do codon usage patterns vary between species?
#
# #### Counting codon usage
# (First with toy dataset)
# IO -> List[Seq]
#
# List[Seq] -> Dict[Codon, Int]
#
# #### Plot usage
#
# Dict[Codon, Int] -> Bar chart
#
# #### Download and analyze metadata
#
# #### Download data
#
# [Species] -> Data
#
# #### Parse data
#
# Data -> Dict[Taxon, List[Seq]]
#
# #### Apply functions to data
# Plots comparing multiple species' usage.

# ## Less-interactive bigger example

# ## Finding NCBI data

# 1. Get summary for all whale genome assemblies
#     - How many assemblies are there?
#     - How many of the assemblies have annotation?
#     - How many assemblies are RefSeq assemblies?
#     - What is the earliest submission date? 2012/01/02
#     - Make a table of Common Name | Sci Name | GenBank Assem. Accn. | RefSeq Assem. Accn.
# 2. Get the summary the genes for the sperm whale
#     - How many genes are there?
#     - Get the gene summary for sperm whale MB
#     - NM_001290722.1
#     - Get transcripts for all whale MB genes (use orthologs in datasets)
#     - Align them using MUSCLE?

# ### API Calls
#
# Getting the ortholog summary for MYH11
# ```bash
# datasets summary ortholog symbol MYH11 > MYH11_orthologs_summary.json
# ```
#
# Orthologs for Human, Elephant shark, Zebrafish, Zenopus
# ```bash
# curl 'https://api.ncbi.nlm.nih.gov/datasets/v1alpha/gene/id/4629/orthologs?taxon_filter=human&taxon_filter=elephant%20shark&taxon_filter=zebrafish&taxon_filter=Xenopus%20tropicalis' > MYH1_orthologs_human_shark_fish_frog.json
# ```
#
# Downloading sequences:
# ```bash
# curl -X POST "https://api.ncbi.nlm.nih.gov/datasets/v1alpha/gene/download" \
#  -H "Accept: application/zip, application/json" \
#  -H "Content-Type: application/json" \
#  -d '{"gene_ids":{"gene_ids":[59067,50615]}}'
#  ```


# # Provide your own gene ids as a list of integers
# input_gene_ids: List[int] = [1, 2, 3, 9, 10, 11, 12, 13, 14, 15, 16, 17]


# def gene_infor_for_geneids(gene_ids: List[int]):
#     if len(gene_ids) == 0:
#         print("Please provide at least one gene-id")
#         return

#     with DatasetsApiClient() as api_client:
#         gene_api = DatasetsGeneApi(api_client)
#         try:
#             gene_reply = gene_api.gene_metadata_by_id(gene_ids)
#             for gene in gene_reply.genes:
#                 print_gene_metadata_by_fields(gene, fields=["gene_id", "symbol"])
#         except DatasetsApiException as e:
#             print(f"Exception when calling GeneApi: {e}\n")


# def example_usage_of_api(gene_ids: List[int]):
#     if len(gene_ids) == 0:
#         print("Please provide at least one gene-id")
#         return

#     with DatasetsApiClient() as api_client:
#         gene_api = DatasetsGeneApi(api_client)

#         # Get just metadata
#         try:
#             gene_reply = gene_api.gene_metadata_by_id(gene_ids)
#             for gene in gene_reply.genes:
#                 print(gene.gene.gene_id)
#         except DatasetsApiException as e:
#             print(f"Exception when calling GeneApi: {e}\n")

#         # Or, download a data package with FASTA files
#         try:
#             print("Begin download of data package ...")
#             gene_ds_download = gene_api.download_gene_package(
#                 gene_ids, include_annotation_type=["FASTA_GENE"], _preload_content=False
#             )
#             gene_reply = gene_api.gene_metadata_by_id(gene_ids)
#             zipfile_name = "gene_ds.zip"

#             with open(zipfile_name, "wb") as f:
#                 f.write(gene_ds_download.data)
#             print(f"Download completed -- see {zipfile_name}")

#         except DatasetsApiException as e:
#             print(f"Exception when calling GeneApi: {e}\n")


# -

# # Begin here

# ## Hello, world!
# It is traditional when learning a new programming language to start with a simple program called "Hello, world!" that prints a greeting. This simple program illustrates many of the important features of Python that we're going to learn about today.
#
# Here is "Hello, world!" in Python. Run the program by selecting the cell and hitting SHIFT+ENTER.

# This line is a comment.
#  (1)  (2)      (3)
greeting = "Hello, world!"
# (4)    (5)
print(greeting)

# - The lines that begin with "#" are "comments". They're a way to write notes about your code that you want Python to ignore. **Exercise**: Try editing the cell to add another line of comments somewhere. Run the cell again. Did anything change about the output?

name_to_greet = "World"
#                    (1)     (2)
greeting2 = "Hello, " + name_to_greet + "!"
print(greeting2)


# +
# (1) (2)   (3)  (4)     (5)
def greet(name: str) -> str:
    # (6)       (7)     (8) (9)
    greeting = "Hello, " + name + "!"
    # (10)    (11)
    return greeting


my_name = "Your name here"
#            (12)  (13)
greeting3 = greet(my_name)
print(greeting3)


# -

# Here's what's going on in that cell:
# - The `def` keyword (1) begins the definition of a *function*. A function takes one or more *arguments* as input, does some computation, and *returns* some output.
# - The name of the function is `greet` (2).
# - `greet` takes one argument called `name` (3). This is the name of the person we would like to greet.
# - `name` should have the "type" `str` (4).

# +
def make_scientific_name(genus: str, species: str) -> str:
    # Fix me (1):
    sci_name = ""
    # Fix me (2):
    return ""


my_genus = "Homo"
my_species = "sapiens"
# Fix me (3):
my_scientific_name = ""
print("The scientific name for humans is: " + my_scientific_name)
# -

list_of_names = ["Charles Darwin", "Alfred Russel Wallace", "Rosalind Franklin"]
for person in list_of_names:
    print(greet(person))

# **Exercise**: add another person or two to the list of names to greet.

list_of_species = [("Bos", "taurus"), ("Homo", "sapiens"), ("Orcinus", "orca")]
for genus, species in list_of_species:
    # Fix me
    print("")

list_of_greetings = [greet(p) for p in list_of_names]

list_of_scientific_names = ["" for genus, species in list_of_species]
print(list_of_scientific_names)


# ## Get GeneID


def get_gene_ids(
    gene_api: DatasetsGeneApi, gene_symbols: List[str], taxon: str
) -> List[int]:
    try:
        md = gene_api.gene_metadata_by_tax_and_symbol(gene_symbols, taxon)
        return [int(g.gene.gene_id) for g in md.genes]
    except DatasetsApiException as e:
        print(f"Exception when calling GeneApi: {e}\n")
    return []


gene_api = DatasetsGeneApi()
gene_ids = get_gene_ids(gene_api, ["MB"], "human")


# ## Ortholog Query


def query_orthologs(gene_api: DatasetsGeneApi, gene_id: int, taxa: List[str]) -> List:
    try:
        response = gene_api.gene_orthologs_by_id(gene_id, taxon_filter=taxa)
        return [item.gene for item in response.genes.genes]
    except DatasetsApiException as e:
        print(f"Exception when calling GeneApi: {e}\n")
    return []


gene_id = gene_ids[0]
taxa = ["human", "whales", "Bos taurus"]
mb_orthologs = query_orthologs(gene_api, gene_id, taxa)
print(mb_orthologs)


for g in mb_orthologs:
    print(f"{g.common_name}:\t{g.taxname}")

orth = mb_orthologs[0]
print(orth.common_name)


def get_longest_transcripts(genes: List) -> Dict:
    ret = {}
    for g in genes:
        longest = max(g.transcripts, key=lambda x: x.length)
        ret[longest.accession_version] = longest
    return ret


ortholog_gene_ids = [int(g.gene_id) for g in mb_orthologs]
print(ortholog_gene_ids)

# ## Download Gene Data


def download_transcripts(
    gene_api: DatasetsGeneApi, gene_ids: List[int], data_dir: Path
) -> str:
    if len(gene_ids) == 0:
        print("Please provide at least one gene-id")
        return ""
    try:
        print("Begin download of data package ...")
        gene_ds_download = gene_api.download_gene_package(
            gene_ids, include_annotation_type=["FASTA_RNA"], _preload_content=False
        )
        with ZipFile(BytesIO(gene_ds_download.data)) as zipfile:
            data_file_name = zipfile.extract("ncbi_dataset/data/rna.fna", path=data_dir)
        print(f"Download completed -- see {data_file_name}")
        return data_file_name
    except DatasetsApiException as e:
        print(f"Exception when calling GeneApi: {e}\n")
    return ""


data_dir = Path("../data")
# !tree {data_dir}

ortholog_fasta = download_transcripts(
    gene_api, ortholog_gene_ids, data_dir / "ortholog_dataset/"
)
# !tree {data_dir}

# +
r = re.compile(r"\[organism=(?P<name>[\w ]+)\]")


def get_organism_name(record: SeqRecord) -> str:
    m = re.search(r, record.description)
    if m:
        return m.group("name")
    else:
        return "NameNotFound"


def get_transcript_dict(rna_fasta: str) -> Dict[str, SeqRecord]:
    return {record.id: record for record in SeqIO.parse(rna_fasta, "fasta")}


# -

transcripts = get_transcript_dict(ortholog_fasta)
for rec_id, rec in transcripts.items():
    print(rec_id, get_organism_name(rec), len(rec.seq))


# ## Get CDS

# - Poly-a tail in transcript
# - CDS starts with Start Codon


def get_cds_region(transcript):
    transcript_range = transcript.cds.range[0]
    return (int(transcript_range.begin), int(transcript_range.end))


longest_transcripts = get_longest_transcripts(mb_orthologs)
cds_regions = {
    accession: get_cds_region(transcript)
    for accession, transcript in longest_transcripts.items()
}
print(cds_regions)

cds_records = {}
for accession, (start, end) in cds_regions.items():
    transcript_record = transcripts[accession]
    cds_record = SeqRecord(
        id=transcript_record.id,
        name=transcript_record.name,
        description=transcript_record.description,
        seq=transcript_record.seq[start - 1 : end],
    )
    cds_records[accession] = cds_record

for rec in cds_records.values():
    print(rec.seq.translate())

# ## Alignment


cds_fasta = data_dir / "cds.fasta"
SeqIO.write(cds_records.values(), cds_fasta, "fasta")

for x in SeqIO.parse(cds_fasta, "fasta"):
    print(x)


muscle_exe = Path("../bin/muscle3.8.31_i86linux64")
muscle_cline = MuscleCommandline(muscle_exe, input=cds_fasta)
print(muscle_cline)
stdout, stderr = muscle_cline()


align = AlignIO.read(StringIO(stdout), "fasta")
print(align)

print(align.substitutions)


scorer = ParsimonyScorer()
searcher = NNITreeSearcher(scorer)
constructor = ParsimonyTreeConstructor(searcher)
pars_tree = constructor.build_tree(align)

print(pars_tree)


Phylo.draw(pars_tree)

Phylo.draw_ascii(pars_tree)


# ## Finding MB in an unannotated assembly

# +
def get_assembly_accessions(
    genome_api: DatasetsGenomeApi,
    taxon: str,
) -> List[str]:
    try:
        assembly_reply = genome_api.assembly_descriptors_by_taxon(taxon)
        return [
            assembly.assembly.assembly_accession
            for assembly in assembly_reply.assemblies
        ]
    except DatasetsApiException as e:
        print(f"Exception when calling GenomeApi: {e}\n")
    return []


def get_wgs_project_accessions(
    genome_api: DatasetsGenomeApi,
    genome_assembly_accessions: List[str],
) -> List[str]:
    try:
        genome_ds_download = genome_api.download_assembly_package(
            genome_assembly_accessions,
            exclude_sequence=True,
            hydrated=AssemblyDatasetRequestResolution.DATA_REPORT_ONLY,
            _preload_content=False,
        )
        with tempfile.NamedTemporaryFile(mode="w+b") as f:
            f.write(genome_ds_download.data)
            f.flush()
            package = dataset.AssemblyDataset(f.name)
            return [
                report.wgs_info.wgs_project_accession
                for report in package.get_data_reports()
            ]
    except DatasetsApiException as e:
        print(f"Exception when calling GenomeApi: {e}\n")
    return []


# -

taxon = "Kogia breviceps"
genome_api = DatasetsGenomeApi()
ga_acc = get_assembly_accessions(genome_api, taxon)
print(ga_acc)
wgs_acc = get_wgs_project_accessions(genome_api, ga_acc)
print(wgs_acc)
k_breviceps_accession = wgs_acc[0]

blast_binary = Path("../bin/blastn_vdb")


# +
sperm_whale_fasta = data_dir / "sperm_whale_mb.fasta"
sperm_whale_name = "Physeter catodon"

for rec_id in longest_transcripts:
    rec = transcripts[rec_id]
    if get_organism_name(rec) == sperm_whale_name:
        SeqIO.write(rec, sperm_whale_fasta, "fasta")
        break
# -

# blastn_vdb binary to blast sperm whale myoglobin against unannotated kagia genome.


def run_blastn_vdb(
    blast_binary: Path, database_accession: str, query: Path, evalue: float = 0.01
) -> Tuple[str, str]:
    command = [
        str(blast_binary),
        "-task",
        "blastn",
        "-db",
        str(database_accession),
        "-query",
        str(query),
        "-evalue",
        str(evalue),
    ]
    print(f"Executing command:\n{' '.join(command)}")
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return (result.stdout.decode("utf-8"), result.stderr.decode("utf-8"))


blast_result, blast_error = run_blastn_vdb(
    blast_binary, k_breviceps_accession, sperm_whale_fasta
)
print(blast_result)

print(blast_error)

# ## TODO:
# - Get whale summary: eg how many genomes does NCBI have? What species? Common names? Submitted?
# - Gotcha: there are older assemblies in it. There are 27 representative genomes for whales. Some don't have annotation.
# - [x] Make an alignment (w/ BioPython)
# - [x] Make a tree from the alignment
# - ~Codon usage (lumped vs. individual species)~
# - [x] Whale myoglobin (with tree). Need to get all the refseq. Dataset is too big.
# - [x] Search with gene symbol ("MB"), rather than gene id. Get Gene ID.

# First thing to do: set the stage by looking up: how many whale genomes, what are the species?
# RefSeq assemblies can get sequences to work with. Whales are diving mammals, myglobin for holding their breath in muscles. Get sequences. Make tree of whale genomes. Include human cow as outgroups.(Cow in same orders). Tree doesn't quite map with the phylogeny. Have unannotated Pygmy sperm whales genome. Do a blast search.
