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

# INTRO: TODO

# In this workshop you will learn how to:
# - Make computational "lab notebooks" with Jupyter
# - Manipulate and analyze data with Python programming
# -

# Motivation:
#
# 1. Do more stuff faster.
# 2. Have a reproducible record of what you did.
# 3. Have more confidence in your results.
# 4. Reuse, revise, extend your work.
#
# Discussion: problems with web-based work flow.

# ### What is a Jupyter notebook and what is it for?
# 1. A Jupyter notebook is a document that allows you to combine code, formatted text,
# and images
# 2. They are displayed and edited in a web browser
# 3. You can run the edit and run the code in place and display the output
# 4. They are useful for:
#     - Exploration: you can quickly test out ideas and see the results
#     - Documentation: a Jupyter notebook constitutes a record of precisely what you
# did. (Think of it as a "lab notebook" for your computational "experiments.")
#     - Communication: Jupyter notebooks make it easy to share what you did with
# colleagues (e.g. reports for your PI, interactive examples to accompany publications)

# ## Getting started with Jupyter

# A jupyter notebook is made up of blocks called *cells*.
# There are two types of cells: Markdown cells and code cells.

# ### Markdown cells
#
# This is a Markdown cell.
# Markdown cells contain formatted text and are useful for
# explaining yourself and organizing your document into sections.
#
# The text is formatted using the "Markdown" formatting language,
# which is an way to specify formatting like:
#
# - unordered lists
# - 1. nested lists
#   2. ordered lists
# - *italics*
# - **bold**
# - `code`
# - [hyperlinks](https://www.markdownguide.org/cheat-sheet/)
# - Latex equations: $E = mc^2$
#
# To see how Markdown works, double click this cell.
# You should see the unformatted text along with the Markdown formatting marks
# that tell Jupyter how to format the text.

# **Exercise**: Below is a Markdown cell for you to play with.
# Double click it to edit.
# When you've entered some text,
# press SHIFT-ENTER or the "Run" button at the top of the notebook
# to render the formatted text.
# Try out some different types of formatting, using the cell above or the
# [Markdown cheatsheet](https://www.markdownguide.org/cheat-sheet/)
# as a guide.

# EDIT ME.

# **Exercise**: Create a new Markdown cell below this one.
#
# - To create a cell, click the "+" button at the top of the notebook.
# - By default a new cell will be a code cell.
#   To make it Markdown, select the cell.
#   Look for the menu at the top that says "Code", click it and select
#   "Markdown" instead.
# - Test that it worked by adding some Markdown text to the cell and running it.
#
# (_Note: you will see that the cell-type menu contains two types of cells we won't
# talk about today, "Raw" and "Heading". "Raw" is rarely used and "Heading" is
# obsolete, so you don't need to worry about them._)

# ### Code cells
# Code cells contain Python code that you can run within the notebook.
# You can think of Python code as a list of instructions for the computer to follow.
# Possible instructions include: doing math, reading the contents of a file,
# searching a database, and displaying output on the screen.
# When you run a Code cell, a program called a "Python interpreter" reads the code
# in the cell, translates the code into instructions that the computer can understand,
# and then carries out the instructions.
# If there is any output to the code in a cell, it will be displayed in the notebook.

# Below are our first two code cells.
# Both cells contain a list of import statements
# (notice how each line contains the keyword `import`).
# Import statements are used to read ("import") Python code that other people
# have written into your notebook so that you can use it.
# This saves us from having to write everything from scratch.
# For example, in the second cell, we import code from the BioPython and
# NCBI Datasets packages.
# The former will let us manipulate sequence data, while the latter will let us
# search for and download data from NCBI.

# We will not be talking about how importing works in detail today.
# You *do not need to understand* how the code in these cells work.
# In case you're curious, we've added comments describing what each imported
# module does.
# (Comments are lines of code that begin with "#".
# Python ignores these lines completely.)

# Select and run both of the cells below.
# (Remember: it's SHIFT-ENTER or the "Run" button to run a cell.)
# Note the text to the left of each cell that says `In [ ]:`.
# Watch what happens when you run the cells.

# Modules from the Python Standard Library:
# "Regular Expressions" for searching and replacing text
import re
# running programs outside of Python, in our case BLAST
import subprocess
# creating temporary files to hold intermediate results that we don't care about saving
import tempfile
# stands for Input/Output
from io import BytesIO, StringIO
# specifying the locations of files on our computer
from pathlib import Path
# specifying the types of function arguments and return values
from typing import Dict, List, Tuple
# reading and writing zip files
from zipfile import ZipFile

# Modules from other packages:
# Matplotlib: plotting figures
import matplotlib.pyplot as plt
# BioPython: reading and writing sequence files,
# manipulating and aligning sequences, making phylogenetic trees
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import (NNITreeSearcher, ParsimonyScorer,
                                        ParsimonyTreeConstructor)
from Bio.SeqRecord import SeqRecord
# NCBI Datasets: searching and downloading NCBI data
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets.openapi import GeneApi as DatasetsGeneApi
from ncbi.datasets.openapi import GenomeApi as DatasetsGenomeApi
from ncbi.datasets.openapi.models import AssemblyDatasetRequestResolution
from ncbi.datasets.package import dataset

# Importing is important, but not very exciting.
# After all, there was no output when we ran those cells above.
# The true value of a Jupyter notebook is being able to perform computations
# and then see the result right away in the same document.
#
# Let's do a little math to illustrate how it works.
# The cell below performs some simple arithmetic.
#
# Run it and see what happens.

(5 + 3 ** 2) * 10

# **Question**: What do you think the `**` does?
#
# **Exercise**: Make a new code cell.
#
# - Create a new cell with the "+" button at the top of the notebook.
# - Double-click it to edit it.
# - Try out your own mathematical expression. Predict what the answer will be.
# - Run the cell.
# - Was the output what you expected?
# - Repeat with a different expression.

# Numerical output is fine, but as scientists, we often want to examine or share
# our data as figures.
# Jupyter lets us make figures that are directly embedded in the notebook next to
# the code that generated them.
# You can think of the code as an automatic "caption" for the figure,
# showing exactly how it was made.
# You can also add markdown cells interpreting the results.
#
# Here is a simple plot. Don't worry about the details yet, just run the cell:

# All the integers from -10 to 10
xs = range(-10, 11)
# y = x ** 2 for each value of x
ys = [x ** 2 for x in xs]
# Plot y vs. x, using blue circles
plt.plot(xs, ys, color="blue", marker="o")
# Set the axis labels and figure title.
plt.xlabel(r"The x-axis: $x$")
plt.ylabel(r"The y-axis: $y = x ^ 2$")
plt.title("A simple plot")

# **Exercise**: Make your own plot:
#
# - Copy the previous cell. Select it and use the "copy" button at the
# top of the notebook.
# - Paste the cell below this one. Select this cell and then click the
# "paste" buton at the top of the notebook.
# - Edit your new cell and run it.
# - Add a markdown cell below the plot explaining what you did.
#
# Some things to try:
# - Change the color to "red".
# - Change the marker to "x".
# - Change the range of x values.
# - Change the y values by changing `x ** 2` to some other mathematical expression.
#
# If you get an error or would like to undo a change, select the cell and
# use `CTRL-Z` (`CMD-Z` on a mac) to undo the most recent change.

# ### Under the hood
# TODO

# The code is run by a program called the kernel
# - Restarting the kernel. Sometimes, you will want the program to "forget" the results
# of the cells you've run and start fresh. EXAMPLE of restarting the kernel
# - Interrupting the kernel. Sometimes a cell is taking too long to run, either because
# you made a mistake or because the task is bigger than you expected.
# EXAMPLE of interrupting an infinite loop

# To run a Jupyter notebook, you need a Jupyter server, either:
# 1. Remotely: i.e., on the internet; what we're doing today LINKS TO ONLINE SERVICES
# 2. Locally: i.e., on your own computer LINK TO INSTALLATION INSTRUCTIONS
# (today we won't cover installation because it may vary based on your system, but it's
# not difficult.)

# ## Getting started with Python

# Variables: giving names to values

# Types: how data is represented
# - Int
# - Float
# - String
# - Boolean

# Functions: transforming input to output

# Data structures: storing lots of values
# - Lists
# - Dicts


# ## Hello, world!
# It is traditional when learning a new programming language to start with a simple program called "Hello, world!" that prints a greeting. This simple program illustrates many of the important features of Python that we're going to learn about today.
#
# Here is "Hello, world!" in Python. Run the program by selecting the cell and hitting SHIFT+ENTER.

# This line is a comment.
#  (1)  (2)      (3)
greeting = "Hello, world!"
# (4)    (5)
print(greeting)

# - Warning: you can run the cells in any order, but you should try not to. When you're done working, you should restart the notebook and run the cells in order to make sure that it works the way you expect. (this may be out of place)

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

# ## Conclusions and follow-up resources
# TODO
# - Git repo + binder link
# - Export pdf?
