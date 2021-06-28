# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Getting started with NCBI data in Python

# Dan Rice and Peter Cooper

# In this workshop you will learn how to:
# - Use Python programming to download, analyze, and visualize data.
# - Use Jupyter to create data analysis 'lab notebooks' that make it easy to reproduce
# and share what you did.
# - Find data that is relevant to your project using the new NCBI Datasets resource.
# - Download NCBI sequence data and manipulate it with the BioPython package.

# Why learn Python/Jupyter as a biologist:
#
# 1. Do more work faster through automation.
# 2. Have a reproducible record of what you did.
# 3. Have more confidence in your results.
# 4. Reuse, revise, extend your work.
# 5. Share your work with colleagues.

# ## Getting started with Jupyter

# ### What is a Jupyter notebook and what is it for?
# 1. A Jupyter notebook is a document that allows you to combine code, formatted text,
# and images.
# 2. Notebooks are displayed and edited in a web browser.
# 3. You can run the edit and run the code in place and display the output.
# 4. They are useful for:
#     - Exploration: you can quickly test out ideas and see the results
#     - Documentation: a Jupyter notebook constitutes a record of precisely what you
# did. (Think of it as a "lab notebook" for your computational "experiments.")
#     - Communication: Jupyter notebooks make it easy to share what you did with
# colleagues (e.g. reports for your PI, interactive examples to accompany publications)
#
# A jupyter notebook is made up of blocks called *cells*.
# There are two types of cells: **Markdown cells** and **code cells**.

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

# Below is our first code cell.
# It contains a list of import statements
# (notice how each line contains the keyword `import`).
# Import statements are used to read ("import") Python code that other people
# have written into your notebook so that you can use it.
# This saves us from having to write everything from scratch.
# For example, at the bottom of of the cell, we import code from the BioPython and
# NCBI Datasets packages.
# The former will let us manipulate sequence data, while the latter will let us
# search for and download data from NCBI.

# We will not be talking about how importing works in detail today.
# You *do not need to understand* how the code in this cell works.
# In case you're curious, we've added comments describing what each imported
# module does.
# (Comments are lines of code that begin with "#".
# Python ignores these lines completely.)

# **Select and run the cell below.**
# (Remember: it's SHIFT-ENTER or the "Run" button to run a cell.)
# Note the text to the left of the cell that says `In [ ]:`.
# Watch what happens when you run the cell.

# +
# Modules from the Python Standard Library:
# "Regular Expressions" for searching and replacing text
import re

# dictionary objects with a default value
from collections import defaultdict

# stands for Input/Output
from io import BytesIO, StringIO

# specifying the locations of files on our computer
from pathlib import Path

# reading and writing zip files
from zipfile import ZipFile

# Modules from other packages:
# Matplotlib: plotting figures
import matplotlib.pyplot as plt

# BioPython: reading and writing sequence files,
# manipulating and aligning sequences, making phylogenetic trees
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import (
    NNITreeSearcher,
    ParsimonyScorer,
    ParsimonyTreeConstructor,
)
from Bio.SeqRecord import SeqRecord

# NCBI Datasets: searching and downloading NCBI data
from ncbi.datasets.openapi import GeneApi

# -

# Importing is important, but not very exciting.
# After all, there was no output when we ran the cell above.
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

# ## Troubleshooting the notebook

# ### Undoing a change

# If you get an error or would like to undo a change, select the cell and
# use `CTRL-Z` (`CMD-Z` on a mac) to undo the most recent change.

# ### Running cells in order

# Cells in a jupyter notebook can be run in any order, but they should be run from top to bottom.
# If you're getting errors, it could be because you forgot to run all the cells above the one
# you're working on.
#
# To rerun all the cells above the one you're working on:
# select the "Cell" drop-down menu and then click "Run all above"..

# ### Interrupting a cell that's taking too long

# The code in a Jupyter notebook is run by a program called the kernel.
# Most of the time, we can ignore it, but if you get stuck, it can help to know
# how to stop or restart the kernel.

# Sometimes a cell is taking too long to run, either because
# you made a mistake or because the task is bigger than you expected.
# If you'd like to stop a cell from running, you can interrupt the kernel by
# hitting the square "stop" button next to "Run" at the top of the notebook.

# ### Restarting the kernel
# Sometimes, you will want Jupyter to "forget" the
# results of the cells you've run and start fresh.
# To restart the kernel:
# 1. Click the "Kernel" drop-down menu at the top and select:
# "Restart and clear output".
# 2. Go to the place in the notebook where you left off and use the "Cell"
# drop-down menu to "Run all above" to run the previous cells and get back on track.

# ### If you would like to start over:

# If you change things in the notebook, can't get it to run, and want to start
# over:
# 1. Select the name of the notebook "workshop" at the top of the page.
# 2. Change the name to anything else like "broken-notebook".
# 3. Press the save button.
# 4. Click the original link that you followed to get to the notebook.
# 2. Go to the place in the notebook where you left off and use the "Cell"
# dropdown menu to "Run all above" to run the previous cells and get back on track.

# **Exercise**: Run the following cell

# ## Application: finding whale myoglobin orthologs

# Our scientific task today is to make a phylogenetic tree of myoglobin ortholog
# sequences in whales.
# Myoglobin is an oxygen-binding protein, similar to hemoglobin, but found primarily in muscle tissue.
# It plays an important role in allowing marine mammals to hold their breath while diving.
# We therefore expect that the myoglbin protein-coding sequence may be under particularly strong natural selection in whales.
# Finding all of the orthologs in the NCBI database, downloading their sequences, and building a tree are the first steps in an envolutionary analysis of the molecule.

# <div>
#     <img src="img/whales.png", width="300">
#     <img src="img/myoglobin.png", width="300">
# </div>

# While this task is possible without knowing any programming,
# we'll aim to show you how working with Python in a Jupyter notebook
# can make your work more reproducible and leave you with a shareable
# record of your work.
#
# We'll break the task down into several stages:
# 1. Getting the gene ID for myglobin and finding its orthologs in whales.
# 2. Downloading the ortholog transcript sequences.
# 3. Extracting the coding sequences and making an alignment.
# 4. Building and plotting a phylogenetic tree.

# ## Looking up gene IDs and an intro to data in Python

# Our first task is to find the gene ID for myglobin.
# We'll need the ID to search the database for orthologs
# and download sequence data.
#
# We are going to need three inputs:
# 1. A reference taxon
# 2. A list of gene symbols
# 3. A way to get data from NCBI over the internet
#
# We'll tackle each of these in turn.
# Along the way we'll get a lightning overview of how data is represented in Python.

# ### Data in Python: variables and objects

# Let's turn to our first input: the taxon in which we'd like to look up the gene IDs.
# For simplicity, we'll use humans.
# We may need to use this information in multiple places around our program,
# so we'll give it a name that we can refer to later.
# In programming jargon, this is called assigning a value to a _variable_.
# In Python it looks like this:

# 1.  2.  3.
taxon = "human"
#  4.   5.
print(taxon)

# Here's what's going on in that cell:
# 1. We define a "variable" called `taxon`. This is a name that give to a value
# so that we can refer to it again later.
# 2. The `=` operator means "assign the value on the right to
# the variable name on the left.
# 3. The value we're assigning to `taxon` is "human" (note the quotation marks).
# 4. The `print` function takes an input (inside the parentheses)
# and "prints" (i.e., displays) it on the screen after the cell.
# 5. We give `print` the variable `taxon` as input.
#
# The lines that begin with "#" are "comments".
# They're a way to write notes about your code that you want Python to ignore.
#
# **Exercise**: Try editing the cell to add another line of comments somewhere.
# Run the cell again. Did anything change about the output?
#
# **Exercise**: What happens if we forget the quotes in `"human"`?
# Try deleting them and running the cell. Then add them back and run it again

# In Python (and many other programming languages), we represent data with _objects_.
# In the cell above, the variable `taxon` is a name that refers to an object that stores
# a bit of data for us:
# the name of the taxon we'll use in our query.
# Every object has a _type_. To see the type of an object,
# we can use the `type` function.
# Run this cell to see the type of `taxon`:

print(type(taxon))

# The "`str`" in the output above is telling us that `taxon`
# is an object of type _string_.
# A string is Python's way of representing a sequence of text characters.
# When we assigned the value `"human"` to the variable `taxon` above,
# the quotes around human told Python that we'd like to make
# a string object with the characters `h`, `u`, `m`, `a`, and `n`.
#
# Run the cell below to check the types of a few more strings:

print(type("this is a string"))
print(type("5"))

# **Exercise**: Remove the quotes around `"5"` above and rerun the cell.
# Does the type change? What do you think is happening?

# ### Storing and using multiple values: Lists and loops

# In programming for biology, we often have more than one item of the same type
# that we'd like to perform a computation on.
# For example, we may have a list of accession numbers that we would like to
# look up in a database.
# Or we may have a list of nucleic acid sequences that we'd like to align.
#
# In Python, we represent lists of objects by surrounding them with square brackets
# and separating them by commas.
# For example:

symbols = ["MB", "BRCA1", "TP53"]
print(type(symbols))
print(symbols)

# Note that the type of a list is "list" and that when we print the list,
# it prints the objects in the list surrounded by square brackets.
#
# **Exercise**: Add another gene symbol or two to the list.
# You can look some up online, use ones you know, or choose from:
# HBB, ALB, CFTR, TNF.

# One of the most common uses of lists is to do the same action
# for every object in the list.
# In programming, this is called a _loop_ and in Python it looks like this:

# 1.  2.   3.        4.
for symbol in symbols:
    # 5.    6.                     7.
    print(symbol, type(symbol), sep="\t")

# 1. The `for` keyword tells Python that we'd like to do a loop.
# 2. We define a new variable `symbol` to represent each object in the list in turn.
# 3. `in symbols` tells Python which list we'd like to loop over.
# 4. The `:` says that we're begining the _body_ of the loop. Everything after the colon
# at the same level of indentation will be executed
# one time for each object in `symbols`.
# 5. The body of this loop is just a single call of the `print` function.
# 6. We can use the loop variable `symbol` in the body of the loop
# to refer to an object in `symbols`.
# 7. The `print` function can take an optional _keyword argument_ `sep`
# (short for separator),
# which is a string to print between the arguments (here the tab character `"\t"`).
#
# **Question**: Why does the loop print "<class 'str'>" every on iteration?
#
# **Exercise**:
# Modify the cell above to print the length of each symbol instead of its type.
# Hint: Replace the `type` function with `len` (length).

# #### Making lists from other lists

# Very often, we have a list of objects and we'd like to make a new list
# by doing the same thing to each object in the list.
#
# For example, say we had our list of gene symbols from humans, which use all caps,
# and we turn them into gene symbols for mice, which only capitalize the first letter.
# To convert between the conventions,
# we can use the `capitalize` method of our string objects.
# (We'll talk about methods in the next section.)
# It works like this:

human_shh = "SHH"
mouse_shh = human_shh.capitalize()
print(mouse_shh)

# This is fine for a single gene name, but we'll want to convert the whole list.
# To do this, we can use a _list comprehension_:

human_symbols = symbols
#                      1.              2.         3.
mouse_symbols = [sym.capitalize() for sym in human_symbols]
print(mouse_symbols)

# List comprehensions look a lot like a loop wrapped in a list:
# 1. We call the `capitalize` method on every symbol.
# 2. `sym` is our loop variable. It refers to each object in the original list in turn.
# 3. `human_symbols` is the original list we're looping over.

# With lists and loops we're starting to see the real power of programming:
# we can store large amounts of data and automate repetitive tasks on it.

# ### Getting gene IDs from NCBI

# Now that we know now to represent data as objects and
# store multiple items and perform repeated computations with lists and loops,
# it's time to get some actual data.
#
# In particular, we'd like to find the gene ids for our list of gene symbols,
# and store them in a dictionary where we can look them up later.
#
# We're going to use the NCBI Datasets Python package, which will allow us
# to download data and metadata from NCBI from within our notebook.
#
# The first step is to create a special object called `GeneApi`.
# An API (Application Programming Interface) is a set of rules
# that define how two programs can talk to one another.
# In our case, our Python program in this notebook will be requesting data
# from the servers at NCBI.
# `GeneApi` is responsible for making those requests for us.
#
# Run the next cell to make a `GeneApi` object:

gene_api = GeneApi()
print(type(gene_api))

# That's a much longer type signature than "str" or "list"!
# The chain of words connected by dots is telling us
# where the GeneApi class of objects is defined.
# You don't need to worry about the details, but recall that
# in the second code cell of the notebook we had this import statement:
#
# ```python
# from ncbi.datasets.openapi import GeneApi, GenomeApi
# ```
#
# This line of code told our program to look in the `ncbi.datasets` package and
# import the objects `GeneAPI` and `GenomeAPI` so we could use them here.

# #### Downloading metadata from NCBI
#
# Now that we have our `GeneApi` object, we can use it to get some data.
# To do so, we need to learn a new aspect of objects in Python: _attributes_.
# Attributes are things that "belong" to an object.
# They are accessed by putting a `.` after the object's name.
#
# You can find detailed information on how to use the NCBI Datasets
# Python package here: https://www.ncbi.nlm.nih.gov/datasets/docs/languages/python/
#
# Some attributes, called _methods_ are also a function that do something.
# Here we will use one of `GeneApi`'s methods to get some metadata for our list of gene symbols:

#                  1.            2.                       3.
metadata = gene_api.gene_metadata_by_tax_and_symbol(symbols, taxon)

# 1. The `.` after `gene_api` means that we are going to access one of the object's _attributes_.
# 2. `gene_metadata_by_tax_and_symbol` is a _method_ attribute.
# That is, it's a function that "belongs" to the `gene_api` object.
# The purpose of the method is in it's name: it looks up gene metadata by taxon and gene symbol.
# 3. This method takes two arguments: a list of gene symbols, and a taxon.

# Let's see what sort of metadata we got:

print(metadata)

# Our `metadata` object is a structured record that holds lots of information
# about the genes we looked up. Scroll through and take a look.

# #### Unpacking information from the metadata object

# Our next task is to "unpack" the structure of the `metadata` object
# and pull out the information that we want.
# We will do this by accessing its attributes.
# At the top level, there is an attribute called `genes`.
# We can access it with `metadata.genes`.
#
# Look at the output above and take a guess what type of object `metadata.genes` is.
# (Hint: Notice the square bracket in the first line of the output.)
# Check your guess by running the next cell.

print(type(metadata.genes))

# We've seen that type before!
#
# Just like we did in the section on lists, we can use a loop do something to each
# element of the list `metadata.genes`.
# Let's print the "description" attribute of each gene:

for gene_data in metadata.genes:
    print(gene_data.gene.description)

# **Exercise**: Scroll up and take a look at the output of `print(metadata)`.
# Pick an attribute you'd like to inspect for each gene and modify the previous cell
# to print it for you.

# We're now in a position to extract the information we want:
# pairs of gene symbols and gene IDs:

for gene_data in metadata.genes:
    symbol = gene_data.gene.symbol
    gene_id = gene_data.gene.gene_id
    #     1.           2.             2.   3.
    print(f"Symbol: {symbol}\tId: {gene_id:>4}")

# **Aside**: We're printing a special formatted string or "f-string",
# which lets us insert the values of objects into a string and format how they are displayed.
# 1. An "f-string" is marked by an `f` before the quotation marks.
# 2. We can use curly braces to insert our variables `symbol` and `gene_id` into the f-string.
# 3. A `:` after a value begins a formatting specification.
# `>4` means right-aligned text that's four character's long.

# #### Looking up values: dictionaries

# Now we have our list of gene symbols and we can print a table of
# symbol-id pairs for each gene.
# This table is nice for humans to look at, but it's not the best format
# for our program to be able access the data.
# What we'd really like do do is be able to take a gene symbol and look up
# the associated gene id.
#
# Python has a type that lets us do exactly that: the `dict` or "dictionary".
# (Note: lists, dicts, and our metadata object are all examples of _data structures_,
# ways of storing data in a logical, structured format.)
#
# A dict is a set of key-value pairs.
# Each key is associated with only one value,
# and we can look up the value associated with each key.
#
# Here's an example of a dict where the keys are letters and the values are numbers:

#              1.  2.
example_dict = {"a": 0, "b": 1, "c": 1}
print(example_dict)
#                   3.
print(example_dict["a"])

# 1. Dicts are created like lists but with curly instead of square brackets.
# 2. Key-value pairs are joined together with colons.
# 3. We can look up the value associated with a key using square brackets.
#
# **Exercise**: Change the last line to look up another value.
#
# **Exercise**: Add another key-value pair to the dict.
#
# **Exercise**: What happens when you try to look up a value that's not in the dict?

# Just like we did with lists, we can create a dict from a list using a comprehension:

#                    key:      value           loop        list
mouse_symbol_dict = {sym: sym.capitalize() for sym in human_symbols}
print(mouse_symbol_dict)
print(f"The mouse symbol for MB is: {mouse_symbol_dict['MB']}")

# Now we can apply the same logic to making a dict that lets us look up
# the gene id for each symbol:

gene_id_dict = {
    #        1.             3.           2.
    gene_data.gene.symbol: int(gene_data.gene.gene_id)
    for gene_data in metadata.genes
}
print(gene_id_dict)


# 1. Our keys will be gene symbols.
# 2. Our values will be gene IDs.
# 3. We use the `int` function to convert the gene IDs from strings to integers.
#
# **Exercise**: In the cell below, use our new dict to look up the gene id myoglobin (MB).


# #### Making it reusable: defining functions

# Let's review what we've done so far:
# 1. We made a _list_ of gene symbols.
# 2. We used NCBI Datasets' `GeneAPI` to look up metadata for those genes.
# 3. We extracted the gene ids from the that metadata.
# 4. We've created a dictionary that associates gene symbols with their IDs.
#
# The only problem is that our code to do this is scattered all over the notebook,
# interspersed with various exercises, print statements, etc.
# If we wanted to change something: add some genes, or change the taxon,
# we'd have to track down all the different pieces and cross our fingers that we
# didn't miss something important.
#
# The practice of programming is all about writing code to solve a particular problem and
# then modifying it to be more re-usable.
# A very useful technique to that end is identifying a common task
# and packaging the code necessary to do that task into a _function_.
#
# A function takes one or more *arguments* as input, does some computation,
# and *returns* some output.
# You use a function by "calling" it with the appropriate arguments.
# We've seen a few examples of functions: `type`, `len`, `print`, and the methods `String.capitalize` and `GeneApi.gene_metadata_by_tax_and_symbol`.
#
# The benefits of writing functions are:
# 1. It simplifies repeated calculations.
# Instead of copying all of the steps individually,
# we just call the function with different arguments.
# 2. It makes it easier to make changes.
# We only have to modify the function definition in one place
# and the changes will propagate throughout our program.
# 3. It makes our code easier to read.
# It's not necessary to understand all the inner workings of a function,
# just what its inputs and outputs are.
#
# Here's a single function that takes a list of symbols and a taxon and
# returns a dictionary of (gene symbol : gene id) pairs.

# +
# 1.      2.        3.      4.
def get_gene_ids(symbols, taxon):
    # 5.
    gene_api = GeneApi()
    #                                                              6.
    gene_metadata = gene_api.gene_metadata_by_tax_and_symbol(symbols, taxon)
    gene_id_dict = {
        gene_data.gene.symbol: int(gene_data.gene.gene_id)
        for gene_data in gene_metadata.genes
    }
    #     7.
    return gene_id_dict


gene_symbols = ["MB", "BRCA1", "TP53"]
#                 8.               9.
gene_ids = get_gene_ids(gene_symbols, "human")
#        10.
print(gene_ids)
# -

# Here's what's going on in that cell:
# 1. The `def` keyword begins the definition of a function.
# 2. The name of the function is `get_gene_ids`.
# 3. The function takes two argument: a list, `symbols` and a string, `taxon`.
# The names of the arguments work like variables
# that are only visible inside the function definition.
# 4. The function body begins with `:`, just like with a `for` loop.
# 5. The _body_ of the function is the code that runs whenever we call a function.
# It is indented like a loop.
# The code in the body here should all be familiar to you from above.
# 6. Inside the body, we can refer to the arguments by name.
# 7. The function `return`s the value of `gene_id_dict`.
# 8. Later, we can "call" our function by using it's name: `get_gene_ids`.
# 9. We give `get_gen_ids` the list `gene_symbols` and the string `"human"` as arguments.
# 10. We store the return value (the dict of gene ids) in a variable `gene_ids`.
#
# When we call `get_gene_ids` with the first argument `gene_symbols`,
# Python substitutes the value of `gene_symbols` for the argument `symbols`
# everywhere in the function body.
# `get_gene_ids` calls the API, gets the metadata, extracts the ids, and returns the dict.
#
# **Exercise**: In the cell below, get the gene ids for a different set of gene symbols.

# FIX ME     "{}" is an empty dictionary.
my_gene_ids = {}
print(my_gene_ids)

# **Congratulations!** You've just completed a crash course in Python programming.
#
# In the interest of time, we're going to stop explaining everything in detail
# as we proceed through the rest our analysis of the whale myoglobin data.
# Our focus is going to switch from showing you how Python works,
# to showing you what you can get done with Python, BioPython, and NCBI Datasets.
#
# We will encounter some complications that are beyond the scope of this tutorial.
# We will package thes complications inside of functions.
# You do not need to understand everything that is going on inside the
# function definitions.
# (But you're welcome to look!)

# ## Finding whale myoglobin orthologs


# Now that we have the gene ID for myoglobin,
# we can use the NCBI Datasets package to search for
# all of the whale myoglobin orthologs.
#
# Here's a function that takes a gene id and a list of taxa and gets
# metadata on the orthologs for that gene in the taxa:


def query_orthologs(gene_id, taxa):
    gene_api = GeneApi()
    response = gene_api.gene_orthologs_by_id(gene_id, taxon_filter=taxa)
    return [item.gene for item in response.genes.genes]


# The components of the function should look familiar to you.
#
# **Questions**:
# 1. What are the _arguments_ to `query_orthologs`?
# 2. What is the name of the _method_ of `GeneApi` that we call in the function body?
# 3. Can you spot the _list comprehension_?
# 4. What is the type of the return value?

# Now lets use our function to get the myglobin orthologs for
# whales and two primate outgroups: humans and rhesus monkeys:

taxa = ["whales", "human", "Macaca mulatta"]
mb_orthologs = query_orthologs(gene_ids["MB"], taxa)
print(mb_orthologs)


# How many orthologs did we find?

print(len(mb_orthologs))

# Just as we did with the gene metadata, we can loop over our list
# of orthologs and extract their attributes using `.`.
# Here are the common and scientific names of the species where we found orthologs:

for g in mb_orthologs:
    print(f"{g.common_name}: {g.taxname}")

# Finally, we will need to get the gene ID for each ortholog
# in order to download their sequences.

ortholog_gene_ids = [int(g.gene_id) for g in mb_orthologs]
print(ortholog_gene_ids)

# By now, this pattern should look familiar:
# 1. Get some data, store it in a data structure like a list or a dict.
# 2. Perform some computations on each piece of data,
# e.g., extracting an attribute or calling a function.
# 3. Store the results in a new data structure.
#
# **Question**: Why are we using the `int` function?
# (Hint, you've seen it before with gene ids.)

# ## Downloading gene sequences


# Armed with our list of ortholog gene ids,
# we're ready to download the gene sequences themselves.
#
# Let's define another function to do that.
# This function takes two arguments:
# 1. The list of gene ids.
# 2. The directory where we'd like to put the data.
#
# It downloads the transcript sequences from NCBI
# and writes them to a fasta file in the data directory.
#
# **Note**: "fasta" is a common data format for storing sequence data.
#
# **Note**: This function is more complicated than the ones you've seen so far.
# Don't worry if you don't understand everything in the definition.


def download_transcripts(gene_ids, data_directory):
    # Create a GeneApi object
    gene_api = GeneApi()
    print("Begin download of data package ...")
    # Use the GeneApi to download a fasta file of gene transcripts.
    gene_ds_download = gene_api.download_gene_package(
        gene_ids, include_annotation_type=["FASTA_RNA"], _preload_content=False
    )
    # The downloaded data package is formatted as a zip file.
    # Extract the fasta file and write it to your hard drive.
    with ZipFile(BytesIO(gene_ds_download.data)) as zipfile:
        data_file_name = zipfile.extract(
            "ncbi_dataset/data/rna.fna", path=data_directory
        )
    print(f"Download completed -- see {data_file_name}")
    # Return the path to the fasta file you just downloaded.
    return Path(data_file_name)


# Let's see it in action:

# We're going to put all our data in a directory called `data`
data_dir = Path("../data")
# Download and save the name of the fasta file.
ortholog_fasta = download_transcripts(ortholog_gene_ids, data_dir)

# We can make sure something happened by checking the contents of the data directory.

print("Data directory contents after download:\n")
# ! ls -R {data_dir}

# The last line of the output shows that we did in fact
# download a fasta file called `rna.fna`
#
# **Note**: The `!` at the beginning of a line is a Jupyter trick to
# run a shell command.
# If you don't know what that means, don't worry, just know that:
# 1. `ls -R` is a command that gets the contents of a directory recursively.
# 2. This won't work in a normal Python program without Jupyter.

# ## Reading and manipulating sequence records

# Now that we have our ortholog transcript sequences downloaded,
# we can import them into our program and start working with them.
#
# We need to get one coding sequence for each species so that we
# can make our multiple sequence alignment and then our tree.
# In order to do that, we will:
# 1. Read the transcript sequences from the fasta file.
# 2. Identify the longest transcript for each species from among the splice variants.
# 3. Extract the coding sequence.
#
# To read the list of sequences, we use the function `SeqIO.parse` from
# the BioPython package.
# This function can read a variety of file formats, so we have to specify
# that our file is in "fasta" format:

#                  Get a list of sequences from the fasta file
ortholog_records = list(SeqIO.parse(ortholog_fasta, "fasta"))
# Print each sequence record
for record in ortholog_records:
    print(record, "\n")


# Our next step is to organize our list of sequences by organism.
# If you look at the output above, you'll notice that "organism"
# is specified in the `description` attribute of each sequence record.
#
# To the organism from a record, we need to:
# 1. Access the `description` attribute of the record using `.description`.
# 2. Match the pattern `"[organism=x]"` and get the value at `"x"`.
#
# We already know how to do the first part.
# For the second we will use "Regular Expressions" (re).
# Regular expressions are a very powerful tool for matching patterns in text.
# We won't explain how they work today.
#
# Here's a function to get the organism name from of a record:

# +
def get_organism_name(record):
    # Use re.search to match the pattern and capture the result.
    match = re.search(r"\[organism=(?P<name>[\w ]+)\]", record.description)
    # If there's a match, return the organism name.
    if match:
        return match.group("name")
    # If not, indicate that no name was found.
    else:
        return "OrganismNameNotFound"


for rec in ortholog_records:
    print(get_organism_name(rec))


# -

# Clearly, we have multiple transcript records per organism.
# This is because of splice variants.
#
# Let's organize our list of records into a dict where we
# can look up all the records for a given organism.
#
# The following function groups records by organism,
# building a dictionary with organism names as keys
# and lists of records as values.


def records_by_organism(records):
    # defaultdict is like a dict but with a default value for missing keys
    # here the default is an empty list, i.e., no records
    organism_dict = defaultdict(list)
    for record in records:
        org = get_organism_name(record)
        # Add the record to the list for its organism
        organism_dict[org].append(record)
    return organism_dict


# We can apply the function and use `len` to get the number of records
# in each organism's record list:

organism_records = records_by_organism(ortholog_records)
# the items method gives a list of key-value pairs
# we can loop over them simultaneously
for org, records in organism_records.items():
    print(f"{org}: {len(records)}")

# This is nice, but we might rather look at the data as a histogram:

counts = [len(records) for records in organism_records.values()]
plt.hist(counts, bins=range(1, 10))
plt.xlabel("Number of records")
plt.ylabel("Number of species")
plt.show()

# We can also make a histogram of the lengths of the transcripts:

transcript_lengths = [
    len(record.seq) for records in organism_records.values() for record in records
]
plt.hist(transcript_lengths)
plt.xlabel("Transcript length")
plt.ylabel("Number of transcripts")
plt.show()


# Now to get the longest transcript for each organism.
# We can use the `max` function to get the largest value in a list.
# `max` takes a keyword argument `key`, which is the function to use
# to determine the size of an object in the list.
# We define a function `record_length` that gets the length of
# a record's sequence, and use it as `key`.
#
# We store the result in a dictionary that maps an organism to the record
# of its longest transcript:

# +
def record_length(rec):
    return len(rec.seq)


def get_longest_transcripts(record_dict):
    return {
        org: max(records, key=record_length)
        for org, records in organism_records.items()
    }


longest_transcripts = get_longest_transcripts(organism_records)
for org, rec in longest_transcripts.items():
    print(org, rec.id)
# -


# ### Extracting the coding sequence


# Our next step is to get the coding sequence (CDS) from each transcript.
# (The transcript includes non-coding RNA upstream of the start codon
# and downstream of the stop codon.)
#
# To find where the CDS starts and stops, we consult our
# gene metadata:

print(mb_orthologs)

# **Exercise**: Look through the output above and try to find the
# start and stop coordinates for the CDS.

# Here are some functions for extracting the CDS start and stop
# coordinates from the metadata.
# We return a dictionary that maps accession numbers to CDS start and stop coordinates.
#
# This is mostly familiar: comprensions, and accessing object attributes.
# There are a few wrinkles, noted in the comments.

# +


def cds_region(transcript):
    # The range object is a list to account for the possibility
    # of multiple ranges. `[0]` takes the first one.
    cds_range = transcript.cds.range[0]
    # `(a, b)` is a pair of values. We're returning two things: begin and end
    return (int(cds_range.begin), int(cds_range.end))


def get_cds_regions(gene_list):
    return {
        transcript.accession_version: cds_region(transcript)
        # Using multiple `for`s loops over nested lists:
        for gene in gene_list
        for transcript in gene.transcripts
    }


cds_regions = get_cds_regions(mb_orthologs)
print(cds_regions)
# -

# With start and end positions in hand, we can can create new records that
# only hold the CDS region, rather than the whole transcript sequence.
# For each organism:
# 1. We look up the CDS region by its accession.
# 2. Create a new SeqRecord object.
# 3. Give the new object the organism as an ID.
# (This will make our eventual tree more interpretable).
# Note that we have to remove the spaces in the organism name to have a valid ID.
# 4. Take the coding sequence. The `:` is "slice" notation for taking part of a sequence
# from one position to another. The `- 1` is because of differences between how
# Python indexes a sequence and how NCBI does.


# +


def get_cds_records(transcript_dict, cds_regions):
    cds_records = []
    for organism, record in transcript_dict.items():
        # 1.
        start, end = cds_regions[record.id]
        #                2.
        cds_record = SeqRecord(
            #                 3.
            id=organism.replace(" ", "_"),
            name=record.name,
            description=record.description,
            #                     4.
            seq=record.seq[start - 1 : end],
        )
        cds_records.append(cds_record)
    return cds_records


cds_records = get_cds_records(longest_transcripts, cds_regions)
# -


# We should check to make sure that we haven't made any mistakes.
# For example, what if we got the tricky indexing in 4. above wrong?
# We can use the `translate` method to get the amino acid sequence
# corresponding to our transcript.
# This lets us check three crucial features of the CDS:
# 1. Its length should be a multiple of 3. If it's not, `translate` will complain.
# 2. It should start with a start codon (translated as "M").
# 3. It should end with a stop codon (translated as "*").
#
# Let's check our work:

for rec in cds_records:
    print(rec.seq.translate())

# This all looks good.
# But we're programming, and the point of programming is to automate things.
# We can write some functions to check that our translations are good.
# Here, we use `startswith` and `endswith`, which return `True`
# if a sequence starts or ends with the argument and `False` otherwise.
# We also use the function `all`, which returns `True` if all its arguments are `True`.

# +


def check_start_codons(proteins):
    return all([p.startswith("M") for p in proteins])


def check_stop_codons(proteins):
    return all([p.endswith("*") for p in proteins])


protein_seqs = [rec.seq.translate() for rec in cds_records]
print(check_start_codons(protein_seqs))
print(check_stop_codons(protein_seqs))
# -

# ## Creating a multiple sequence alignment and building a tree


# BioPython lets us call the multiple sequence alignment tool MUSCLE
# from within our Python program.
#
# MUSCLE requires that our sequences be written to a file,
# so we'll make a fasta file of our CDS sequences
# using `SeqIO.write`.

cds_fasta = data_dir / "cds.fna"
SeqIO.write(cds_records, cds_fasta, "fasta")
# Check that we created a file `cds.fna` in our data directory
# ! ls {data_dir}

# Now we can call MUSCLE.
# (Note that we have to provide the path to our local installation of MUSCLE.)

# +


def align_with_muscle(input_fasta):
    muscle_exe = Path("../bin/muscle3.8.31_i86linux64")
    muscle_cline = MuscleCommandline(muscle_exe, input=input_fasta)
    # The variable `stdout` ("standard out") captures the output from MUSCLE
    # `stderr` ("standard error") captures any errors.
    stdout, stderr = muscle_cline()
    # `AlignIO` reads an alignment
    # `StringIO` lets BioPython treat a string as though it were a file
    return AlignIO.read(StringIO(stdout), "fasta")


alignment = align_with_muscle(cds_fasta)
print(alignment)


# -


# We have our alignment!
#
# We're finally ready to make a tree using functions from BioPython's `Phylo` module.
# We'll make a maximum parsimony tree.
# In the cell below, `scorer` and `searcher` are helper objects that tell
# `ParsimonyTreeConstructor` how to build the tree.

# +
def build_parsimony_tree(alignment):
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor = ParsimonyTreeConstructor(searcher)
    return constructor.build_tree(alignment)


tree = build_parsimony_tree(alignment)
print(tree)
# -


# Printing the tree gives us a peak at how Python represents the tree data structure.
# But it's not very nice to look at.
#
# `draw_ascii` gives a simple asii art rendering of the tree.
# It's good enough for our purposes:

Phylo.draw_ascii(tree)


# We could make a prettier tree with some of BioPython's more complicated tree plotting
# functions, but we don't have time today.
# To learn more, see the
# [Phylo module's documentation](https://biopython.org/wiki/Phylo).


# ## Putting it all together

# It's time to take a step back and appreciate what we've done.
# The following cell contains our whole program from start to finish.
# That is, from a gene symbol and a list of taxa
# to a tree build from ortholog coding sequences.
#
# Because we took the time to package our code into functions that
# perform discrete tasks, we're left with a compact
# and fairly readable block of code.
# Moreover, if we wanted to repeat this analysis with a different
# gene or set of taxa, it's clear what we'd need to change.

# +
# Inputs:
gene_symbol = "MB"
reference_taxon = "human"
search_taxa = ["whales", "human", "Macaca mulatta"]
data_dir = Path("../data")

# Getting orthologs:
gene_ids = get_gene_ids([gene_symbol], reference_taxon)
orthologs = query_orthologs(gene_ids[gene_symbol], search_taxa)
ortholog_gene_ids = [int(g.gene_id) for g in orthologs]
ortholog_fasta = download_transcripts(ortholog_gene_ids, data_dir)
ortholog_records = list(SeqIO.parse(ortholog_fasta, "fasta"))
organism_records = records_by_organism(ortholog_records)

# Finding the longest transcript and extracting the CDS
longest_transcripts = get_longest_transcripts(organism_records)
cds_regions = get_cds_regions(orthologs)
cds_records = get_cds_records(longest_transcripts, cds_regions)

# Making an alignment:
cds_fasta = data_dir / "cds.fna"
SeqIO.write(cds_records, cds_fasta, "fasta")
alignment = align_with_muscle(cds_fasta)

# Building and plotting the tree:
tree = build_parsimony_tree(alignment)
Phylo.draw_ascii(tree)
# -

# ## Where to go now
#
# If we've piqued your interest about learning to program for
# biology applications, here are a few resources to keep learning:
#
# ### [Run this notebook online](https://mybinder.org/v2/gh/drice-codeathons/workshop-ncbi-data-with-python/main?filepath=notebooks%2Fworkshop.py)
#
# If you'd like to explore this notebook on your own,
# you can use the free service MyBinder to run the notebook remotely
# by following [this link](https://mybinder.org/v2/gh/drice-codeathons/workshop-ncbi-data-with-python/main?filepath=notebooks%2Fworkshop.py).
#
# ### Installing Jupyter on your computer
#
# If you're ready to take the plunge and install Jupyter on your
# own computer, there are installation instructions on the
# [Jupyter project website](https://jupyter.org).
#
# ### NCBI Datasets
#
# To learn more about the NCBI Datasets tools, you can visit the
# [Datasets homepage](https://www.ncbi.nlm.nih.gov/datasets/).
# There is lots of information about the web interface,
# command line tools, and Python and R libraries.
#
# ### The BioPython package
#
# To learn what else you can do with BioPython, you can see their
# [documentation here](https://biopython.org).
# Particularly useful are their
# [tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
# and [cookbook](https://biopython.org/wiki/Category%3ACookbook).
#
# ### Learning Python
#
# There are lots of resources online to help you learn Python in more depth.
# One good place to start is the
# [Python beginner's guide](https://www.python.org/about/gettingstarted/).
