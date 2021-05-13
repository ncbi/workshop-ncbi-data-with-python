[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/drice-codeathons/workshop-ncbi-data-with-python/main?filepath=notebooks%2Fworkshop.py)

# Introduction to NCBI Data with Python
Created by Daniel Rice and Peter Cooper, Spring 2021

## Learning objectives:
In this workshop you will learn how to:
- Use Python programming to download, analyze, and visualize data.
- Use Jupyter to create data analysis 'lab notebooks' that make it easy to reproduce
and share what you did.
- Find data that is relevant to your project using the new NCBI Datasets resource.
- Download NCBI sequence data and manipulate it with the BioPython package.

## Materials:

We will be using [The Littlest JupyterHub](https://tljh.jupyter.org/en/latest/index.html) to serve Jupyter notebooks to a class of 30--50 students.

Resource usage:
- During BLAST memory usage spikes to 500 MB, stabilizes at 215 MB, as measured by Jupyter's usage log extension.
- The downloaded data takes up only 65K of storage.
- The repository itself is about 40MB, mostly in the `bin` directory that holds the blast and muscle binaries.

## Instructions:

1. Start AWS instance with TLJH. Copy the `tljh-setup.sh` from this repo
into the instances "User data" field.
2. Connect the instance to the ElasticIP
3. Wait ~10 min for the Hub to be set up.
4. Install SSL certificate and configure hub to use it
5. Log in by going to:
https://jupyterhub01.ncbi.nlm.nih.gov/hub?next=%2Fuser-redirect%2Fgit-pull?repo%3Dhttps%253A%252F%252Fgithub.com%252Fdrice-codeathons%252Fworkshop-ncbi-data-with-python%26branch%3Dmain%26urlpath%3Dtree%252Fworkshop-ncbi-data-with-python%252F
6. Configure jupyterhub by opening a terminal and running `tljh-config.sh`.
   Restart your server if it's running.
7. Add the students to the user accounts under: Control panel -> Admin -> Add users.
8. Send students to:
https://jupyterhub01.ncbi.nlm.nih.gov/hub?next=%2Fuser-redirect%2Fgit-pull?repo%3Dhttps%253A%252F%252Fgithub.com%252Fdrice-codeathons%252Fworkshop-ncbi-data-with-python%26branch%3Dmain%26urlpath%3Dtree%252Fworkshop-ncbi-data-with-python%252Fnotebooks%252Fworkshop.py

`jupytext` will automatically create the notebook from `workshop.py`

### Recovering a broken notebook
If a student changes their notebook so it won't run and can't undo the changes:
- Rename the notebook by clicking on the title "workshop" at the very top of the page.
  It doesn't matter what name they choose.
- Click the workshop link again to download a fresh copy.
