[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/drice-codeathons/workshop-ncbi-data-with-python/main?filepath=notebooks%2Fworkshop.py)

# Introduction to NCBI Data with Python
Created by Daniel Rice and Peter Cooper, Spring 2021

## Learning objectives:
- Objective 1
- Objective 2
- Objective 3

## Materials:

We will be using [The Littlest JupyterHub](https://tljh.jupyter.org/en/latest/index.html) to serve Jupyter notebooks to a class of 30--50 students.

Resource usage:
- During BLAST memory usage spikes to 500 MB, stabilizes at 215 MB.

## Instructions:

1. Start AWS instance with TLJH. Copy the `jupyterhub_setup.sh` from this repo
into the instances "User data" field.
2. Connect the instance to the ElasticIP
3. Wait ~10 min for the Hub to be set up.
4. Install SSL certificate and configure hub to use it
5. Log in by going to:
https://jupyterhub01.ncbi.nlm.nih.gov/hub?next=%2Fuser-redirect%2Fgit-pull?repo%3Dhttps%253A%252F%252Fgithub.com%252Fdrice-codeathons%252Fworkshop-ncbi-data-with-python%26branch%3Dmain%26urlpath%3Dtree%252Fworkshop-ncbi-data-with-python%252F
6. Configure jupyterhub by running `./config.sh`
7. Add the students to the user accounts
8. Send students to:
https://jupyterhub01.ncbi.nlm.nih.gov/hub?next=%2Fuser-redirect%2Fgit-pull?repo%3Dhttps%253A%252F%252Fgithub.com%252Fdrice-codeathons%252Fworkshop-ncbi-data-with-python%26branch%3Dmain%26urlpath%3Dtree%252Fworkshop-ncbi-data-with-python%252Fworkshop.py

`jupytext` will automatically create the notebook from `workshop.py`
