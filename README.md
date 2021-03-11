# Introduction to NCBI Data with Python
Created by Daniel Rice and Peter Cooper, Spring 2021

## Learning objectives:
- Objective 1
- Objective 2
- Objective 3

## Materials:

We will be using [The Littlest JupyterHub](https://tljh.jupyter.org/en/latest/index.html) to serve Jupyter notebooks to a class of 30--50 students.

## Instructions:

1. Start AWS instance with TLJH.
2. Wait ~10 min for the Hub to be set up.
3. Get public IP address for AWS instance.
4. Go to:
http://<INSTANCE.IP.ADDRESS>/hub?next=%2Fuser-redirect%2Fgit-pull?repo%3Dhttps%253A%252F%252Fgithub.com%252Fdrice-codeathons%252Fworkshop-ncbi-data-with-python%26branch%3Dmain%26urlpath%3Dtree%252Fworkshop-ncbi-data-with-python%252F
5. Log in as admin, open a console and run `sudo -E pip install -r requirements.txt`
6. Send students to:
http://<INSTANCE.IP.ADDRESS>/hub?next=%2Fuser-redirect%2Fgit-pull?repo%3Dhttps%253A%252F%252Fgithub.com%252Fdrice-codeathons%252Fworkshop-ncbi-data-with-python%26branch%3Dmain%26urlpath%3Dtree%252Fworkshop-ncbi-data-with-python%252Fnotebooks%252F<NOTEBOOK>.py

`jupytext` will automatically create the notebook from `<NOTEBOOK>.py`
