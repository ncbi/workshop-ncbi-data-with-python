#!/bin/bash
# This file should be in the "User data" of the EC2 instance
curl -L https://tljh.jupyter.org/bootstrap.py \
  | sudo python3 - \
    --admin ricedap \
    --user-requirements-txt-url https://raw.githubusercontent.com/drice-codeathons/workshop-ncbi-data-with-python/main/requirements.txt
