#!/bin/bash

# Turn off resource usage monitor
sudo -E jupyter nbextension disable jupyter_resource_usage/main --sys-prefix

# Set timeout to 4hrs
sudo tljh-config set services.cull.timeout 14400

# Set server limits
sudo tljh-config set limits.memory 8G
sudo tljh-config set limits.cpu 2

sudo tljh-config reload
