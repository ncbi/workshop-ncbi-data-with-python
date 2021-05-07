#!/bin/bash

# Set timeout to 4hrs
sudo tljh-config set services.cull.timeout 14400
sudo tljh-config reload

# Turn off resource usage monitor
sudo -E jupyter nbextension disable jupyter_resource_usage/main --sys-prefix
