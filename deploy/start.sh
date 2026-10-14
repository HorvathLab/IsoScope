#!/bin/bash
set -e
cd /code
python3 -c "
import os
from huggingface_hub import snapshot_download
snapshot_download('HorvathLab/IsoScopeDataset', repo_type='dataset', local_dir='/code/shiny_data', token=os.environ['IsoScope'])
"
R --quiet -e "shiny::runApp(host='0.0.0.0', port=7860)"
