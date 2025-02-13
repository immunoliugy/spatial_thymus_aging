#!/usr/bin/env bash
set -ex

# Function to run notebooks in separate processes
run_notebook() {
    local notebook_name="$1"
    
    # Print start time
    echo "Starting $notebook_name at $(date)"
    local start_time=$(date +%s)
    
    python -m jupyter nbconvert --to html --execute --ExecutePreprocessor.allow_errors=True --ExecutePreprocessor.timeout=-1 --FilesWriter.build_directory=../results "$notebook_name"
    
    # Print end time
    local end_time=$(date +%s)
    echo "Finished $notebook_name at $(date)"
    
    # Calculate and print duration
    local duration=$((end_time - start_time))
    echo "Execution time for $notebook_name: ${duration} seconds"
    
    # Kill any lingering Jupyter processes to free memory
    pkill -f jupyter || true
    sleep 5  # Short wait before the next execution
}

# WORKS
run_notebook Reproduce_Fig_1.ipynb # 89 seconds
run_notebook Reproduce_Fig_2.ipynb # 112 seconds
run_notebook Reproduce_Fig_3.ipynb # 536 seconds
run_notebook Reproduce_Fig_4.ipynb # 96 seconds
run_notebook Reproduce_Fig_4_pt2.ipynb # 443 seconds
run_notebook Reproduce_Fig_4_pt3.ipynb # 1114 seconds
run_notebook Reproduce_Fig_5.ipynb # 141 seconds
run_notebook Reproduce_Fig_5_pt2.ipynb # 89 seconds 
run_notebook Reproduce_Fig_5_pt3.ipynb # 43 seconds
run_notebook Reproduce_Fig_5_pt4.ipynb # 75 seconds 

# # Uncomment if needed
# run_notebook Reproduce_Fig_S1.ipynb # 469 seconds
# run_notebook Reproduce_Fig_S2.ipynb # 107 seconds
# run_notebook Reproduce_Fig_S2_pt2.ipynb # 124 seconds
# run_notebook Reproduce_Fig_S2_pt3.ipynb # 76 seconds
# run_notebook Reproduce_Fig_S3.ipynb # 77 seconds
# run_notebook Reproduce_Fig_S3_pt2.ipynb # 106 seconds
# run_notebook Reproduce_Fig_S4.ipynb # 82 seconds
# run_notebook Reproduce_Fig_S5.ipynb # 260 seconds 
# run_notebook Reproduce_Fig_S6.ipynb # 57 seconds
# run_notebook Reproduce_Fig_S8.ipynb # 52 seconds
# run_notebook Reproduce_Fig_S9_pt1.ipynb # 152 seconds 
# run_notebook Reproduce_Fig_S9_pt2.ipynb # 121 seconds
# run_notebook Reproduce_Fig_S9_pt3.ipynb # 113 seconds
# run_notebook Reproduce_Fig_S10.ipynb # 324 seconds
