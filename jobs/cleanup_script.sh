#!/bin/bash

###########################
# USAGE WITH RUN_JOB.ENV:
# source run_cleanup.env <home_directory_path> # NOT WORKING YET....#
# WHICH RUNS THIS SCRIPT:
# ./cleanup_script.sh <home_directory_path>
###########################

## Check if HOME_DIRECTORY_PATH environment variable is set
#if [ -z "$HOME_DIRECTORY_PATH" ]; then
#    echo "HOME_DIRECTORY_PATH is not set. Make sure to run this script through run_cleanup.env."
#    exit 1
#fi

#!/bin/bash

home_dir_path="/data/xchem-fragalysis/kfieseler/D68EV3CPROA/elabs/2_step"
suffixes_to_keep=( '.minimised.json' '.holo_minimised.pdb' '.minimised.mol' )

# Function to check if a file has a suffix to keep
function has_suffix_to_keep {
    local file="$1"
    for suffix in "${suffixes_to_keep[@]}"; do
        if [[ "$file" == *"$suffix" ]]; then
            return 0  # Found a suffix to keep
        fi
    done
    return 1  # No matching suffix found
}

# Function to remove files and directories without user confirmation
function remove_files_and_directories {
    local path="$1"
    if [[ -d "$path" ]]; then
        for file in "$path"/*; do
            if [[ -f "$file" ]]; then
                if ! has_suffix_to_keep "$file"; then
                    rm "$file"
                    echo "Deleted file: $file"
                fi
            elif [[ -d "$file" ]]; then
                remove_files_and_directories "$file"
            fi
        done
    fi
}

# Start removing files and directories without user confirmation
remove_files_and_directories "$home_dir_path"
echo "Cleanup complete."
