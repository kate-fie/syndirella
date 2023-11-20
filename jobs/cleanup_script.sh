#!/bin/bash

# USAGE:
# ./cleanup_script.sh <home_directory_path>

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

# Function to remove files and directories with user confirmation
function remove_files_and_directories {
    local path="$1"
    if [[ -d "$path" ]]; then
        for file in "$path"/*; do
            if [[ -f "$file" ]]; then
                if ! has_suffix_to_keep "$file"; then
                    read -p "Delete file '$file'? (y/n): " choice
                    if [[ "$choice" == "y" ]]; then
                        rm "$file"
                        echo "Deleted file: $file"
                    else
                        echo "Skipped file: $file"
                    fi
                fi
            elif [[ -d "$file" ]]; then
                remove_files_and_directories "$file"
            fi
        done
    fi
}

# Check if the home directory path argument is provided
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <home_directory_path>"
    exit 1
fi

home_dir_path="$1"
suffixes_to_keep=( '.minimised.json' '.holo_minimised.pdb' '.minimised.mol' )

# Start removing files and directories
read -p "This script will remove files and directories. Proceed? (y/n): " confirm
if [[ "$confirm" == "y" ]]; then
    remove_files_and_directories "$home_dir_path"
    echo "Cleanup complete."
else
    echo "Cleanup aborted."
fi
