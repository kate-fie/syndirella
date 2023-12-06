#!/bin/bash

FOLDER_PATH="/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/D68EV3CPROA/elabs/1_step"
TARGET1_PATH="/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/D68EV3CPROA/elabs/1_step_1-1"
TARGET2_PATH="/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/D68EV3CPROA/elabs/1_step_1-2"
TARGET3_PATH="/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/D68EV3CPROA/elabs/1_step_1-3"
TARGET4_PATH="/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/D68EV3CPROA/elabs/1_step_1-4"

# Get a list of directories
dirs=($(find "$FOLDER_PATH" -maxdepth 1 -type d | tail -n +2))

# Calculate quarter points for directory distribution
quarter=$(( (${#dirs[@]} + 3) / 4 ))  # +3 to round up and divide more evenly

# Move directories to the four targets
for ((i=0; i < $quarter; i++)); do
    mv "${dirs[$i]}" "$TARGET1_PATH"
done

for ((i=$quarter; i < 2*$quarter; i++)); do
    mv "${dirs[$i]}" "$TARGET2_PATH"
done

for ((i=2*$quarter; i < 3*$quarter; i++)); do
    mv "${dirs[$i]}" "$TARGET3_PATH"
done

for ((i=3*$quarter; i < ${#dirs[@]}; i++)); do
    mv "${dirs[$i]}" "$TARGET4_PATH"
done

echo "Directories have been moved."
