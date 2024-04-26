#!/bin/bash

# Directories
SRC_DIR="F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/results"
DEST_DIR="F:/My git/paper_code/2023_Imaging_Transcriptomics/results"

# Navigate to the destination directory
cd "$DEST_DIR" || exit

# Loop through files larger than 10MB in the source directory
find "$SRC_DIR" -type f -size +10M | while read -r file; do
    # Copy the large file to the destination directory
    cp "$file" "$DEST_DIR"

    # Extract the filename for commit message
    FILENAME=$(basename "$file")

    # Add, commit, and push the file to Git
    git add .
    git commit -m "Adding file $FILENAME"
    git push

    # Optional sleep to avoid bombarding the server
    # sleep 10s
done
