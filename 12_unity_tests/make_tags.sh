#!/bin/bash
# Script to generate TAGS file for Emacs

echo "Generating TAGS file for Emacs..."

# Find all C and H files and generate TAGS
find . -type f \( -name "*.c" -o -name "*.h" \) \
    ! -path "*/build*/*" \
    ! -path "*/.git/*" \
    ! -path "*/out/*" \
    -print | etags -

echo "TAGS file created successfully!"
echo "Use M-. (find-tag) in Emacs to navigate to definitions"
