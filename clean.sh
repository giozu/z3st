#!/bin/bash
# --.. ..- .-.. .-.. ---
# Z3ST cleanup script
# --.. ..- .-.. .-.. ---

echo "Cleaning caches, build artifacts, and junk..."

# Remove Python bytecode and cache
find . -type d -name "build" -exec rm -rf {} +
find . -type d -name "__pycache__" -exec rm -rf {} +
find . -type f -name "non-regression.json" -exec rm -f {} +
find . -type f -name '*:Zone.Identifier' -exec rm -f {} \;

echo "Cleanup complete."
