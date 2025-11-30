#!/bin/bash
# Simple script to activate venv and run main.py

# Activate virtual environment
source venv/bin/activate

# Run main.py with any arguments passed to the script
python src/main.py "$@"
