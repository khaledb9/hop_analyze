#!/usr/bin/env python3
import argparse
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Rename atom indices and extract hopping data from a file."
    )
    parser.add_argument(
        "file",
        nargs="?",
        default="out.dat",
        help="Input data file (default: out.dat)"
    )
    args = parser.parse_args()
    input_file = args.file

    # Verify the file exists
    try:
        with open(input_file, 'r'):
            pass
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)

    # Define the number of each element
    elements = {
        "Cr": 8,
        "Se": 16
    }

    # Apply the sed replacements in-place
    current_index = 0
    for element, count in elements.items():
        for i in range(count):
            # On macOS/BSD sed, change `-i` to `-i ''` if needed.
            sed_cmd = (
                f"sed -i 's/atom  {current_index} /{element}{i + 1} /g' {input_file}"
            )
            subprocess.run(sed_cmd, shell=True, check=True)
            current_index += 1

    # Extract hopping sections with grep
    for el1 in elements:
        for el2 in elements:
            grep_cmd = (
                f"grep -P '{el1}[0-9]+\\s*\\(000\\)<-->"
                f"{el2}[0-9]+' {input_file} -A7 > hop.{el1}-{el2}.dat"
            )
            subprocess.run(grep_cmd, shell=True, check=True)

    print(f"Processing of '{input_file}' completed successfully.")

if __name__ == "__main__":
    main()
