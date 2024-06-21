import re
import numpy as np
import sys

def parse_hopping_data(filename):
    with open(filename, 'r') as file:
        content = file.read()

    # Split the content into blocks
    blocks = content.split('--\n')

    # Compile regex for extracting values and atoms
    matrix_regex = re.compile(r'\((-?\d+\.\d+)\s*([+-]\s*\d+\.\d+)i\)')
    atoms_regex = re.compile(r'between\s+(\w+\d+)\s*\(.*?\)\s*<-->\s*(\w+\d+)\s*\(.*?\)')
    radius_regex = re.compile(r'radius\s+(\d+\.\d+)')
    results = []

    for block in blocks:
        if not block.strip():
            continue

        # Find the atoms
        atoms_match = atoms_regex.search(block)
        if atoms_match:
            atoms = f"{atoms_match.group(1)}-{atoms_match.group(2)}"
            atom_pair = f"{atoms_match.group(1)}-{atoms_match.group(2)}"
        else:
            atoms = "Unknown"  # Default to Unknown if no match is found
            atom_pair = "Unknown"
            # Debugging: Print a message if atom pairs are not found
            print(f"Atom pairs not found in block: {block}")

        # Find the radius and decide the group
        radius_match = radius_regex.search(block)
        if radius_match:
            radius = float(radius_match.group(1))
        else:
            # If radius is not found, skip this block
            print(f"Radius not found in block: {block}")
            continue

        if 1 <= radius <= 3.5:
            group = '1NN'
        elif 4 <= radius <= 5.9:
            group = '2NN'
        elif 6 <= radius <= 6.7:
            group = '3NN'
        elif radius <= 3:
            group = 'On-Site'
        else:
            continue  # skip the blocks that do not fit any group criteria

        # Find all matrix elements
        elements = matrix_regex.findall(block)
        max_abs_value = 0
        max_phase = 0

        # Calculate the max absolute value and corresponding phase
        for elem in elements:
            real_part = float(elem[0])
            imag_part = float(re.sub(r'\s+', '', elem[1]))  # Remove spaces in the imaginary part
            abs_value = np.abs(complex(real_part, imag_part))
            phase = np.angle(complex(real_part, imag_part)) / np.pi
            if abs_value > max_abs_value:
                max_abs_value = abs_value
                max_phase = phase

        results.append((atoms, max_abs_value, max_phase, group, radius, atom_pair))

    # Sort results by the group
    results.sort(key=lambda x: ('1NN', '2NN', '3NN', 'On-Site').index(x[3]))
    return results

def main():
    if len(sys.argv) < 2:
        print("Usage: python script_name.py filename")
        sys.exit(1)
    filename = sys.argv[1]
    result_list = parse_hopping_data(filename)
    print("Atoms Pair     Max Abs Value (Hopping)    Phase              Group       Radius")
    for atoms, value, phase, group, radius, atom_pair in result_list:
        print(f"{atom_pair:<15} {value:.4f}                    {phase:.4f}π             {group}    {radius}")

if __name__ == "__main__":
    main()
