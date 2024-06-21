# Define the number of each element
elements = {
    "Mo": 9,
    "Nb": 9,
    "W": 9,
    "Ta": 9
}

# Generate the sed commands
current_index = 0
for element, count in elements.items():
    for i in range(count):
        print(f"sed -i 's/atom  {current_index} /{element}{i + 1} /g' out.dat")
        current_index += 1

# Generate the grep commands
for element1 in elements.keys():
    for element2 in elements.keys():
        print(f"grep -P '{element1}[0-9]+\\s*\\(000\\)<--\\>{element2}[0-9]+' out.dat -A7 > hop.{element1}-{element2}.dat")
