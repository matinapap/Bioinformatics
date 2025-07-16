import random
import sys

# Ensure the output encoding supports UTF-8
sys.stdout.reconfigure(encoding='utf-8')

def random_symbols(n):
    alphabet = ['A', 'C', 'G', 'T']  # Alphabet
    return ''.join(random.choices(alphabet, k=n)) # Select n random symbols from the alphabet

def synthesize_string():
    num_symbols = random.randint(1, 3) # Choose randomly between 1 to 3 symbols
    symbols = random_symbols(num_symbols) # Get the random symbols
    return symbols

def main_code_to_generate_a_string():
    # Test with the given patterns
    patterns = ["AATTGA", "CGCTTAT", "GGACTCAT", "TTATTCGTA"]

    string = synthesize_string() 

    version = 1
    for pattern in patterns: 
        num_symbols = random.randint(1, 2) # How many symbols to replace (up to 2)
        for i in range(1, num_symbols + 1):
            x = int(random.randint(1, len(pattern) - 1 ))
            choices = ['A', 'C', 'G', 'T', '']
            choices.remove(pattern[x])  # Ensure we don't replace a letter with itself 
            extend_string = pattern[:x-1] + random.choice(choices) + pattern[x:] # Replace with a randomly selected different symbol or with an empty string (deletion)
            pattern = extend_string
        string += extend_string
        version += 1

    for i in range(1, len(patterns)): 
        string = string + random_symbols(random.randint(1, 2))
    return string

strings = []

for i in range(50):
    strings.append(main_code_to_generate_a_string())

# Generate 50 synthetic strings
random.shuffle(strings)

# Split the list into two parts: datasetA and datasetB
datasetA = strings[:15]
datasetB = strings[15:]

# Create the FullDataset.txt file
with open("auxiliary2024/FullDataset.txt", "w", encoding="utf-8") as file:
    file.write("\n".join(strings))
    print(" - - Δημιουργήθηκε το αρχείο FullDataset - - ")

# Create the datasetA.txt file
with open("auxiliary2024/datasetA.txt", "w", encoding="utf-8") as file:
    file.write("\n".join(datasetA))
    print(" - - Δημιουργήθηκε το αρχείο datasetA - - ")

# Create the datasetB.txt file
with open("auxiliary2024/datasetB.txt", "w", encoding="utf-8") as file:
    file.write("\n".join(datasetB))
    print(" - - Δημιουργήθηκε το αρχείο datasetB - - ")
