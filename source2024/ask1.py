import random
import sys


sys.stdout.reconfigure(encoding='utf-8')

def random_symbols(n):
    alphabet = ['A', 'C', 'G', 'T'] # Αλφάβητο
    return ''.join(random.choices(alphabet, k=n)) # Επιλέγουμε n τυχαία σύμβολα από το αλφάβητο

def synthesize_string():
    num_symbols = random.randint(1, 3) # Επιλέγουμε τυχαία 1 έως 3 σύμβολα
    symbols = random_symbols(num_symbols) # Παίρνουμε τα τυχαία σύμβολα
    return symbols

def main_code_to_generate_a_string():
    # Δοκιμή με τα δεδομένα patterns
    patterns = ["AATTGA", "CGCTTAT", "GGACTCAT", "TTATTCGTA"]

    string = synthesize_string() # erotima (i-a)

    version = 1
    for pattern in patterns: # erotima (i-b)
        num_symbols = random.randint(1, 2) # πόσα σύμβολα θα αντικατστήσουμε(το πολύ 2)
        for i in range(1, num_symbols + 1):
            x = int(random.randint(1, len(pattern) - 1 ))
            choices = ['A', 'C', 'G', 'T', '']
            choices.remove(pattern[x])  #έτσι είμαστε σίγουροι ότι δεν θα αντικαταστήσει ένα γράμμα με τον εαυτό του 
            extend_string = pattern[:x-1] + random.choice(choices) + pattern[x:] # αντικατάσταση με ένα άλλο τυχαία επιλεγμένο σύμβολο είτε με κενη συμβολοσειρα(διαγραφή)
            pattern = extend_string
        string += extend_string
        version += 1

    for i in range(1, len(patterns)): # erotima (i-c)
        string = string + random_symbols(random.randint(1, 2))
    return string

strings = []

for i in range(50):
    strings.append(main_code_to_generate_a_string())

# Ανακατεύουμε την λίστα με τα 50 strings
random.shuffle(strings)

# Χωρίζουμε τη λίστα σε δύο μέρη: datasetA και datasetB
datasetA = strings[:15]
datasetB = strings[15:]

# Δημιουργία του αρχείου FullDataset.txt
with open("auxiliary2024/FullDataset.txt", "w", encoding="utf-8") as file:
    file.write("\n".join(strings))
    print(" - - Δημιουργήθηκε το αρχείο FullDataset - - ")

# Δημιουργία του αρχείου datasetA.txt
with open("auxiliary2024/datasetA.txt", "w", encoding="utf-8") as file:
    file.write("\n".join(datasetA))
    print(" - - Δημιουργήθηκε το αρχείο datasetA - - ")

# Δημιουργία του αρχείου datasetB.txt
with open("auxiliary2024/datasetB.txt", "w", encoding="utf-8") as file:
    file.write("\n".join(datasetB))
    print(" - - Δημιουργήθηκε το αρχείο datasetB - - ")