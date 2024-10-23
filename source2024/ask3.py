from collections import Counter
import numpy as np

file_name = 'auxiliary2024/multiple_alignment_result.txt'

try:
    with open(file_name, 'r') as file:
        # Read all lines in the file
        lines = file.readlines()
        
        # Strip whitespace characters and newline characters from each line, then store in datasetA
        multiple_alignment = [line.strip() for line in lines]
        
except FileNotFoundError:
    print(f"The file {file_name} does not exist.")
except Exception as e:
    print(f"An error occurred: {e}")


file_name = 'auxiliary2024/datasetB.txt'

try:
    with open(file_name, 'r') as file:
        # Read all lines in the file
        lines = file.readlines()
        
        # Strip whitespace characters and newline characters from each line, then store in datasetA
        sequences = [line.strip() for line in lines]
        
except FileNotFoundError:
    print(f"The file {file_name} does not exist.")
except Exception as e:
    print(f"An error occurred: {e}")



threshold = 70
num_rows = len(multiple_alignment)
num_cols = len(multiple_alignment[0])
alphabeta = ['A', 'C', 'T', 'G']
conserved_states = []
deletion_states = []
insertion_states = []
hmmstates = []
states = []


def check_region(align):
    clean_align = align.replace("-", "")
    symbol_count = len(clean_align)
    conserved_regions = False
    deletions = False
    insertions = False

    # Calculate the percentage of letters (non-dash characters)
    if (symbol_count / len(align)) * 100 >= threshold:
        conserved_regions = True
        if symbol_count < len(align):
            deletions = True
    else:
        insertions = True

    return conserved_regions, deletions, insertions

def search_list(lst, element):
    i = 0
    while i < len(lst):
        if lst[i] == element:
            return True
        i += 1
    return False

def take_the_column(n):
    col = ""
    for i in range(len(multiple_alignment)):
        col += multiple_alignment[i][n]
    return col

def hmm_profile(multiple_alignment):
    for i in range(num_cols):
        conserved, deletion, insertion = check_region(take_the_column(i))
        states.append(i)

        if (conserved): conserved_states.append(i)
        if (deletion): deletion_states.append(i)
        if (insertion): insertion_states.append(i)

    for i in range(num_cols):
        if not (search_list(deletion_states, states[i])) and search_list(conserved_states, states[i]) :
            hmmstates.append("match") 
        elif search_list(deletion_states, states[i]) and search_list(conserved_states, states[i]) :
            hmmstates.append("delete")
        elif search_list(insertion_states, states[i] ):
            hmmstates.append("insert")
    
hmm_profile(multiple_alignment)


#Emmision Propability Table for HMM Profile
def create_Emmision_Prob_table():

    emission_table = {symbol: [0 for _ in range(len(conserved_states))] for symbol in alphabeta}
    
    conserved_index = 0  
    for i in range(num_cols):
        if search_list(conserved_states, states[i]):
            col = take_the_column(i)
            col_no_gaps = col.replace("-", "")  
            symbol_counts = Counter(col_no_gaps)
            total_symbols = len(col_no_gaps)
            
            for symbol in alphabeta:
                if symbol in symbol_counts:
                    prob = symbol_counts[symbol] / total_symbols
                else:
                    prob = 0
                emission_table[symbol][conserved_index] = round(prob, 3)
            conserved_index += 1  

    print("\nEmission Probability Table:")
    for symbol in alphabeta:
        print(f"{symbol}: {emission_table[symbol]}")


create_Emmision_Prob_table()
print("")


transition_probs = {
    'M': {},
    'I': {},
    'D': {},
}


def create_Transition_Prob_table():
    j = 0
    prob_table = [[0 for _ in range(len(conserved_states)-1)] for _ in range(9)]
    print("HMM Profile: ", hmmstates, "\n")
    count = [0 for _ in range(num_rows)] 

    for i in range(len(hmmstates) - 1): 

        if hmmstates[i] == "match" and hmmstates[i + 1] == "delete":
            count_letters = take_the_column(i+1)
            count_letters = count_letters.replace("-", "")
            count_letters = len(count_letters)/num_rows
            prob_table[0][j] = round(count_letters, 2)
            prob_table[2][j] = round(1 - count_letters, 2)
            j += 1
        elif hmmstates[i] == "match" and hmmstates[i + 1] == "insert":
            counter_of_inserts = 0
            prob_counter = 0
            i += 1   
            
            while(hmmstates[i] == "insert"): 
                inserts = take_the_column(i)

                for k in range(len(inserts)):
                    if inserts[k] == "-": count[k] += 1 
                counter_of_inserts += 1
                i += 1
            
            for k in range(len(inserts)):
                if count[k] == 0: prob_counter += 1  
            

            if counter_of_inserts > 2:  
                prob_table[1][j] = round(prob_counter / num_rows, 2)
                prob_table[3][j] = round(prob_counter / counter_of_inserts, 2)  
                prob_table[0][j] = round(1 - prob_table[1][j], 2)
                if(i+1 < len(hmmstates)) and  hmmstates[i + 1] != "delete": 
                    prob_table[3][j] = round(prob_counter / counter_of_inserts, 2)  
                    prob_table[4][j] = round(1 - prob_table[3][j], 2) 
                if(i+1 < len(hmmstates)) and  hmmstates[i + 1] == "delete": 
                    prob_table[4][j] = round((1 - prob_table[3][j])/2, 2) 
                    prob_table[7][j] = round((1 - prob_table[3][j])/2, 2)
            else:
                prob_table[0][j] = round((1 - prob_table[3][j])/3, 2)
                prob_table[4][j] = round((1 - prob_table[3][j])/3, 2) 
                prob_table[7][j] = round((1 - prob_table[3][j])/3, 2)
            j += 1
        elif hmmstates[i] == "match" and hmmstates[i + 1] == "match":
            prob_table[0][j] = 1
            j += 1
        elif hmmstates[i] == "delete" and hmmstates[i + 1] == "delete":
            count_letters = take_the_column(i+1)
            count_letters = count_letters.replace("-", "")
            count_letters = len(count_letters)/num_rows
            prob_table[0][j] = round(count_letters, 2)
            prob_table[5][j] = round(1 - count_letters, 2)
            if i + 2 < len(hmmstates) and hmmstates[i + 2] == "insert": prob_table[8][j] = round(prob_table[5][j], 2)
            j += 1
        elif hmmstates[i] == "delete" and hmmstates[i + 1] == "match":
            prob_table[6][j] = prob_table[2][j-1]
            prob_table[0][j] = prob_table[0][j-1] 
            j += 1
        elif hmmstates[i] == "delete" and hmmstates[i + 1] == "insert":
            counter_of_inserts = 0
            prob_counter = 0

            while i + 1 < len(hmmstates) and hmmstates[i+1] == "insert":
                inserts = take_the_column(i)

                for k in range(len(inserts)):
                    if inserts[k] == "-": 
                        count[k] += 1
                counter_of_inserts += 1
                i += 1

            for k in range(len(inserts)):
                if count[k] == 0: 
                    prob_counter += 1

            if j < len(prob_table[0]):
                if counter_of_inserts > 0:
                    if counter_of_inserts < prob_counter :prob_table[3][j] = round(counter_of_inserts / prob_counter, 2) #i -> i
                    if counter_of_inserts >= prob_counter :prob_table[3][j] = round(prob_counter / counter_of_inserts, 2) #i -> i
                    
                    if j+1 < len(prob_table[0]):
                        prob_table[0][j+1] = round(prob_table[4][j], 2) #m -> m
                        prob_table[4][j] = round((1 - prob_table[3][j])/3, 2)  #i -> m
                        prob_table[7][j] = round(prob_table[4][j], 2)  #i -> d
                else:
                    prob_table[7][j] = 0.5
                    prob_table[4][j] = 0.5
                    
            j += 1

    print("------Transition Probabilities------")
    for i in range(len(prob_table[0])):
        if prob_table[0][i] != 0: print("Transition: M" + str(i + 1) + " to M" + str(i + 2) + ": " + str(prob_table[0][i]))
        if prob_table[2][i] != 0: print("Transition: M" + str(i + 1) + " to D" + str(i + 1) + ": " + str(prob_table[2][i]))
        if prob_table[1][i] != 0: print("Transition: M" + str(i + 1) + " to I" + str(i + 1) + ": " + str(prob_table[1][i]))
        if prob_table[3][i] != 0: print("Transition: I" + str(i + 1) + " to I" + str(i + 1) + ": " + str(prob_table[3][i]))
        if prob_table[4][i] != 0: print("Transition: I" + str(i + 1) + " to M" + str(i + 2) + ": " + str(prob_table[4][i]))
        if prob_table[5][i] != 0: print("Transition: D" + str(i) + " to D" + str(i + 1) + ": " + str(prob_table[5][i]))
        if prob_table[6][i] != 0: print("Transition: D" + str(i) + " to M" + str(i + 2) + ": " + str(prob_table[6][i]))
        if prob_table[7][i] != 0: print("Transition: I" + str(i + 1) + " to D" + str(i + 2) + ": " + str(prob_table[7][i]))
        if prob_table[8][i] != 0: print("Transition: D" + str(i+1) + " to I" + str(i + 2) + ": " + str(prob_table[8][i]))

        print("--------------------------------")
    results = []
    for col in range(len(prob_table)):
        temp = 0
        for j in range(len(prob_table[col])):
            temp += prob_table[col][j]
        results.append(round(temp/len(prob_table[col]), 3))
    
    global transition_probs  
    transition_probs = {
        'M': {'M': results[0], 'I': results[1], 'D': results[2]},
        'I': {'M': results[3], 'I': results[4], 'D': results[7]},  
        'D': {'M': results[5], 'I': results[6], 'D': results[8]},  
    }


create_Transition_Prob_table()



def Emmisions_Prob_for_Viterbi(): 
    emission_probs = {'M': [], 'I': [], 'D': []}

    total_A, total_C, total_T, total_G = 0, 0, 0, 0

    for col in conserved_states:
        for row in multiple_alignment:
            if row[col] == 'A': total_A += 1
            if row[col] == 'C': total_C += 1
            if row[col] == 'T': total_T += 1
            if row[col] == 'G': total_G += 1

    total = (len(conserved_states))*num_rows - 1

    overall_prob_dict = {
        'A': round(total_A / total, 2),
        'C': round(total_C / total, 2),
        'G': round(total_G / total, 2),
        'T': round(total_T / total, 2)
    }

    emission_probs['M'].append(overall_prob_dict)

    total_A, total_C, total_T, total_G = 0, 0, 0, 0

    for col in insertion_states:
        for row in multiple_alignment:
            if row[col] == 'A': total_A += 1
            if row[col] == 'C': total_C += 1
            if row[col] == 'T': total_T += 1
            if row[col] == 'G': total_G += 1

    total = total_A + total_C + total_T + total_G

    if total > 0:
        overall_prob_dict = {
            'A':  round(total_A / total, 2),
            'C':  round(total_C / total, 2),
            'G':  round(total_G / total, 2),
            'T':  round(total_T / total, 2)
        }
    else:
        overall_prob_dict = {
            'A': 0.0,
            'C': 0.0,
            'G': 0.0,
            'T': 0.0
        }

    emission_probs['I'].append(overall_prob_dict)

    total_A, total_C, total_T, total_G, total = 0, 0, 0, 0, 0
    for col in deletion_states:
        for row_idx, row in enumerate(multiple_alignment):
            if row[col] == '-':
                total += 1
                back_col = col
                while back_col >= 0 and hmmstates[back_col] != 'match':
                    back_col -= 1
                if back_col >= 0:
                    match_nucleotide = multiple_alignment[row_idx][back_col]
                    if match_nucleotide == 'A': total_A += 1
                    if match_nucleotide == 'C': total_C += 1
                    if match_nucleotide == 'T': total_T += 1
                    if match_nucleotide == 'G': total_G += 1

    if total > 0:
        overall_prob_dict = {
            'A':  round(total_A / total, 2),
            'C':  round(total_C / total, 2),
            'G':  round(total_G / total, 2),
            'T':  round(total_T / total, 2)
        }
    else:
        overall_prob_dict = {
            'A': 0.0,
            'C': 0.0,
            'G': 0.0,
            'T': 0.0
        }
    emission_probs['D'].append(overall_prob_dict)

    return emission_probs

emission_probs = Emmisions_Prob_for_Viterbi() 
print("Emission Probabilities for Viterbi: \n", emission_probs, "\n")
print("Emission Probabilities for Viterbi: \n", transition_probs, "\n")
