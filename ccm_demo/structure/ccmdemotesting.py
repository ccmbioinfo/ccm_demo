import string
def protein_prediction_complex(sequences, stoichiometry): 
    if type(stoichiometry) == list and len(stoichiometry) != len(sequences):
                raise ValueError("Enter one stoichiometry value for every sequence.")
        
    # 1. Generate alphabetical IDs
    def generate_alpha_IDs(n):
        alpha_IDs = list(string.ascii_uppercase)
        if n>26:
            for first in string.ascii_uppercase:
                for second in string.ascii_uppercase: 
                    alpha_IDs.append(first + second)
                    if len(alpha_IDs) >= n:
                        return alpha_IDs
        return alpha_IDs[:n]

    # Calculate number of chains per complex, generate alpha IDs accordingly
    number_of_chains = 0
    index = 0 
    for value in stoichiometry: 
        number_of_chains = number_of_chains + int(stoichiometry[index])
    alpha_ID_list = generate_alpha_IDs(number_of_chains)

    # assign alpha IDs to each sequence in the list
    
    assigned_alpha_IDs = []
    for seq in sequences:
        assigned_IDs = alpha_ID_list[index : index+stoichiometry[index]]
        index += stoichiometry[index]
        assigned_alpha_IDs.append({
            'Sequence': seq,
            'Stoichiometry': stoichiometry[index],
            'Assigned Alpha IDs': ",".join(assigned_IDs)
        })
    print(assigned_alpha_IDs)    


protein_prediction_complex(sequences = ["ACC", "CGGA", "GTAG"], stoichiometry=[2, 3, 1])