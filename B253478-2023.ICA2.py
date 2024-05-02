#!/usr/bin/python3

# Import necessary modules for the script
import os            # To interact with the operating system
import re            # For regular expression operations
import sys           # For system-specific parameters and functions
import subprocess    # To run external commands and programs

# ↓ Start: Public function area

def ask_output_folder(folder_name):
    '''
    This function prompts the user to input a path for saving output data.

    Parameters:
    folder_name:the name of the folder for which the path is being asked.
    '''

    while True:
        output_folder = input(f"Please enter the output path you want to save the {folder_name} to: ").rstrip(" /\\")
        output_folder = output_folder.replace(" ", "_").replace("/", "_").replace("\\", "_")
        try:
            os.makedirs(output_folder, exist_ok=True)
        except Exception as e:
            print(f"An error occurred while creating the folder: {e}. Please try again.")
        else:
            print(f'The output path for {folder_name} has been set: {output_folder}')
            print('#---------------------------')
            return output_folder

# ↑ End: Public function area

# ↓ Start: Obtain the sequences dataset

def get_data():
    """
    This function guides users to enter search criteria and obtain corresponding data from NCBI, and save the data to a file.
    It also requires users to confirm whether they wish to continue processing the extracted data.

    Returns:
    the path to the fetched FASTA file, taxonomic group, and protein family.
    """

    print('#===========================\n')
    while True:
        # Obtain user input for search criteria
        taxonomic_group, protein_family, partial, maxnum = get_user_input()

        # Fetch data from NCBI based on the input criteria
        fasta_data = run_edirect(taxonomic_group, protein_family, partial, maxnum)

        # Analyze FASTA data to get sequence and species counts
        sequence_count, species_count, species_set, processed_sequences = parse_fasta(fasta_data, maxnum)

        # Inform the user about the number of sequences and species in the dataset
        print(f"The dataset contains {sequence_count} sequences from {species_count} species.")
        if sequence_count == 0:
            print('Please choose a different dataset.')
            continue
        # Ask the user for the directory to save the fetched protein data
        fasta_output = ask_output_folder("protein data")

        # Save the fetched FASTA data to a file for the user to review
        fasta_file = save_fasta_to_file(processed_sequences, protein_family, taxonomic_group, fasta_output)

        # Ask the user if they are satisfied with the data, loop back if not
        if ask_continue(fasta_file):
            return fasta_file, taxonomic_group, protein_family


def get_user_input():
    """
    Prompts the user to input search criteria for fetching protein data, including the taxonomic group,
    protein family, whether to include partial sequences, and the maximum number of sequences to fetch.

    Returns:
    the taxonomic group, protein family, partial sequence inclusion flag, and maximum number.
    """

    # Request user input for the taxonomic group of interest
    taxonomic_group = input(
        "Please enter the taxonomic group you are interested in (e.g., Aves, Mammalia, Rodentia, Vertebrata): ")

    # Request user input for the protein family of interest
    protein_family = input(
        "Please enter the protein family you want to analyze (e.g., glucose-6-phosphatase, kinases, cyclases, transporters): ")

    # Loop to ensure correct input for including partial sequences
    while True:
        partial = input(
            "Do you want to keep the partial sequence? (Y/N): ").strip().upper()
        if partial in ["Y", "N"]:
            break
        else:
            print("Unable to recognize your input, please follow the prompts to enter!")

    # Loop to ensure a valid positive integer for the maximum number of sequences
    while True:
        try:
            maxnum = int(input(
                "Please enter the maximum number of sequences you want to search for (Recommended: 1000): "))
            if maxnum > 0:
                break
            else:
                print("Please enter a positive integer.")
        except ValueError:
            print("Invalid input. Please enter a valid positive integer.")

    # Inform the user about the search criteria based on their input
    if partial == "Y":
        print(f"We will search in NCBI: {protein_family}[Title] AND {taxonomic_group}[Organism]，max amount = {maxnum}...")
    elif partial == "N":
        print(f"We will search in NCBI: {protein_family}[Title] AND {taxonomic_group}[Organism] NOT partial，max amount = {maxnum}...")

    print('#---------------------------')
    return taxonomic_group, protein_family, partial, maxnum

def run_edirect(taxonomic_group, protein_family, partial, maxnum):
    """
    Runs the NCBI EDirect command-line tool to fetch protein data based on specified criteria.

    Parameters:
    taxonomic_group (str): The taxonomic group to filter the data.
    protein_family (str): The protein family of interest.
    partial (str): Flag to include ('Y') or exclude ('N') partial sequences.
    maxnum (int): Maximum number of sequences to fetch.

    Returns:
    The fetched FASTA formatted data as a string.
    """

    # Construct the NCBI query command based on user preferences for partial sequences
    if partial == "Y":
        cmd = f'esearch -db protein -query "{protein_family}[Title] AND {taxonomic_group}[Organism]" -retmax {maxnum} | efetch -format fasta'
    elif partial == "N":
        cmd = f'esearch -db protein -query "{protein_family}[Title] AND {taxonomic_group}[Organism] NOT partial" -retmax {maxnum} | efetch -format fasta'

    # Execute the constructed command and handle potential errors
    try:
        result = subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        return None
    else:
        return result.stdout


def parse_fasta(fasta_data, maxnum):
    """
    Parses the fetched FASTA data to count the number of sequences and species present.

    Parameters:
    fasta_data (str): The FASTA formatted data as a string.
    maxnum (int): The maximum number of sequences to consider.

    Returns:
    the sequence count, species count, and a set of species.
    """

    # Regular expression to match FASTA headers and extract species information.
    species_pattern = re.compile(r'\[([^\]]+)\]')
    species_set = set()
    sequence_count = 0
    processed_sequences = ""  # Initialize an empty string to store sequences

    # Splitting the FASTA data into individual sequences. Each sequence in a FASTA file is preceded by '>'.
    # The strip and split operations isolate each sequence for individual analysis.
    sequences = fasta_data.strip().split('>')

    # Iterating through each sequence in the split FASTA data.
    for seq in sequences:

        if seq:  # Check if the sequence is not an empty string
            sequence_count += 1  # Increment the sequence count for each valid sequence found.
            header, sequence = seq.split('\n', 1)  # Splitting header and sequence
            processed_sequences += ">" + header + "\n" + sequence + "\n"  # Append sequence

            # Use the regular expression to search for species name in the sequence header.
            species_match = species_pattern.search(seq)
            if species_match:
                species = species_match.group(1)  # Extract the species name.
                species_set.add(species)  # Add the species name to the set, ensuring uniqueness.

            # Stop processing if the number of sequences processed reaches the user-defined maximum.
            if sequence_count == maxnum:
                break

    return sequence_count, len(species_set), species_set, processed_sequences


def save_fasta_to_file(fasta_data, protein_family, taxonomic_group, directory):
    """
    Saves the fetched FASTA data to a file in the specified directory.
    Formats the filename based on the protein family and taxonomic group.

    Parameters:
    fasta_data (str): The FASTA formatted data to be saved.
    protein_family (str): The protein family of interest.
    taxonomic_group (str): The taxonomic group of interest.
    directory (str): The directory where the FASTA file will be saved.

    Returns:
    The path to the saved FASTA file.
    """

    # Format the filename using the protein family and taxonomic group
    file_name = f"{protein_family}_in_{taxonomic_group}.fasta"
    # Replace characters that might cause file system issues
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    # Create the full path for the file
    fasta_path = os.path.join(directory, file_name)

    # Writing the FASTA data to the file
    with open(fasta_path, 'w') as file:
        file.write(fasta_data)

    print(f"Fasta data has been saved to: {fasta_path}. You can check it now.")

    return fasta_path

def ask_continue(fasta_file):
    """
    Asks the user whether they wish to continue with the current dataset.
    Provides an option to either continue with the fetched data or choose a different dataset.

    Parameters:
    fasta_file (str): The path to the FASTA file containing the current dataset.

    Returns:
    bool: True if the user chooses to continue with the current dataset, False otherwise.
    """

    while True:
        # Prompt the user to confirm if they want to proceed with the current dataset
        choice = input(f"Do you want to continue with this dataset: {fasta_file} ? (Y/N): ").strip().upper()

        # If the user chooses to continue, confirm and exit the loop
        if choice == 'Y':
            print("Continuing with the dataset...")
            print('#===========================')
            return True

        # If the user chooses not to continue, inform them to select a different dataset and exit the loop
        elif choice == 'N':
            print('#===========================')
            print("Please choose a different dataset.")
            return False

        # Handle invalid input and prompt the user again
        else:
            print("Unable to recognize your input, please follow the prompts to enter!")

# ↑ End: Obtain the sequences dataset

# ↓ Start: Conservation analysis

def conservation_analysis(protein_family, taxonomic_group, fasta_file):
    """
    Performs conservation analysis on protein sequences. This includes aligning sequences,
    running infoalign for analysis, filtering sequences based on conservation, and visualizing
    conservation levels.

    Parameters:
    protein_family (str): Name of the protein family.
    taxonomic_group (str): Taxonomic group of interest.
    fasta_file (str): File path to the FASTA file.

    Returns:
    Path to the aligned sequence file.
    """

    print('#===========================')
    print("\nWe will determine and plot the level of conservation between the protein sequences！")
    con_output = ask_output_folder('Conservation Analysis')

    # Align the sequences using Clustal Omega
    aligned_file = align_sequence(protein_family, taxonomic_group, fasta_file, 'aligned', con_output)

    # Run infoalign and capture the output
    infoalign_output = run_infoalign(protein_family, taxonomic_group, aligned_file, con_output)

    if infoalign_output:
        # Parse the output from infoalign to extract sequence information
        sequences_info = parse_infoalign_output(infoalign_output)
        while True:
            # Ask the user for criteria to filter sequences based on conservation
            filtered_sequence_ids = ask_filter(sequences_info)
            if filtered_sequence_ids != []:
                # Extract the filtered sequences from the original FASTA file
                filtered_sequences = extract_sequences(fasta_file, filtered_sequence_ids)
                break
            else:
                print('The filtering criteria are too strict, please re-enter!')

        # Save the filtered sequences to a new FASTA file
        selected_fasta_file = save_selected_fasta(filtered_sequences, protein_family, taxonomic_group, con_output)
        # Align the filtered sequences for further analysis
        aligned_filtered_file = align_sequence(protein_family, taxonomic_group, selected_fasta_file, 'aligned_filtered', con_output)
        # Visualize the conservation levels
        plot_con(protein_family, taxonomic_group, aligned_filtered_file, con_output)
    else:
        print("Error running infoalign.")

    print('#===========================')
    return aligned_file


def align_sequence(protein_family, taxonomic_group, fasta_file, name, directory):
    """
    Aligns sequences using Clustal Omega. Essential for conservation analysis.

    Parameters:
    protein_family (str): Protein family name.
    taxonomic_group (str): Taxonomic group.
    fasta_file (str): FASTA file path.
    name (str): Descriptor for output file.
    directory (str): Output directory.

    Returns:
    Path to the aligned sequence file.
    """

    file_name = f'{protein_family}_in_{taxonomic_group}_{name}.fasta'
    # Replace characters that might cause file system issues
#    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    aligned_file = os.path.join(directory, file_name)

    try:
        # Execute the Clustal Omega command for sequence alignment
        print('Clustal Omega is working for you to align the sequences. Please be patient and wait...')
        subprocess.run(["clustalo", "-i", fasta_file, "-o", aligned_file, "--force"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        print('#---------------------------')
    else:
        print(f"Aligned sequences have been saved to {aligned_file}")
        print('#---------------------------')
        return aligned_file


def run_infoalign(protein_family, taxonomic_group, alignment_file, directory):
    """
    Runs infoalign on aligned sequences for detailed analysis.

    Parameters:
    protein_family (str): Protein family name.
    taxonomic_group (str): Taxonomic group.
    alignment_file (str): File path to aligned sequences.
    directory (str): Output directory.

    Returns:
    Path to the infoalign output file.
    """

    file_name = f'{protein_family}_in_{taxonomic_group}_aligned_info'
    # Replace characters that might cause file system issues
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    infoalign_file = os.path.join(directory, file_name)
    try:
        # Run the infoalign command
        print('Infroalign is working for you to analyse the sequences...')
        subprocess.run(['infoalign', alignment_file, '-out', infoalign_file], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        print('#---------------------------')
    else:
        print(f'Infoalign analyses result has been saved to {infoalign_file}. You can check it now.')
        print('#---------------------------')
        return infoalign_file


def parse_infoalign_output(output_file):
    """
    Extracts sequence identity information from infoalign output.

    Parameters:
    output_file (str): Infoalign output file path.

    Returns:
    Dictionary of sequence IDs and their percent identity.
    """

    print("Parsing infoalign output...")
    sequences_info = {}

    with open(output_file, 'r') as file:
        lines = file.readlines()

    for line in lines[1:]:  # Skipping the header line in the output file
        parts = line.split()  # Splitting each line by whitespace to separate columns
        seq_id = parts[1]  # Extracting the sequence ID, which is in the second column

        parts2 = line.split('\t')  # Splitting the line by tab to access the percent identity
        percent_identity = float(parts2[-3])  # Parsing the percent identity value, which is in the third last column

        sequences_info[seq_id] = percent_identity  # Mapping each sequence ID to its percent identity

    print(f"Alignment sequences have been parsed.")
    print('#---------------------------')
    return sequences_info


def ask_filter(sequences_info):
    """
    Filters sequences based on user-defined percent identity criteria.

    Parameters:
    sequences_info (dict): Sequence IDs and their percent identity.

    Returns:
    List of filtered sequence IDs.
    """

    while True:
        filter = input('Please enter the criteria you want to use to filter for changing values (format example: >70): ')
        if re.match(r'^>[0-9]+', filter):
            # Parse threshold for filtering sequences with identity greater than a specified value
            identity_threshold = float(filter.strip().split()[0][1:])
            filtered_sequence_ids = [seq_id for seq_id, percent_identity in sequences_info.items() if
                                     percent_identity >= identity_threshold]
            print('#---------------------------')
            return filtered_sequence_ids
        elif re.match(r'^<[0-9]+', filter):
            # Parse threshold for filtering sequences with identity less than a specified value
            identity_threshold = float(filter.strip().split()[0][1:])
            filtered_sequence_ids = [seq_id for seq_id, percent_identity in sequences_info.items() if
                                     percent_identity <= identity_threshold]
            print('#---------------------------')
            return filtered_sequence_ids
        else:
            # Handle invalid input format
            print("Invalid filter format. It should start with '>' or '<' followed by one or more digits!")


def extract_sequences(fasta_file, filtered_sequence_ids):
    """
    Extracts sequences from FASTA file based on filtered IDs.

    Parameters:
    fasta_file (str): Original FASTA file path.
    filtered_sequence_ids (list): IDs of sequences to extract.

    Returns:
    Dictionary of filtered sequence data.
    """

    filtered_sequences = {}
    with open(fasta_file, 'r') as f:
        record = None
        for line in f:
            if line.startswith('>'):
                # Identify the header line of each sequence
                record = line.strip().split()[0][1:]  # Extract the sequence ID from the header line

                if record in filtered_sequence_ids:
                    # Initialize a new entry in the dictionary for the filtered sequence
                    filtered_sequences[record] = ''
            elif record and record in filtered_sequence_ids:
                # Append sequence lines to the corresponding sequence ID in the dictionary
                filtered_sequences[record] += line.strip()
    return filtered_sequences


def save_selected_fasta(filtered_sequences, protein_family, taxonomic_group, directory):
    """
    Saves filtered sequences to a new FASTA file.

    Parameters:
    filtered_sequences (dict): Filtered sequence data.
    protein_family (str): Protein family name.
    taxonomic_group (str): Taxonomic group.
    directory (str): Output directory.

    Returns:
    Path to the new FASTA file.
    """

    file_name = f'{protein_family}_in_{taxonomic_group}_filtered.fasta'
    # Replace characters that might cause file system issues
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    filtered_file = os.path.join(directory, file_name)

    with open(filtered_file, 'w') as f:
        for seq_id, sequence in filtered_sequences.items():
            # Write the header line of each sequence
            f.write(f'>{seq_id}\n{sequence}\n')
    print(f"Filtered sequences have been saved to: {filtered_file}. You can check it now.")
    print('#---------------------------')
    return filtered_file


def plot_con(protein_family, taxonomic_group, selected_aligned_file, directory):
    """
    Visualizes conservation levels using plotcon.

    Parameters:
    protein_family (str): Protein family name.
    taxonomic_group (str): Taxonomic group.
    selected_aligned_file (str): File path to aligned sequences.
    directory (str): Output directory.

    Returns:
    None, but saves the plot to a file.
    """
    # Compose the conplot png file name and its path
    file_name = f'{protein_family}_in_{taxonomic_group}_conplot.png'
    # Replace characters that might cause file system issues
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    plotcon_output = os.path.join(directory, file_name)

    try:
        # Run the plotcon tool to generate a conservation plot
        print('Plotcon is working for you to plot the conservation level.')
        subprocess.run(["plotcon", "-sequence", selected_aligned_file, "-graph", "png", "-goutfile", plotcon_output], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
    else:
        print(f'Plot has been saved to {plotcon_output}. You can check it now.')

# ↑ End: Conservation analysis

# ↓ Start: Scan motifs in Prosite

def scan_prosite_motifs(protein_family, taxonomic_group, fasta_file):
    """
    Executes motif scanning on protein sequences using the PROSITE database. This process includes:
    Checking and validating the EMBOSS environment configuration.
    Segmenting the FASTA file into individual sequences.
    Iteratively running patmatmotifs on each sequence to identify motifs.
    Summarizing and extracting motif hits from patmatmotifs output files.

    Parameters:
    protein_family (str): The protein family of interest.
    taxonomic_group (str): The taxonomic group of interest.
    fasta_file (str): Path to the FASTA file with the sequences.

    Returns:
    None, but creates output files with motif scan results.
    """

    print('#===========================')
    print("\nWe will scan protein sequences with motifs from the PROSITE database！")

    # Prompt user for output directory and set it up
    patmatmotifs_output = ask_output_folder("Patmatmotifs")

    # Check if the EMBOSS environment is properly configured
    check_emboss_environment(fasta_file)

    # Prepare for running patmatmotifs by segmenting sequences from the FASTA file
    sequences = cut_fasta(fasta_file)
    # Prepare for display progress bar
    total_sequences = len(sequences)
    processed_sequences = 0

    print(f'There are {total_sequences} sequences we are going to scan.')

    for seq_id, sequence in sequences.items():
        # Run patmatmotifs for each sequence
        run_patmatmotifs(sequence, seq_id, patmatmotifs_output)
        processed_sequences += 1
        # Update and display progress bar after each sequence is processed
        print_progress_bar('Patmatmotifs', processed_sequences, total_sequences)

    print(f"All sequences have been processed, results saved in {patmatmotifs_output}.")
    print('#---------------------------')

    # Extract motif hits from patmatmotifs output
    extract_motif_hits(patmatmotifs_output, protein_family, taxonomic_group)
    print('#===========================')


def check_emboss_environment(fasta_file):
    """
    Checks the EMBOSS environment settings required for motif scanning.

    Parameters:
    fasta_file (str): Path to the FASTA file used for testing the setup.

    Returns:
    None, but prompts user to set EMBOSS_DATA path if not correctly set.
    """

    # Get the current setting of the EMBOSS_DATA environment variable
    emboss_data = os.environ.get('EMBOSS_DATA')
    # Name for a temporary output file used in testing the patmatmotifs command
    test_output = 'test_output'

    while True:

        try:
            # Displaying a message to the user about the check being performed
            print(f'Checking the EMBOSS_DATA path...')

            # Run the patmatmotifs command as a test. This checks if the EMBOSS_DATA environment variable is set correctly.
            result = subprocess.run(["patmatmotifs", "-sequence", fasta_file, "-outfile", test_output],
                                    check=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True)

        except subprocess.CalledProcessError as e:
            # Print the error message indicating a problem with the EMBOSS_DATA path
            print(f"EMBOSS_DATA path isn't correct: {e.stderr}")
            print('#---------------------------')

            # Prompt the user to enter the correct path for EMBOSS_DATA
            if emboss_data is not None:  # If emboss_data is not None, it means the variable was set but incorrectly.
                print(f"Currently EMBOSS_DATA: {emboss_data}")
                emboss_data = input(f'Please enter the correct EMBOSS_DATA path: ')
            else:  # If it is None, it indicates the environment variable was not set at all.
                print("EMBOSS_DATA: Environment variables not declared")
                emboss_data = input(f'Please enter the EMBOSS_DATA path: ')

            # Update the EMBOSS_DATA environment variable with the new path provided by the user
            os.environ['EMBOSS_DATA'] = f'{emboss_data}'

        else:
            # If the patmatmotifs command runs successfully, inform the user that the setup is correct
            print(f'Congratulation! It is correct!')
            print('#---------------------------')

            # Break out of the loop as the EMBOSS environment is correctly set
            break

        finally:
            # Remove the temporary output file created during the test
            os.remove(test_output)


def cut_fasta(fasta_file):
    """
    Segments a FASTA file into individual sequences.

    Parameters:
    fasta_file (str): Path to the FASTA file.

    Returns:
    dict: A dictionary of sequences keyed by their IDs.
    """

    sequences = {}
    with open(fasta_file, 'r') as file:
        sequence_id = ''
        sequence = ''

        for line in file:
            # Check if the line is a header line (starts with '>')
            if line.startswith('>'):
                # Save the current sequence (except for the first header) before processing the next one
                if sequence_id:
                    sequences[sequence_id] = sequence

                # Extract the sequence ID from the header
                sequence_id = line.strip()[1:]
                # Reset sequence for the new header
                sequence = ''
            else:
                # Append the line to the sequence string
                sequence += line.strip()

        # Ensure the last sequence is also added to the dictionary
        sequences[sequence_id] = sequence

    return sequences


def run_patmatmotifs(sequence, seq_id, patmatmotifs_output):
    """
    Runs patmatmotifs for each sequence to scan for motifs.

    Parameters:
    sequence (str): The protein sequence.
    seq_id (str): Sequence identifier.
    patmatmotifs_output (str): Directory for saving motif scan results.

    Returns:
    None, but generates an output file for each sequence.
    """

    # Extract the ID part from the sequence descriptor
    seq_id_part = seq_id.split()[0]

    # Create a temporary file for the individual sequence
    temp_sequence_file = f'{seq_id_part}.fasta'
    with open(temp_sequence_file, 'w') as temp_file:
        temp_file.write(f'>{seq_id}\n{sequence}\n')

    # Define the output file path for patmatmotifs results
    output_file = os.path.join(patmatmotifs_output, f'patmatmotifs_output_{seq_id}.txt')

    # Run the patmatmotifs command
    command = ['patmatmotifs', '-sequence', temp_sequence_file, '-outfile', output_file]
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        print(f"\nError running patmatmotifs for sequence {seq_id}")

    # Clean up by removing the temporary file
    os.remove(temp_sequence_file)


def extract_motif_hits(output_folder, protein_family, taxonomic_group):
    """
    Extracts and summarizes motif hits from patmatmotifs output files.

    Parameters:
    output_folder (str): Directory containing patmatmotifs output files.
    protein_family (str): Protein family name.
    taxonomic_group (str): Taxonomic group.

    Returns:
    None, but creates a summary file of motif hits.
    """

    # Compose the summary file name and its path
    file_name = f'{protein_family}_in_{taxonomic_group}_motif_hits_summary.txt'
    # Replace characters that might cause file system issues
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    summary_file = os.path.join(output_folder, file_name)

    # Initialize counters for processed files and total file count
    processed_file = 0
    count = sum(1 for file_name in os.listdir(output_folder) if file_name.startswith('patmatmotifs_output_'))

    with open(summary_file, 'w') as summary:
        # Iterate over each file in the output folder
        for file_name in os.listdir(output_folder):
            if file_name.startswith('patmatmotifs_output_'):
                with open(os.path.join(output_folder, file_name), 'r') as file:
                    content = file.read()
                    # Check for motif hits in the file content
                    if "HitCount: 0" not in content:
                        summary.write(f"File: {file_name}\n")
                        # Find the start and end of the relevant sequence information
                        start = content.find("# Sequence:")
                        end = content.find("#--", start)
                        # Write the extracted information to the summary file
                        summary.write(content[start:end])
                        summary.write("----------------------------------\n\n")

                processed_file += 1
                # Update the progress bar after processing each file
                print_progress_bar('Summary', processed_file, count, bar_length=50)

    print(f"The relevant motif information has been saved to: {summary_file}.")


def print_progress_bar(progress, iteration, total, bar_length=50):
    """
    Prints a progress bar to the console.

    Parameters:
    progress (str): Description of the progress being tracked.
    iteration (int): Current iteration (or progress) of the task.
    total (int): Total iterations (or the completion point) of the task.
    bar_length (int, optional): Length of the progress bar in characters. Defaults to 50.

    Returns:
    None, but prints the progress bar to the console.
    """

    # Calculate the percentage completed
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    # Determine the length of the filled portion of the bar
    filled_length = int(bar_length * iteration // total)
    # Create the bar string with filled and unfilled portions
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    # Print the progress bar with the percentage completion
    sys.stdout.write(f'\r{progress}: |{bar}| {percent}% Complete')
    sys.stdout.flush()  # Ensure the output is displayed immediately

    # Print a newline character once the task reaches completion
    if iteration == total:
        print()  # Print the final newline to move to the next line in the console

# ↑ End: Scan motifs in Prosite

# ↓ Start: Blast analysis

def blast_analysis(protein_family, taxonomic_group, fasta_file):
    """
    Conducts BLAST analysis using a specified FASTA file as the database.
    It involves creating a BLAST database, choosing BLAST type, and running BLAST queries.

    Parameters:
    protein_family (str): Name of the protein family.
    taxonomic_group (str): Name of the taxonomic group.
    fasta_file (str): Path to the FASTA file for database creation.
    """

    print('#===========================')
    print(f"\nWe will use the previously obtained sequence: {fasta_file} as the blast database!")

    # Inform the user to have the sequence ready for BLAST analysis
    print(
        "Before conducting blast analysis, please ensure that you have obtained the sequence you want to use for blast analysis in advance!")

    # Confirm with the user to proceed with the current FASTA file
    while True:
        choice = input("Do you want to continue with it? (Y/N): ").strip().upper()
        if choice == 'Y':  # User chooses to continue with the current FASTA file
            print("Continuing!")
            print('#---------------------------')
            break
        elif choice == 'N':  # User chooses not to continue, exits the function
            print("Returning to Options Menu")
            print('#===========================')
            return
        else:  # User input is not recognized, prompting again
            print("Unable to recognize your input, please follow the prompts to enter!")

    # Ask the user for the output folder for BLAST analysis
    blast_output = ask_output_folder('Blast Analysis')

    # Create a BLAST database from the provided FASTA file
    blastdb = run_mkdb(protein_family, taxonomic_group, fasta_file, blast_output)

    # Get the BLAST query type from the user
    while True:
        blast_type = input("Please enter BLAST query type (blastp, blastx): ").strip()
        # Validating the input for BLAST type
        if blast_type != "blastp" and blast_type != "blastx":
            # Input is not a valid BLAST type, prompting the user again
            print("Unable to recognize your input, please follow the prompts to enter!")
            print('#---------------------------')
        else:  # Valid input, breaking out of the loop
            break

    # Run the BLAST analysis
    run_blast(protein_family, taxonomic_group, blast_type, blastdb, blast_output)
    print('#===========================')


def run_mkdb(protein_family, taxonomic_group, fasta_file, directory):
    """
    Creates a BLAST database from a FASTA file.

    Parameters:
    protein_family (str): The protein family name.
    taxonomic_group (str): The taxonomic group.
    fasta_file (str): Path to the FASTA file.
    directory (str): Directory where the BLAST database will be saved.

    Returns:
    Path to the created BLAST database.
    """

    # Define the database name
    file_name = f'{protein_family}_in_{taxonomic_group}_blast_db'
    # Replace characters that might cause file system issues
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    # Set the path for the BLAST database output
    blastdb_output = os.path.join(directory, file_name)

    try:
        # Informing the user that the BLAST database creation is in progress
        print(f'Makeblastdb is working for you to make {fasta_file} as a blast database.')
        # Running the makeblastdb command to create a BLAST database from the given FASTA file
        result = subprocess.run(['makeblastdb', '-in', fasta_file, '-dbtype', 'prot', '-out', blastdb_output],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        # Catch and display any errors during database creation
        print(f"Error occurred: {e}\n{e.stderr}")
    else:
        # Inform the user about the successful creation of the BLAST database
        print(f'Blast database has been saved to {blastdb_output}. You can check it now.')
        print('#---------------------------')

    return blastdb_output


def run_blast(protein_family, taxonomic_group, blast_type, db, directory):
    """
    Runs the BLAST analysis using the specified query type (blastp or blastx) against the created
    BLAST database. The results are saved in the provided directory.

    Parameters:
    protein_family (str): The protein family name.
    taxonomic_group (str): The taxonomic group.
    blast_type (str): Type of BLAST query (blastp or blastx).
    db (str): Path to the BLAST database.
    directory (str): Directory for saving the BLAST results.

    Returns:
    None, but saves the BLAST results to a file.
    """

    # Define the output file name for BLAST results
    file_name = f'{protein_family}_in_{taxonomic_group}_{blast_type}_result'
    # Replace characters that might cause file system issues
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    blast_output = os.path.join(directory, file_name)

    while True:
        # Prompting user to input the path of the query file
        query_file = input("Please enter the path to query the file: ").strip()

        # Check if the query file exists
        if not os.path.exists(query_file):
            print("The file does not exist. Please enter a valid file path.")
            continue

        # Verifying if the file is a valid FASTA file
        if not (query_file.endswith('.fasta') or query_file.endswith('.fa')) or not is_fasta(query_file):
            print("The file does not appear to be a FASTA file.")
            continue

        # Extract the first sequence for type checking
        with open(query_file, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    continue  # Skip header lines
                first_sequence = line.strip()
                break

        # Validating the first sequence against the expected sequence type for the BLAST query
        if blast_type == 'blastp' and not is_protein_sequence(first_sequence):
            print("The file does not contain valid protein sequences for blastp.")
            continue

        if blast_type == 'blastx' and not is_nucleotide_sequence(first_sequence):
            print("The file does not contain valid nucleotide sequences for blastx.")
            continue

        break

    try:
        # Run the BLAST command
        blast_command = [blast_type, '-query', query_file, '-db', db, '-out', blast_output]
        subprocess.run(blast_command, check=True)
    except subprocess.CalledProcessError as e:
        # Handle any errors during the BLAST run
        print(f"Error occurred: {e}")
    else:
        # Inform the user that the BLAST results have been saved
        print(f"{blast_type} result has been saved to {blast_output}. You can check it now.")
        print('#---------------------------')
        return


def is_fasta(filename):
    """
    Checks if the given file is in FASTA format.

    Parameters:
    filename (str): Path to the file to be checked.

    Returns:
    bool: True if the file is in FASTA format, False otherwise.
    """

    with open(filename, 'r') as file:
        # Read the first line of the file
        first_line = file.readline()
        # Check if the first line starts with '>'
        return first_line.startswith('>')

def is_protein_sequence(sequence):
    """
    Determines if a given sequence consists only of protein-specific amino acids.

    Parameters:
    sequence (str): The sequence to be checked.

    Returns:
    bool: True if the sequence is a valid protein sequence, False otherwise.
    """

    # Check if every character in the sequence is one of the 20 standard amino acids
    return all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in sequence.upper())

def is_nucleotide_sequence(sequence):
    """
    Verifies if a given sequence is composed of nucleotide bases.

    Parameters:
    sequence (str): The sequence to be checked.

    Returns:
    bool: True if the sequence is a valid nucleotide sequence, False otherwise.
    """

    # Checking if every character in the sequence is a valid nucleotide base
    return all(c in 'ACGTN' for c in sequence.upper())

# ↑ End: Blast analysis

# ↓ Start: View statistical data

def stats_data(protein_family, taxonomic_group, fasta_file):
    """
    Initiates the calculation of various protein properties using pepstats.

    Parameters:
    protein_family (str): Name of the protein family being analyzed.
    taxonomic_group (str): Taxonomic group to which the proteins belong.
    fasta_file (str): Path to the FASTA file containing protein sequences.

    """

    print('#===========================')
    print("\nWe will calculate statistics of protein properties！ e.g. Molecular weight, Number of residues, Average residue weight.")

    # Requesting the user to specify the output folder for saving pepstats results
    pepstats_output = ask_output_folder('Pepstats')

    # Running pepstats analysis with the specified FASTA file and output directory
    run_pepstats(protein_family, taxonomic_group, fasta_file, pepstats_output)


def run_pepstats(protein_family, taxonomic_group, fasta_file, directory):
    """
    Executes the pepstats command to calculate statistics of protein properties
    such as Molecular weight, Number of residues, and Average residue weight.

    Parameters:
    protein_family (str): Name of the protein family.
    taxonomic_group (str): Name of the taxonomic group.
    fasta_file (str): Path to the FASTA file containing protein sequences.
    directory (str): Directory where the pepstats output will be saved.

    Returns:
    None, but saves a file with the calculated protein statistics.
    """

    # Defining the output file name based on the protein family and taxonomic group
    file_name = f'{protein_family}_in_{taxonomic_group}.pepstats'
    # Replace characters that might cause file system issues
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    # Constructing the full path for the pepstats output file
    pepstats_output = os.path.join(directory, file_name)

    try:
        # Preparing and running the pepstats command with specified parameters
        command = ["pepstats", "-sequence", fasta_file, "-outfile", pepstats_output]
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        # Catching and reporting any errors encountered during the pepstats execution
        print(f"Error occurred: {e}")
    else:
        # Confirming successful completion and location of pepstats results
        print(f'Pepstats result has been saved to {pepstats_output}. You can check it now.')
        print('#===========================')

# ↑ End: View statistical data

# ↓ Start: Predict Tmap

def tmap(protein_family, taxonomic_group, fasta_file):
    """
    Initiates the prediction and plotting of transmembrane segments in protein sequences.

    Parameters:
    protein_family (str): Name of the protein family being analyzed.
    taxonomic_group (str): Taxonomic group to which the proteins belong.
    fasta_file (str): Path to the FASTA file containing protein sequences.

    """

    print('#===========================')
    # Inform the user about the upcoming transmembrane segment prediction task
    print("\nWe will predict and plot transmembrane segments in protein sequences.")

    # Asking the user for the output folder where Tmap results will be saved
    tmap_output = ask_output_folder('Tmap')

    # Running the Tmap analysis with the specified FASTA file and output directory
    run_tmap(protein_family, taxonomic_group, fasta_file, tmap_output)

def run_tmap(protein_family, taxonomic_group, fasta_file, directory):
    """
    Executes the Tmap command to predict and plot transmembrane segments in protein sequences.

    Parameters:
    protein_family (str): Name of the protein family.
    taxonomic_group (str): Name of the taxonomic group.
    fasta_file (str): Path to the FASTA file containing protein sequences.
    directory (str): Directory where the Tmap output will be saved.

    Returns:
    None, but saves a file with the predicted transmembrane segments and their graphical representation.
    """

    # Define the output file name based on the protein family and taxonomic group
    file_name = f'{protein_family}_in_{taxonomic_group}.tmap'
    # Replace characters that might cause file system issues
    file_name = file_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
    # Set the full path for saving the Tmap output
    tmap_output = os.path.join(directory, file_name)

    try:
        # Preparing and executing the Tmap command with specified parameters
        command = ['tmap', '-sequence', fasta_file, '-graph', 'png', '-goutfile', tmap_output, '-outfile', tmap_output]
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        # Catch and report any errors encountered during Tmap execution
        print(f"Error occurred: {e}")
    else:
        # Confirm the successful completion and location of Tmap results
        print(f'Tmap result has been saved to {tmap_output}. You can check it now.')
        print('#===========================')

# ↑ End: Predict Tmap

def main():
    """
    The main function of the script. It handles the user interface, allowing users to select various functionalities
    related to bioinformatics data analysis. These functionalities include changing datasets, conservation analysis,
    scanning with motifs, blast analysis, viewing statistical data, and predicting transmembrane segments.
    """

    print('\nWelcome to our sequence analysis program!'
          '\nIt is designed for bioinformatics data analysis.'
          '\nIt can provide you with many powerful sequence analysis tools.'
          '\nEnjoy it!')

    # Step 1: Get user input to identify the desired data and perform a basic analysis
    fasta_file, taxonomic_group, protein_family, = get_data()

    # ↓ Start: Options Menu

    # Step 2: Interactive interface for users to choose functionalities and loop back after each operation
    while True:
        print(f"\nWhat do you want to do next with your data in {fasta_file}?")
        print("1. Change dataset: \n\tChange the search criteria to get new data (subsequent outputs will all be based on the new dataset)")
        print("2. Conservation analysis: \n\tMultiple Sequence Alignment, Conditional Screening, Plotting Conserved Levels")
        print("3. Scan with motifs: \n\tScan protein sequence(s) of interest with motifs from the PROSITE database")
        print("4. Blast analysis: \n\tUsing the dataset as a blast database, perform blast analysis on the specified sequences")
        print("5. View statistical data: \n\tCalculate statistics of protein properties")
        print("6. Tmap: \n\tPredict and plot transmembrane segments")
        print("7. Exit program")

        choice = input('\nPlease enter the number before the option (e.g., 1, 2, 3): ')

        if choice == '1':
            # Allows user to change the search criteria and fetch new data
            fasta_file, taxonomic_group, protein_family = get_data()
        elif choice == '2':
            # Perform conservation analysis including multiple sequence alignment and plotting conservation levels
            aligned_file = conservation_analysis(protein_family, taxonomic_group, fasta_file)
        elif choice == '3':
            # Scan for motifs in protein sequences using the PROSITE database
            scan_prosite_motifs(protein_family, taxonomic_group, fasta_file)
        elif choice == '4':
            # Perform BLAST analysis using the dataset as a database
            blast_analysis(protein_family, taxonomic_group, fasta_file)
        elif choice == '5':
            # Calculate and display statistics related to protein properties
            stats_data(protein_family, taxonomic_group, fasta_file)
        elif choice == '6':
            # Predict and visualize transmembrane segments in proteins
            tmap(protein_family, taxonomic_group, fasta_file)
        elif choice == '7':
            # Exit the program with a closing message
            print("Exiting the program."
                  "\nThank you for your using!")
            print('#===========================')
            exit()
        else:
            # Handle invalid input and prompt the user to enter a valid option
            print("Unable to recognize your input, please follow the prompts to enter!")

    # ↑ End: Options Menu

if __name__ == "__main__":
    main()
