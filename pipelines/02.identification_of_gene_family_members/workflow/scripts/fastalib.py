
#!/usr/bin python3

import os
from natsort import natsorted
from xopen import xopen
import argparse
import re
from collections import defaultdict


def isfile ( input_file, delimiter = '\t' ):
    try:
        with open ( input_file, "r" ) as file:
            lines = file.readlines()
            
            instances = []
            
            for line in lines:
                instances.append(line)

    except FileNotFoundError:    
        instances = sequence

    return instances



class FASTA:

    __slots__ = ['id','seq']
   
    def __init__ (self, header = "", sequence = ""):
        self.id = str(header)
        self.seq = str(sequence)
        

    def __str__ ( self ) :
        return f">{self.id}\n{self.seq}"
    

    def fasta_parser(fasta_file, sprint = False) -> list:

        '''

        Parse a FASTA file and return a list of FASTA instances.

        Parameters:
        - fasta_file (str): Path to the FASTA file.
        - sprint (bool): If True, print the parsed sequences.

        Returns:
        - list: List of FASTA instances.

        '''      
        
        fasta_instances = []

        with xopen(fasta_file, 'r') as fasta:
            lines = fasta.readlines()
            seqid = ""
            sequence = ""
            entries = defaultdict()
            encountered_ids = set()

            for line in lines:
                line = line.strip('\n')  # Remove newline character

                if line.startswith('>'):
                    seqid = re.sub("\s+", "_", line[1:])
                    seqid = re.sub("\r", "\n", line[1:])
                    sequence = ""

                    if seqid in encountered_ids:
                        print(f"Error: Duplicate seqid found - '{seqid}'.")

                    else:
                        encountered_ids.add(seqid)

                    if len(seqid) > 50: # Check if seqid length is greater than 50
                        print(f"Error: Seqid '{seqid}' has length greater than 50.")
                            
                elif not re.search("--|#", line):  # Ignore lines containing "--" or "#"
                    sequence += re.sub("\*$|\.$", "", re.sub("\W", "", line))

                    if seqid and sequence:
                        entries[seqid] = sequence

        for seqids, sequences in entries.items():
            fasta_instance = FASTA(seqids, sequences)
            fasta_instances.append(fasta_instance)

        return fasta_instances
    

    def output_one_by_one ( fasta_file ):

        '''
        Option to output each sequence one-by-one in a new FASTA file.

        '''
        
        fasta_instances = FASTA.fasta_parser ( fasta_file )

        for fasta_instance in fasta_instances:
            with open ( fasta_instance.id, "w" ) as out:
                out.write ( str(fasta_instance) )

    
    def fasta_sizer(fasta_file, sequence=""):
        lst = []
        selected_lst = []
        fasta_instances = FASTA.fasta_parser(fasta_file)

        if sequence is not None and not sequence == "":
            sequence_p = re.compile(sequence)
        else:
            sequence_p = None

        for fasta_instance in fasta_instances:

            if sequence_p is not None:
                if re.search(sequence_p, fasta_instance.id):
                    selected_fasta_size = '\t'.join([fasta_instance.id, str(len(fasta_instance.seq))])
                    selected_lst.append(selected_fasta_size)
                else:
                    print(f"The requested {sequence} does not exist in the FASTA file. Please check your input.")
                    return
            else:
                fasta_size = '\t'.join([fasta_instance.id, str(len(fasta_instance.seq))])
                lst.append(fasta_size)

        if len(selected_lst) > 0:
            sorted_list = natsorted(selected_lst, key=lambda x: x.split('\t')[1])
            print('Printing sizes of selected sequences in the provided FASTA file.')
            return sorted_list

        elif len(lst) > 0:
            sorted_list = natsorted(lst, key=lambda x: x.split('\t')[1])
            print('Printing sizes of all sequences in the provided FASTA file.')
            return sorted_list

       
    def replace_names ( fasta_file, sequence="", ncbi_tsa_submission = False ) -> list:

        '''
        Replaces original IDs of fasta sequences with new associated IDs.
        
        Parameters:
        - fasta_file (str): Path to the FASTA file.
        - sequence (str): Path to the list of user-specified sequence IDs.

        Returns:
        - fasta: Fasta file with replace IDs for specified sequences.
        
        '''
       
        fasta_instances = FASTA.fasta_parser ( fasta_file )

        new_fasta_instances = []
        
        for fasta_instance in fasta_instances:

            if ncbi_tsa_submission == True:
                fasta_instance.id = fasta_instance.id + ' [moltype=transcribed_RNA]'
                new_fasta_instances.append(fasta_instance)

            elif sequence:
                for g in isfile( sequence ):
                    if re.search ( fasta_instance.id , g ):
                        fasta_instance.id = g
                        new_fasta_instances.append(fasta_instance)

        return new_fasta_instances

    
    def translate ( fasta_file ) -> list:

        '''
        Translate nucleotide sequences to protein sequences.

        Parameters:
        - fasta_file (str): Path to the FASTA CDS file.

        Returns:
        - list: List of FASTA instances with translated protein sequences.
        
        '''
        
        codon2aa = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
            'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W'
        }

       
        fasta_instances = FASTA.fasta_parser(fasta_file)

        protein_fasta = []
        
        for fasta_instance in fasta_instances:
            fasta_instance.seq = fasta_instance.seq.upper() ## Convert CDS to uppercase

            if not re.search ("A|T|G|C", fasta_instance.seq):
                print ( fasta_instance.id + " is not a nucleotide sequence. Found other than A,T,G,C." )

            else:
                aminoacids = []
                
                for i in range(0,len(fasta_instance.seq),3):
                    codon = fasta_instance.seq[i:i+3]
                    
                    if codon in codon2aa:
                        aminoacids.append(codon2aa[codon])

                    else:
                        aminoacids.append("N")

            protein_sequence = FASTA ( fasta_instance.id, ''.join(aminoacids) )   
            protein_fasta.append( protein_sequence )

        if len(protein_fasta ) > 0 :
            return protein_fasta

        else:
            print ( "Provided fasta file is empty.")

            
    def extract_subsequences (fasta_file, sequence, start_position="", end_position="", extract = True, ncbi_tsa_submission = False) -> list:

        '''
        Extracts or removes sequence(s) and subsequence(s) from a fasta file.
        Checks each sequence instance independently and if start or end position for removing the sequences or subsequences are not provided, it will remove the entire sequence.

        Parameters:
          - fasta_file (str): Path to the FASTA file.
          - sequence (str): Name of file with provided sequence(s).
          - start_position(int): Start position in the sequence. If no start_position is provided, 1 will be used as default.
          - end_position(int): End position in the sequence. If no start_position is provided, length of sequence will be used as default.
          - extract(bool): Extract or remove subsequence(s). Default: extract=True
          - ncbi_tsa_submission(bool): Option to make checks for NCBI TSA submission.
 
        Returns:
          - list: List of sequences based on options to extract or remove entire sequence(s) and subsequence(s).

        '''
        
        fasta_instances = FASTA.fasta_parser(fasta_file)
        matches = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
        intact_sequences = []
        intact_sequences_r = []
        subsequences = []
        start_pos = int()
        end_pos = int()
        
        for fasta_instance in fasta_instances: ## Find matches between the fasta file and the provided sequence list, which is either a file or a string, and add them into a dictionary
            if os.path.isfile(sequence):
                for sequence_file in isfile(sequence):                   
                    if len(sequence_file.split('\t')) == 3:
                        if re.search(sequence_file.strip('\n').split('\t')[0], fasta_instance.id):
                            matches[fasta_instance.id][sequence_file.split('\t')[1]][sequence_file.split('\t')[2]] = fasta_instance.seq
                    elif len(sequence_file.strip('\n').split('\t')) == 1:
                        if re.search( sequence_file.strip('\n').split('\t')[0], fasta_instance.id):
                            matches[fasta_instance.id][1][len(fasta_instance.seq)] = fasta_instance.seq
                    else:
                        print (f"Please check your file {sequence_file} format.")
                        
            elif isinstance(sequence, str): ## Provided sequence list is a string, NOT a file
                if re.search(sequence, fasta_instance.id):
                    if start_position == "": ## If no start and/or end positions are provided, return the entire sequence
                        start_position = 1
                    if end_position == "":
                        end_position = len(fasta_instance.seq)
                    matches[fasta_instance.id][start_position][end_position] = fasta_instance.seq
            else:
                print ('Provided file is neither a FILE nor a STRING. Please check your input.')
                return fasta_instances

        if len(matches) > 0: ## If matches are found, extract or remove entire sequence(s) or subsequence(s) based on user-provided specifications
            intact_sequences = [fasta_instance for fasta_instance in fasta_instances if fasta_instance.id not in matches.keys()] ## Creates a list with all sequences NOT matching the sequence list ONCE. Isolates the sequences we want to remove from the file.
                    
            for matching_seq, positions in matches.items(): ## Create dictionaries to save the fasta ID, start, end positions and sequence of matching sequences.
                if isinstance(positions, dict):                   
                    for start_positions, end_positions in positions.items():
                        start_position = int(start_positions)

                    for end_positions, sequences in end_positions.items():
                        end_position = int(end_positions)
                        sequence_f = str(sequences)
        
                if start_position and end_position and sequence_f: 
                    if 1 <= start_position <= len(sequence_f) and 1 <= end_position <= len(sequence_f):
                        if extract:  ## Options to extract sequence(s) and subsequence(s)
                            subsequence = sequence[start_position:end_position]                          
                            subsequences.append( FASTA(matching_seq, subsequence) )

                            print ( f"Extracted {start_position} - {end_position} from {matching_seq}" )

                        else: ## Options to remove sequence(s) and subsequence(s)

                            if start_position == 1 and end_position == len(sequence_f):
                                print ( f"Removed entire sequence {matching_seq}" )
                                continue

                            else:

                                if not start_position == 1 and not end_position == len(sequence_f): ## If start and end position are in the middle of the sequence, collapse the subsequence and merge the flanking regions.
                                    subsequence = sequence[:start_position-1] + sequence_f[end_position+1:]
                                    intact_sequences_r.append( FASTA(matching_seq, subsequence) )

                                elif start_position == 1 and not end_position == len(sequence_f): ## If start position is 1 and end position is in the middle, keep sequence from start to the provided end position
                                    subsequence = sequence_f[end_position+1:]
                                    intact_sequences_r.append( FASTA(matching_seq, subsequence) )
                                    
                                elif not start_position == 1 and end_position == len(sequence_f): ## If start is in the middle and end position is last position in sequence, keep only up to the start
                                    subsequence = sequence_f[:start_position-1]
                                    intact_sequences_r.append( FASTA(matching_seq, subsequence) )

                                print ( f"Removed {start_position} - {end_position} from {matching_seq}" )
                                
                    else:
                        print(f"Requested range {start_position} - {end_position} not present in {matching_seq}")
                else:
                    print (f"Start and end positions not provided for {matching_seq}.")
        else:
            print ('No matching sequences found in fasta file.')

            if ncbi_tsa_submission == True:
                fasta_sequences = []
                 
                for fasta_instance in fasta_instances:
                    if len(fasta_instance.seq) >= 200:
                        fasta_sequences.append(fasta_instance)
                return fasta_sequences

        
        ## Size filtering here, if ncbi_tsa_submission is enabled to avoid doing multiple times in the loop

        extracted_sequences = []
        kept_sequences = []
        
        if subsequences: ## Extract option enabled and subsequence(s) were extracted
            
            if ncbi_tsa_submission == True:                
                for fasta_sequence in subsequences:
                    if len(fasta_sequence.seq) >= 200:
                        extracted_sequences.append(fasta_sequence)
                        return extracted_sequences
                    
            else:    
                return subsequences

            
        elif len(intact_sequences_r) > 0 and not subsequences: ## Remove option enabled and subsequence(s) were removed           
            intact_sequences_r = intact_sequences_r + intact_sequences

            if ncbi_tsa_submission == True:
                for fasta_sequence in intact_sequences_r:
                    if len(fasta_sequence.seq) >= 200:
                        kept_sequences.append(fasta_sequence)
                return kept_sequences
            
            else:
                return intact_sequences_r
                        
        
        elif intact_sequences and not subsequences and not intact_sequences_r: ## Remove option enabled for entire sequence(s)
            
            if ncbi_tsa_submission == True:
                for fasta_sequence in intact_sequences:
                    if len(fasta_sequence.seq) >= 200:
                        kept_sequences.append(fasta_sequence)
                return kept_sequences
                
            else:
                return intact_sequences

        

    def check_position (fasta_file, sequence, position_number, length = 0) -> str:

        '''
        Returns the sequence from a fasta file, based on position number.

        Parameters:
          - fasta_file (str): Path to the FASTA file.
          - sequence (str): Name of the sequence.
          - position_number (int): Position number in the sequence.
          - length (int): Optional length parameter.

        Returns:
          - str: Subsequence based on the provided position and length.

        '''
        
        fasta_instances = FASTA.fasta_parser(fasta_file)

        sequence = None
        sequences = []

        for fasta_instance in fasta_instances:
            if re.search(sequence, fasta_instance.id):
                sequence = fasta_instance
                break

        if sequence:
            position_number -= 1 ## counting starts by default at zero
            
            if not position_number > len(sequence.seq) :

                if length > 0:
                    end_position = position_number + length
                    sequences.append(sequence.seq[position_number:end_position])
                    return ''.join(sequences)

                elif length < 0:
                    start_position = position_number - length
                    sequences.append(sequence.seq[start_position:position_number])
                    return ''.join(sequences)
                
                else:
                    return sequence.seq[position_number + length]
                
            else:
                print ( f"Requested position number not present in {sequence.id}" )
                return None
            
        else:
            print ( f"Requested sequence {sequence} does not exist in {fasta_file}" )
            return None
        

        
## Implementation as a main script ##

def parse_arguments():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating fasta files.')

    ## options to parse and convert multi-line fasta to one-line fasta
    parser.add_argument('-f','--fasta', type=str, help='FASTA file')
    parser.add_argument('-o','--one_line', action="store_true", help='Option to extract FASTA sequences.')

    ## option to compute size of fasta sequences
    parser.add_argument('-s','--size', action="store_true", help='Option to print sizes of FASTA sequences.')

    ## options to extract or remove sequences
    parser.add_argument('-seq','--sequence', type=str, help='User-provided string or list of sequence IDs to extract, remove or replace with new in FASTA file.')
    parser.add_argument('-obo','--one_by_one', action="store_true", help='Print all fasta sequences in individual files one-by-one.')
    parser.add_argument('-e','--extract', action="store_true", help='Option to extract FASTA sequences.')
    parser.add_argument('-r','--remove', action="store_true", help='Option to remove FASTA sequences.')
    parser.add_argument('-st','--start_position', type=str, required = False, help='Start position to extract FASTA sequences.')
    parser.add_argument('-end','--end_position', type=str, required = False, help='End position to extract FASTA sequences.')
    parser.add_argument('-ncbi','--ncbi_tsa_submission', action="store_true", help='Option to perform checks for NCBI TSA submissions.')

    ## replace sequence names; original sequence ids should be associated with new names
    parser.add_argument('-nn', '--new_names', action="store_true", help='Option to replace names/IDs of FASTA sequences.')

    ## option to convert CDS to PEP 
    parser.add_argument('-tr','--cds2pep', action="store_true", help='Option to print sizes of FASTA sequences.')

    ## get sequence based on position number
    parser.add_argument('-pos','--position', action="store_true", help='Option to return requested position number of FASTA sequence.')
    parser.add_argument('-num','--number', type=int, help='Position number to requested sequence from a FASTA sequence.')
    parser.add_argument('-len','--length', type=int, help='Length to add to position number to requested sequence from a FASTA sequence. e.g. 280 + length')

    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")

    return parser.parse_args()


def main():
    
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating fasta files.')
    args = parse_arguments()

    threads = []

    if args.fasta:

        if re.search(".fsa", args.fasta):
            if args.ncbi_tsa_submission:
                print (f"Renaming {args.fasta} file because .fsa suffix already exists and would cause an error!")
                os.system(f"rename 's/.fsa/.fasta/' {args.fasta}")
            
        inp = re.sub (".aa$|.fa$|.faa$|.fna$|.fsa|.fasta$|.1l$","",args.fasta)

        if args.one_line:
            with open(f"{inp}.fasta.1l", "w") as f:
                for out in FASTA.fasta_parser(args.fasta):
                    print ( out, file = f )
                    
        elif args.one_by_one:
            FASTA.output_one_by_one(args.fasta)

        elif args.size:
            if args.sequence:
                print ( FASTA.fasta_sizer ( args.fasta, sequence = args.sequence ) )
            else:
                print ( "No specific sequence ID was provided. Will print the sizes of all sequences in FASTA file.")
                
                with open(f"{inp}.fasta.sizes", "w") as f:
                    for out in FASTA.fasta_sizer(args.fasta, sequence = args.sequence):
                        print ( out, file = f )

        elif args.extract:
            if args.sequence:
                if args.ncbi_tsa_submission:
                    with open(f"{inp}.fsa", "w") as f:
                        for out in FASTA.extract_subsequences(fasta_file = args.fasta, sequence = args.sequence, ncbi_tsa_submission = True, extract = True):
                            print ( out, file = f )
                else:
                    with open(f"{inp}.extracted.{args.sequence}.fasta", "w") as f:
                        for out in FASTA.extract_subsequences(fasta_file = args.fasta, sequence = args.sequence, extract = True):
                            print ( out, file = f )
            else:
                print ( "Please provide a sequence ID or list of sequence IDs to extract (sub)sequence(s).")
 
        elif args.remove:
            if args.sequence:
                if args.ncbi_tsa_submission:
                    with open(f"{inp}.fsa", "w") as f:
                        for out in FASTA.extract_subsequences(fasta_file = args.fasta, sequence = args.sequence, ncbi_tsa_submission = True, extract = False):
                            print ( out, file = f )
                else:
                    with open(f"{inp}.removed.{args.sequence}.fasta", "w") as f:
                        for out in FASTA.extract_subsequences(fasta_file = args.fasta, sequence = args.sequence, extract = False):
                            print ( out, file = f )
            else:
                print ( "Please provide a sequence ID or list of sequence IDs to remove (sub)sequence(s).")
                
        elif args.new_names:
            if args.ncbi_tsa_submission:
                with open(f"{inp}.ncbi_tsa_submission.fsa", "w") as f:
                    for out in FASTA.replace_names(args.fasta, ncbi_tsa_submission = True):
                        print ( out, file = f )

            elif args.sequence:
                with open(f"{inp}.new_names.fasta", "w") as f:
                    for out in FASTA.replace_names(args.fasta, args.sequence):
                        print ( out, file = f )
                        
            else:
                print ( 'Please select a valid option for renaming fasta sequences: --ncbi_tsa_submission or --sequence.')        
                        
        elif args.cds2pep:
            with open(f"{inp}.pep.fasta", "w") as f:
                for out in FASTA.translate(args.fasta):
                    print ( out, file = f )

        elif args.position:
            if args.sequence:
                if args.number >=0:
                    if FASTA.check_position ( args.fasta, args.sequence, args.number, args.length if args.length else 0 ):
                        print ( "Requested site number", args.number, "in", args.sequence, "is", args.end==" " )
                        print ( FASTA.check_position ( args.fasta, args.sequence, args.number, args.length if args.length else 0 ) )
                        
                    else:
                        print ("Please provide a position number to return.")

                else:
                    print ( "Please provide the name of the sequence to search for returning the position number.")                                                
        else:
            print("Please select a potential FASTA operation.")
    else:
        print ( "Please provide the input fasta file.")

if __name__ == "__main__":
    main()
