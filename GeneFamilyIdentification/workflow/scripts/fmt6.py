
#!/usr/bin/env python3

import argparse
import re
from natsort import natsorted
import os 

class FMT6:

    """
    Custom class object for outfmt6 blast output format. 
    
    Includes functions for:
    1. Parsing and sorting fmt6 files.
    2. Extracting user-specified number of hits (from best to worst) from fmt6 file
    3. Filtering fmt6 files based on identity cutoff.
    4. Extracting user-specified hits from fmt6 files.
    5. Converting fmt6 to gff3 file format.
    6. Converting fmt6 to bed file format.

    """
    

    __slots__ = [ 'query', 'subject', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore' ]
    

    def __init__(self, qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore):
        self.query = qseqid
        self.subject = sseqid
        self.pident = pident
        self.length = length
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.evalue = evalue
        self.bitscore = bitscore


    def __str__ ( self ):
        return f"{self.query}\t{self.subject}\t{self.pident}\t{self.length}\t{self.mismatch}\t{self.gapopen}\t{self.qstart}\t{self.qend}\t{self.sstart}\t{self.send}\t{self.evalue}\t{self.bitscore}"


    @staticmethod
    def parse_fmt6 ( input_file ): ## opens fmt6 file format and creates new instances of FMT6 class
        fmt6_instances = []

        with open ( input_file, "r" ) as file:
            lines = file.readlines()
        
            for line in lines:
                line = line.strip('\n')              
                cols = line.split ('\t')
                fmt6_instance = FMT6( cols[0], cols[1], float(cols[2]), int(cols[3]), int(cols[4]), int(cols[5]), int(cols[6]), int(cols[7]), int(cols[8]), int(cols[9]), cols[10], float(cols[11]) )
                fmt6_instances.append(fmt6_instance)
            
        return fmt6_instances 

    
    def sort_fmt6(input, field = None):
        fmt6_instances = FMT6.parse_fmt6 ( input )
       
        if field == None:
            def custom_sort_key(element):
                return (element.subject, element.sstart, element.send)                   

            sorted_list = natsorted(fmt6_instances, key=custom_sort_key)

        else:
            sorted_list = natsorted( fmt6_instances, key=lambda instance: getattr(instance, field) )

        return sorted_list


    def hits ( fmt6_file, nhits ):
        seen = []
        hits = []
        fmt6_instances = FMT6.parse_fmt6 ( fmt6_file )
       
        for fmt6_instance in fmt6_instances:
            i = 0

            if not fmt6_instance.query in seen:
                seen.append ( fmt6_instance.query )
                i += 1

                if i <= nhits:
                    hits.append(fmt6_instance)
                    
        hits = list(set(hits))

        return hits


    def filter_fmt6(fmt6_file, cutoff):

        """
        Parse fmt6 format, filter based on identity cutoff, and sort based on the selected field.

        """

        lst = []
        fmt6_instances = FMT6.sort_fmt6(fmt6_file)
        
        for fmt6_instance in fmt6_instances:
            if fmt6_instance.pident >= float(cutoff):
                lst.append(fmt6_instance)

        return lst


    def extract_genes ( fmt6_file, sequences, query = True   ):

        hits = []
        
        glist = isfile(sequences)
        fmt6_instances = FMT6.parse_fmt6 ( fmt6_file )

        for fmt6_instance in fmt6_instances:
            if query == True:
                if fmt6_instance.query in glist:
                    hits.append(fmt6_instance)
            else:
                if fmt6_instance.subject in glist:
                    hits.append(fmt6_instance)
                    

        return hits      
      

    def blast2gff ( fmt6_file, gff_file, cutoff ):
        fmt6_instances = FMT6.sort_fmt6(fmt6_file)  # Get the list of instances from sort_fmt6

        for fmt6_instance in fmt6_instances:

            out = str()
            
            source = "BLAST"
            stype = "GENE"
            phase = "1"
            score = "."

            if fmt6_instance.pident >= cutoff:               

                if fmt6_instance.sstart < fmt6_instance.send:
                    fmt6_instance.strand = "+"
                    out = '\t'.join(["ID="+fmt6_instance.subject, source, stype, fmt6_instance.sstart, fmt6_instance.send, phase, fmt6_instance.strand, score, "ID="+fmt6_instance.query] )
                                
                else:
                    fmt6_instance.strand = "-"
                    out = '\t'.join(["ID="+fmt6_instance.subject, source, stype, fmt6_instance.send, fmt6_instance.sstart, phase, fmt6_instance.strand, score, "ID="+fmt6_instance.query] )

                write_file ( out, args.fmt6 )


    def blast2bed ( fmt6_file, bed_file, cutoff): ## opens fmt6 file, creates new instances of the fmt6 class and outputs a bed file
        fmt6_instances = FMT6.sort_fmt6(fmt6_file)  # Get the list of instances from sort_fmt6
        beds = []
        out = str()

        for fmt6_instance in fmt6_instances:
            
            source = "BLAST"
            stype = "GENE"
            phase = "1"
            score = "."

            if fmt6_instance.pident >= cutoff:

                if fmt6_instance.sstart < fmt6_instance.send:
                    fmt6_instance.strand = "+"
                    out = '\t'.join(["ID="+fmt6_instance.subject, fmt6_instance.sstart, fmt6_instance.send, "ID="+fmt6_instance.query, score, fmt6_instance.strand] )
                                
                else:
                    fmt6_instance.strand = "-"
                    out = '\t'.join(["ID="+fmt6_instance.subject, fmt6_instance.send, fmt6_instance.sstart, "ID="+fmt6_instance.query, score, fmt6_instance.strand] )

                beds.append(out)

            return beds



## Implementation ##

def parse_arguments():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating blast outfmt6 files.')
    parser = argparse.ArgumentParser()

    ## options for parsing, sorting and filtering fmt6 output files
    parser.add_argument('--fmt6', type=str, help='FMT6 input file')
    parser.add_argument('-s','--sort', action="store_true", help='Option to sort the FMT6 file.')
    parser.add_argument('-c','--cutoff', type=str, help='Identity cutoff to filter fmt6 results.')
    parser.add_argument('--field', type=str, help='Field used for sorting FMT6 file.')
    parser.add_argument('-nh','--nhits', type=int, help='Number of blast hits to keep in reverse order.')

    ## options for extracting gene list
    parser.add_argument('-e','--extract', action="store_true", help='Option to extract range from FMT6 file.')
    parser.add_argument('-s', '--sequences', type=str, help='List with gene IDs that you would like to extract.')
    parser.add_argument('-q','--query', type=str, help='Compare gene list with query or subject names of FMT6 file.')
        
    parser.add_argument('-r','--range', action="store_true", help='Option to extract range from FMT6 file.')
    parser.add_argument('-chr','--chromosome', type=str, help='Chromosome whose range you would like to extract.')
    parser.add_argument('-st','--start', type=str, help='Start position of the range you would like to extract.')
    parser.add_argument('-end','--end', type=str, help='End position of the range you would like to extract.')

    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")

    return parser.parse_args()



def main():
    
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating fmt6 files.')
    args = parse_arguments()
    
    if args.fmt6:
        
        inp = re.sub (".fmt6$|.b1$","",args.fmt6)
        
        if args.sort:
            with open(f"{inp}.{args.field}.sorted.fmt6", "w") as f:
                print ( "Sorting", args.fmt6, "based on", args.field,"!" )
                for out in FMT6.sort_fmt6(args.fmt6, args.field):
                    print ( out, file = f )

        elif args.cutoff:
            with open(f"{inp}.{args.cutoff}.fmt6", "w") as f:
                print ( "Filtering", args.fmt6, "with", args.cutoff+"% identity cutoff!" )
                for out in FMT6.filter_fmt6(args.fmt6, args.cutoff):
                    print ( out, file = f )

        elif args.nhits:
            with open(f"{args.fmt6}.b{args.nhits}", "w") as f:
                print ( "Keeping", args.nhits, "from", args.fmt6,"!" )
                for out in FMT6.hits(args.fmt6, args.nhits):
                    print ( out, file = f )

        elif args.extract:
            if args.sequences:
                with open(f"{inp}.{args.sequences}.fmt6", "w") as f:
                    print ( "Extracting", args.sequences, "from", args.fmt6 )

                    for out in FMT6.extract_genes(args.fmt6, args.sequences, query=args.query):
                        print ( out, file = f )
            else:
                print ( "Please provide a gene list." )

        else:
            print("Please provide an option.")

    else:
        print ("Please provide a fmt6 file as input.")

            
            
if __name__ == "__main__":
    main()

