
import os
from collections import defaultdict

def trinity_samples_file(samples, data_dir):
    species = defaultdict(list)

    for sample in samples:
        line = f"{sample}\t{os.path.join(data_dir, f'{sample}_1.fastq.gz')}\t{os.path.join(data_dir, f'{sample}_2.fastq.gz')}"
        species[sample[:-1]].append(line)

    for sp, lines in species.items():
        with open(f"{sp}.samples", "w") as outfile:
            outfile.write("\n".join(sorted(set(lines))))

data_directory = config["data_directory"]
samples = [str(file) for file in sorted([f[:-11] for f in os.listdir(data_directory) if f.endswith(".fastq.gz")])]
trinity_samples_file(samples, data_directory)

species = samples[:-3]


#============================================================== PIPELINE ============================================================#

rule all:
     input: 'evigene/trinity_unfiltered_spadeshard/okayset/{species}.trinityunfiltered.spadeshard.okay.tr'

rule make_dir:
     output: directory({species})
     shell: "mkdir {species}"
                    
rule trinity:
     input: species_dir = rules.make_dir.output, 
            species_samples = '{species}.samples'

     output: trinity_fasta = '{species}/trinity_out_dir/trinity_out_dir.Trinity.fasta'
     
     threads: config['threads']

     shell: """ trinity \
     	    	--seqType fq \
		--SS_lib_type RF \
		--max_memory 350G \
		--samples_file {input.species_samples} \
		--CPU {threads} \
		--trimmomatic """

     
rule rnaspades:
     input: '{species}/trinity_out_dir/trinity_out_dir.Trinity.fasta'
                       
     output: '{species}_spadesrna_outdir/hard_filtered_transcripts.fasta'

     threads: config['rnaspades_threads']

     shell: """ cat {species}/trinity_out_dir/*_1.fastq.gz.PwU.qtrim.fq > {species}_1.all.trimmed.fastq && \
       	    	cat {species}/trinity_out_dir/*_2.fastq.gz.PwU.qtrim.fq > {species}_2.all.trimmed.fastq && \

                rnaspades.py --ss rf  \
                             -1 {species}_1.all.trimmed.fastq  \
                             -2 {species}_2.all.trimmed.fastq  \
                             --checkpoints all  \
                             --threads {threads}  \
                             -o {species}_spadesrna_outdir """

rule evigene:
     input:
         trinity_transcripts = '{species}/trinity_out_dir/trinity_out_dir.Trinity.fasta',
         spades_hard_filtered ='{species}_spadesrna_outdir/hard_filtered_transcripts.fasta'

     output:
         okay_transcripts= 'evigene/trinity_unfiltered_spadeshard/okayset/{species}.trinityunfiltered.spadeshard.okay.tr'

     threads: config['evigene_threads']

     shell: """ cat {input.trinity_transcripts} {input.spades_hard_filtered} > {species}.trinityunfiltered.spadeshard.fasta && \
                perl ~/Programs/evigene/scripts/prot/tr2aacds.pl -cdnaseq {species}.trinityunfiltered.spadeshard.fasta -NCPU {threads} -logfile """

                       
rule transdecoder:
     input: '{species}.trinityunfiltered.spadeshard.fasta'

     output:

     shell: """ cd $i/functional_annotation  /data/iasonas/Programs/TransDecoder/TransDecoder.LongOrfs -t ../evigene/trinity_unfiltered_spadeshard/okayset/trinityunfiltered.spadeshard.okay.tr  \
     	    	
     diamond blastp --query trinityunfiltered.spadeshard.okay.tr.transdecoder_dir/longest_orfs.pep --db /data/panos/db/swissprot_Metazoa/swissprot_Metazoa.dmnd --outfmt 6 --out $i.longest_orfs.trinityunfiltered.spadeshard.okay.transdecoder.pep.vs.uniref50.fmt6 --more-sensitive --threads 20 && 
		hmmsearch --cpu 30 --domtblout $i.longest.orfs.pfam.domtblout /data/iasonas/db/Pfam-A.hmm trinityunfiltered.spadeshard.okay.tr.transdecoder_dir/longest_orfs.pep
		
    TransDecoder.Predict -t ../evigene/trinity_unfiltered_spadeshard/okayset/trinityunfiltered.spadeshard.okay.tr --retain_pfam_hits $i.longest.orfs.pfam.domtblout --retain_blastp_hits $i.longest_orfs.trinityunfiltered.spadeshard.okay.transdecoder.pep.vs.uniref50.fmt6
