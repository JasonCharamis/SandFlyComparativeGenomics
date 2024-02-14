species = { f[:-21] for f in os.listdir("/home/iasonas/genomics/sandflies/sandfly_transcriptomes/functional_annotation/genesets/") if f.endswith(".final_genesets.pep.
f") }

rule all:
     input:
        expand('/home/iasonas/genomics/sandflies/sandfly_transcriptomes/functional_annotation/genesets/{species}.dmnd', species = species),
	expand('papatasi.genomic.CCEs.VS.{species}.fmt6', species = species ),
	expand ('papatasi.genomic.CCEs.VS.{species}.fmt6.b1', species = species ),
	expand ("{species}.blast.pfam.id", species = species ),
	expand ('{species}.CCEs.pep.fasta', species = species )
	
rule makedb:
     input: '/home/iasonas/genomics/sandflies/sandfly_transcriptomes/functional_annotation/genesets/{species}.final_genesets.pep.f'
     output: '/home/iasonas/genomics/sandflies/sandfly_transcriptomes/functional_annotation/genesets/{species}.dmnd'
     shell: " diamond makedb --in {input} --db {output} "

rule diamond:
      input: '/home/iasonas/genomics/sandflies/sandfly_transcriptomes/functional_annotation/genesets/{species}.dmnd'
      output: 'papatasi.genomic.CCEs.VS.{species}.fmt6'
      shell: " diamond blastp --query ../papatasi.genomic.CCEs.pep.fasta --db {input} --out {output} --ultra-sensitive "

rule best_hit:
     input: 'papatasi.genomic.CCEs.VS.{species}.fmt6'
     output: 'papatasi.genomic.CCEs.VS.{species}.fmt6.b1'
     shell: " perl /home/iasonas/bin/parse_m8_blast_keep_kbest_hsp_per_query.pl {input} > {output} "

rule common_with_pfam:
     input: blast = 'papatasi.genomic.CCEs.VS.{species}.fmt6.b1', pfam = '../{species}.cce.id'
     output: '{species}.blast.id', '{species}.blast.pfam.id'
     shell: " cut -f2 {input.blast} > {output[0]} && fgrep -f {output[0]} {input.pfam} | sort -u > {output[1]} "

rule extract_sequence:
     input: rules.makedb.input, rules.common_with_pfam.output[1]
     output: '{species}.CCEs.pep.fasta'
     shell: " perl /home/iasonas/bin/jason/extract_fasta.pl {input[0]} {input[1]} > {output} "

