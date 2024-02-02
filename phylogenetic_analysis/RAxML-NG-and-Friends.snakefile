import os
import re
import subprocess

include: "snakemake_utils.smk"

###=============================================================== Snakemake Pipeline =========================================================================###


seqs = {f for f in os.listdir(".") if f.endswith((".fasta", ".fa", ".faa"))}

rule all:
    input:
        "all_genes.trimmed.aln.phy.raxml.support.tree",
        "all_genes.trimmed.aln.phy.raxml.support.tree.svg"


rule concatenate:
    input: seqs=seqs
    output: all_fasta="all_genes.fasta"
    shell: """ cat {input.seqs} > {output.all_fasta} """


rule mafft:
    input: all_fasta="all_genes.fasta"
    output: aln="all_genes.aln"
    conda: "envs/phylo.yaml"
    shell: "mafft --auto {input.all_fasta} > {output.aln}"


rule trimal:
    input: aln="all_genes.aln"
    output: trm="all_genes.aln.trimmed"
    params: thr=config['trimal_threshold']
    conda: "envs/phylo.yaml"
    shell: """ trimal -in {input.aln} -out {output.trm} -fasta -gt {params.thr} """


rule convert:
    input: trm="all_genes.aln.trimmed"
    output: phy="all_genes.aln.trimmed.phy"
    message: "Converting trimmed alignment to phy"
    shell: "python3 /home/iasonas/snakemake/RAxML-NG-and-Friends/workflow/scripts/ETElib.py --alignment {input.trm}"


rule pythia:
    input: phy="all_genes.aln.trimmed.phy"
    output: "all_genes.aln.pythia.out"
    conda: "envs/phylo.yaml"
    shell: """ pythia --msa {input.phy} -r raxml-ng --removeDuplicates -o {output} """


rule raxml:
    input: pythia="all_genes.aln.pythia.out"
    output: "all_genes.trimmed.aln.phy.raxml.support"
    threads: config['raxmlng_threads']
    conda: "envs/phylo.yaml"
    params:
        rtn=config['random_tree_number'],
        ptn=config['parsimony_tree_number'],
        outgroup=config['outgroup']

    shell: """ best_model='LG+G4' && outgroup_arg="" && [ -z {params.outgroup} ] || outgroup_arg="--outgroup {params.outgroup}" && raxml-ng --all --msa all_genes.aln.trimmed.phy --model $best_model --tree rand{{params.rtn}}, pars{{params.ptn}} --threads 40 --workers auto $outgroup_arg """
   
rule midpoint_root:
    input: "all_genes.trimmed.aln.phy.raxml.support"
    output: "all_genes.trimmed.aln.phy.raxml.support.tree"
    run:
        if not re.search ("\w", config['outgroup']):
            subprocess.run(["python3", "/home/iasonas/snakemake/RAxML-NG-and-Friends/workflow/scripts/ETElib.py", "--tree", "all_genes.trimmed.aln.phy.raxml.support", "--midpoint"])
        else:
            subprocess.run(["cp", "all_genes.trimmed.aln.phy.raxml.support", "all_genes.trimmed.aln.phy.raxml.support.tree"])


rule visualize_tree:
    input: "all_genes.trimmed.aln.phy.raxml.support.tree"
    output: "all_genes.trimmed.aln.phy.raxml.support.tree.svg"
    shell: """ python3 /home/iasonas/snakemake/RAxML-NG-and-Friends/workflow/scripts/ETElib.py --tree {input} --visualize """
