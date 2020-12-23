configfile: 'config.yaml'

from re import escape
wildcard_constraints:
    index = '|'.join([escape(x) for x in config['sample']])

rule dl:
    output:
        'data/{index}.fastq.gz'
    conda: 
        'envs/1_alignment.yaml'
    shell:
        """
        cd data
        fastq-dump --gzip {wildcards.index} 
        cd ..
        """

rule cut:
    input: 
        rules.dl.output
    output:
        temp('data/{index}.cu.fastq')
    params:
        adaptor = 'TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG'
    conda: 
        'envs/1_alignment.yaml'
    shell:
        """
        cutadapt -a {params.adaptor} -j 6 --discard-untrimmed -m 22 -M 30 -o {output} {input}
        """
        
rule align:
    input:
        rules.cut.output
    output:
        temp('data/{index}.cu.filtered.bam')
    params:
        ref = config['reference'],
        n = lambda wc: 'data/' + wc.index
    conda: 
        'envs/1_alignment.yaml'
    shell:
        """
        bwa aln -t 16 {params.ref} {input} > {params.n}.cu.sai
        bwa samse {params.ref} {params.n}.cu.sai {input} | samtools view -bS | samtools view -bF 4 -q 1 > {output}
        rm -f {params.n}.cu.sai
        """

rule bamindex:
    input:
        rules.align.output
    output:
        'data/{index}.bam'
    conda: 
        'envs/1_alignment.yaml'
    shell:
        """
        picard SortSam I={input} O={output} SORT_ORDER=coordinate; 
        picard BuildBamIndex I={output} 
        """


rule all:
    input:
        expand(rules.bamindex.output, index = config['sample'])
