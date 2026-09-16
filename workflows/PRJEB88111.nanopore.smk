nanopore_runs_df = pl.read_csv('inputs.PRJEB88111.nanopore/filereport_read_run_PRJEB88111-md5.tsv', separator='\t')

NANOPORE_ACCS = nanopore_runs_df['run_accession'].to_list()

wildcard_constraints:
    acc="[^.]+"

rule nanopore:
    input:
        expand('outputs.PRJEB88111.nanopore/datasets/{acc}.fastq.gz',
               acc=NANOPORE_ACCS),
        expand('outputs.PRJEB88111.nanopore/sigs/{acc}.sig.zip',
               acc=NANOPORE_ACCS),
        'outputs.PRJEB88111.nanopore/PRJEB88111.mf.csv',
        expand('outputs.PRJEB88111.nanopore/bams/{acc}.x.all-core.depth.csv',
               acc=NANOPORE_ACCS),
        'outputs.PRJEB88111.nanopore/contigs_to_species.csv',

rule dl:
    params:
        url=lambda w: nanopore_runs_df.filter(pl.col('run_accession') == w.acc)['fastq_ftp'].head(1).item()
    output:
        'outputs.PRJEB88111.nanopore/datasets/{acc}.fastq.gz'
    shell:
        'curl -L {params.url} -o {output}'

rule sketch:
    input:
        'outputs.PRJEB88111.nanopore/datasets/{acc}.fastq.gz',
    output:
        'outputs.PRJEB88111.nanopore/sigs/{acc}.sig.zip',
    shell: '''
        sourmash sketch dna {input} -p k=21,k=31,k=51,dna,scaled=1000,abund \
           --name {wildcards.acc} -o {output}
    '''

rule mf:
    input:
        expand('outputs.PRJEB88111.nanopore/sigs/{acc}.sig.zip',
               acc=NANOPORE_ACCS),
    output:
        'outputs.PRJEB88111.nanopore/PRJEB88111.mf.csv',
    shell: """
        sourmash sig collect -F csv {input} -o {output}
    """

rule manysearch_nanopore:
    input:
        db="outputs.cds3/mag+gtdb.cds.singleton.sig.zip",
        mf="outputs.PRJEB88111.nanopore/PRJEB88111.mf.csv",
    output:
        "outputs.PRJEB88111.nanopore/PRJEB88111.x.mags+gtdb.manysearch.csv",
    threads: 24
    shell: """
        sourmash scripts manysearch -k 21 --scaled=1000 --threshold=0 \
           -c {threads} \
           {input.db} {input.mf} -o {output}
    """

rule concat_refs:
    input:
        "inputs.PRJEB88111.nanopore/GCA_000434935.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCA_002299675.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCA_002320005.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCA_004552595.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCA_004557565.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCA_022774325.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCA_022785155.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCA_034089285.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCA_900546925.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCA_902472425.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCF_002251295.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCF_002706375.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCF_009695905.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCF_009696175.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCF_019052365.1_genomic.fna.gz",
        "inputs.PRJEB88111.nanopore/GCF_021531895.1_genomic.fna.gz",
    output:
        ref="outputs.PRJEB88111.nanopore/all-core.fa",
    shell: """
        gunzip -c {input} > {output}
    """

rule make_contigs_to_species_csv:
    input:
        genomes=[
            "inputs.PRJEB88111.nanopore/GCA_000434935.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCA_002299675.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCA_002320005.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCA_004552595.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCA_004557565.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCA_022774325.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCA_022785155.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCA_034089285.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCA_900546925.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCA_902472425.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCF_002251295.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCF_002706375.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCF_009695905.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCF_009696175.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCF_019052365.1_genomic.fna.gz",
            "inputs.PRJEB88111.nanopore/GCF_021531895.1_genomic.fna.gz",
            ],
        gtdb_reps_csv='inputs.PRJEB88111.nanopore/reps.csv',
    output:
        'outputs.PRJEB88111.nanopore/contigs_to_species.csv'
    shell: """
        scripts/make-contigs-to-species-csv.py {input.gtdb_reps_csv} \
           --genome-dir inputs.PRJEB88111.nanopore -o {output}
    """

rule map:
    input:
        fastq="outputs.PRJEB88111.nanopore/datasets/{acc}.fastq.gz",
        ref="outputs.PRJEB88111.nanopore/all-core.fa",
    output:
        bam="outputs.PRJEB88111.nanopore/bams/{acc}.x.all-core.bam",
        bai="outputs.PRJEB88111.nanopore/bams/{acc}.x.all-core.bam.bai",
    threads: 8
    shell: """
        minimap2 -ax map-ont -t {threads} {input.ref} {input.fastq} | \
            samtools view -b -F 4 - | samtools sort - > {output.bam}
        samtools index {output.bam}
    """

rule depth:
    input:
        bam="outputs.PRJEB88111.nanopore/bams/{acc}.x.all-core.bam",
        fa='outputs.PRJEB88111.nanopore/all-core.fa',
    output:
        "outputs.PRJEB88111.nanopore/bams/{acc}.x.all-core.depth.txt",
    shell: """
        samtools depth -aa {input.bam} {input.fa} > {output}
    """

rule depth_csv:
    input:
        "outputs.PRJEB88111.nanopore/bams/{acc}.x.all-core.depth.txt"
    output:
        "outputs.PRJEB88111.nanopore/bams/{acc}.x.all-core.depth.csv",
    shell: """
        scripts/summarize-mapping-depth-to-csv.py {input} -o {output}
    """
