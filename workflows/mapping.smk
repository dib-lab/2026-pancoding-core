rule do_mapping:
    input:
        expand('outputs.mapping/cds/{species}.cds.fa.gz', species=CORE_NAMES),
        expand('outputs.mapping/cds-dedup-ident/{species}.cds.fa.gz', species=CORE_NAMES),
        expand('outputs.mapping/minsig/{s}.cds3.min50.sig.zip', s=CORE_NAMES),
        expand('outputs.mapping/cds-min50/{s}.cds3.min50.fa', s=CORE_NAMES),
        'outputs.mapping/cds-min50-singleclust/all-dedup-95.fa',
        expand('outputs.mapping/bams.cds3.min50.rand/{m}.x.all-dedup-95.depth.txt', m=RAND_METAG),

rule concat_cds:
    input:
        mag_lin_csv=MAG_LIN_CSV,
    params:
        mag_prodigal='outputs.cds3/prodigal/',
        gtdb_prodigal='~/scratch3/2026-gtdb-dl-prodigal',
        gtdb_lincsv=GTDB_LIN_CSV,
        outdir='outputs.mapping/cds',
    output:
        expand('outputs.mapping/cds/{species}.cds.fa.gz', species=CORE_NAMES),
    shell: """
        scripts/extract-species-cds-from-prodigal.py \
            --prodigal {params.mag_prodigal} {params.gtdb_prodigal} \
           --tax {input.mag_lin_csv} {params.gtdb_lincsv} \
           --species inputs.branchwater/names.list -o {params.outdir}

        gzip -9 outputs.mapping/cds/*.cds.fa
    """

rule dedup_ident_cds:
    input:
        'outputs.mapping/cds/{species}.cds.fa.gz'
    output:
        'outputs.mapping/cds-dedup-ident/{species}.cds.fa.gz',
    threads: 16
    shell: """
        cd-hit -c 1.0 -i {input:q} -o {output:q} -T {threads} -M 5000
    """

rule rand_metags_mf_csv:
    input:
        expand(GRIST_RAND100 + "sigs/{m}.trim.sig.zip", m=RAND_METAG)
    output:
        "outputs.mapping/rand-metags.mf.csv",
    shell: """
        sourmash sig collect --abspath -F csv -o {output} {input}
    """

rule screen_min50:
    input:
        sig='inputs.branchwater/queries/{s}.cds3.sig.zip',
        metags='outputs.mapping/rand-metags.mf.csv',
    output:
        'outputs.mapping/minsig/{s}.cds3.min50.sig.zip',
    shell: '''
        scripts/screen-sigs-x-metags.py {input.sig:q} {input.metags:q} \
            -m 50 -o {output:q} -k 21
    '''

# get all the sequences correspondingg to the k-mers present in
# 50% of the rand collection of metagenomes: cds3.min50.fa.
rule kmers_wc:
    input:
        sig='outputs.mapping/minsig/{s}.cds3.min50.sig.zip',
        seqs='outputs.mapping/cds-dedup-ident/{s}.cds.fa.gz',
    output:
        fa=touch('outputs.mapping/cds-min50/{s}.cds3.min50.fa'),
    shell: """
        sourmash sig kmers --sig {input.sig:q} --sequences {input.seqs:q} \
            --save-sequences {output.fa:q} -k 21 || true
    """

# combine all of the cds-min50 into a single file
rule combine_minsig:
    input:
        expand('outputs.mapping/cds-min50/{s}.cds3.min50.fa',
               s=CORE_NAMES),
    output:
        'outputs.mapping/cds-min50/all.fa',
    shell: """
        cat {input:q} > {output:q}
    """

rule cluster_all_minsig:
    input:
        'outputs.mapping/cds-min50/all.fa',
    output:
        fa='outputs.mapping/cds-min50-singleclust/all-dedup-95.fa',
        clstr='outputs.mapping/cds-min50-singleclust/all-dedup-95.fa.clstr',
    threads: 16
    shell: """
        cd-hit -c 0.95 -i {input:q} -o {output.fa:q} -T {threads} -M 5000
    """

rule map_index_rand_cds3_min50:
    input:
        fa='outputs.mapping/cds-min50-singleclust/all-dedup-95.fa',
        metag=GRIST_RAND100 + "trim/{m}.trim.fq.gz",
    output:
        bam='outputs.mapping/bams.cds3.min50.rand/{m}.x.all-dedup-95.bam',
        bai='outputs.mapping/bams.cds3.min50.rand/{m}.x.all-dedup-95.bam.bai',
    threads: 8
    shell: """
        minimap2 -ax sr -t {threads} {input.fa:q} {input.metag:q} | samtools view -b -F 4 - | samtools sort - > {output.bam:q}
        samtools index {output.bam:q}
    """

# calculate coverage txt files for all min50 dedup mappings.
# this gives us direct read-to-gene breadth & depth info, to
# be processed in a notebook.
rule map_rand_cds3_min50_depth:
    input:
        bam='outputs.mapping/bams.cds3.min50.rand/{m}.x.all-dedup-95.bam',
        bai='outputs.mapping/bams.cds3.min50.rand/{m}.x.all-dedup-95.bam.bai',
        fa='outputs.mapping/cds-min50-singleclust/all-dedup-95.fa',
    output:
        'outputs.mapping/bams.cds3.min50.rand/{m}.x.all-dedup-95.depth.txt',
    threads: 8
    shell: """
        samtools depth -aa {input.bam:q} {input.fa:q} > {output:q}
    """
