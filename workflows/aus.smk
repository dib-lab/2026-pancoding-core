# explore PRJNA526405 study
rule do_aus:
    input:
        expand('outputs.aus/{s}.isect.abund.csv', s=CORE_NAMES),
        'outputs.aus/combined-pangenomes.x.aus.csv',

rule aus_intersection_explosion:
    input:
        queries=expand('inputs.branchwater/queries/gtdb+mags.{s}.cds3.sig.zip',
                       s=CORE_NAMES),
        metag='PRJNA526405.aus/accs.txt',
    output:
        directory('outputs.aus/aus_isect.d')
    shell: """
        scripts/extract-intersection-sigs.py --query {input.queries:q} \
             --metagenomes {input.metag} --output-dir {output}
    """

rule aus_intersection_collect:
    input:
        'outputs.aus/aus_isect.d'
    output:
        expand('outputs.aus/aus_isect.collect.d/{s}.collected.sig.zip',
               s=CORE_NAMES)
    shell: """
        scripts/collect-intersection-sigs-by-species.py {input} \
             -o outputs.aus/aus_isect.collect.d/
    """

rule aus_intersections_kmer_stats:
    input:
        expand('outputs.aus/aus_isect.collect.d/{s}.collected.sig.zip',
               s=CORE_NAMES)
    output:
        'outputs.aus/aus_kmer_stats.csv'
    shell: """
        scripts/isect-extract-kmer-stats.py {input:q} -o {output:q}
    """

rule aus_intersection_do_fits_for_abund:
    input:
        'outputs.aus/aus_isect.collect.d/{s}.collected.sig.zip'
    output:
        abund='outputs.aus/{s}.isect.abund.csv',
        kmer='outputs.aus/{s}.kmer_stats.csv'
    shell: """
        #scripts/fit-isect-abund-experimental.py {input:q} -o {output:q}
        scripts/isect-estimate-species-abundances.py {input:q} \
            --output-kmer-stats-csv {output.kmer:q} \
            --output-abund-estimates {output.abund:q}
    """

rule proportion_aus:
    input:
        queries='outputs.proportion/pangenome-many-queries.sig',
        against='PRJNA526405.aus/accs.txt',
    output:
        'outputs.aus/combined-pangenomes.x.aus.csv'
    threads: 16
    shell: """
        sourmash scripts manysearch {input.queries} {input.against} -o {output} \
            -k 21 -t 0 -c {threads}
    """
