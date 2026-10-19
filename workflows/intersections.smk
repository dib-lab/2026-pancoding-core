rule do_intersections:
    input:
        expand('outputs.isect/{s}.isect.abund.csv', s=CORE_NAMES),

rule intersection_explosion:
    input:
        queries=expand('inputs.branchwater/queries/gtdb+mags.{s}.sig.zip',
                       s=CORE_NAMES),
        metag='3216-sigs.txt',
    output:
        directory('outputs.isect/3216_isect.d')
    shell: """
        scripts/extract-intersection-sigs.py --query {input.queries:q} \
             --metagenomes {input.metag} --output-dir {output}
    """

rule intersection_collect:
    input:
        'outputs.isect/3216_isect.d'
    output:
        expand('outputs.isect/3216_isect.collect.d/{s}.collected.sig.zip',
               s=CORE_NAMES)
    shell: """
        scripts/collect-intersection-sigs-by-species.py {input} \
             -o outputs.isect/3216_isect.collect.d/
    """

rule intersections_kmer_stats:
    input:
        expand('outputs.isect/3216_isect.collect.d/{s}.collected.sig.zip',
               s=CORE_NAMES)
    output:
        'outputs.isect/3216_kmer_stats.csv'
    shell: """
        scripts/isect-extract-kmer-stats.py {input:q} -o {output:q}
    """

rule intersection_do_fits_for_abund:
    input:
        'outputs.isect/3216_isect.collect.d/{s}.collected.sig.zip'
    output:
        abund='outputs.isect/{s}.isect.abund.csv',
        kmer='outputs.isect/{s}.kmer_stats.csv'
    shell: """
        #scripts/fit-isect-abund-experimental.py {input:q} -o {output:q}
        scripts/isect-estimate-species-abundances.py {input:q} \
            --output-kmer-stats-csv {output.kmer:q} \
            --output-abund-estimates {output.abund:q}
    """
