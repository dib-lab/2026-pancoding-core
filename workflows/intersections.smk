rule do_intersections:
    input:
        expand('outputs.isect/{s}.isect.abund.csv', s=CORE_NAMES),

rule intersection_explosion:
    input:
        queries=expand('inputs.branchwater/queries/{s}.sig.zip',
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

rule intersection_do_fits_for_abund:
    input:
        'outputs.isect/3216_isect.collect.d/{s}.collected.sig.zip'
    output:
        'outputs.isect/{s}.isect.abund.csv'
    shell: """
        scripts/fit-isect-abund.py {input:q} -o {output:q}
    """

