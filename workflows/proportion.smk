rule do_proportion:
    input:
        'outputs.proportion/combined-pangenomes.x.3216.csv'

rule extract_species_sketchx1:
    input:
        dbs='gtdb-foo/gtdb-rs226-k21.species.sig.zip',
    output:
        'outputs.proportion/queries/{species}.pangenome-entire.sig.zip',
    params:
        exact=lambda w: f"{w.species}"
    shell: """
        sourmash sig grep {params.exact:q} {input:q} -o {output:q}
    """

rule extract_species_sketchx2:
    input:
        db='gtdb-foo/gtdb-rs226-k21.species.singlehash.sig.zip',
    output:
        'outputs.proportion/queries/{species}.pangenome-specific.sig.zip',
    params:
        exact=lambda w: f"{w.species}"
    shell: """
        sourmash sig grep {params.exact:q} {input:q} -o {output:q}
    """

rule make_query_pangenome_entire:
    input:
        expand('outputs.proportion/queries/{species}.pangenome-entire.sig.zip',
               species=CORE_NAMES)
    output:
        'outputs.proportion/combined-pangenome-entire.sig',
    shell: """
        sourmash sig merge {input:q} -o {output:q} --set-name='entire'
    """

rule make_query_pangenome_specific:
    input:
        expand('outputs.proportion/queries/{species}.pangenome-specific.sig.zip',
               species=CORE_NAMES)
    output:
        'outputs.proportion/combined-pangenome-specific.sig',
    shell: """
        sourmash sig merge {input:q} -o {output:q} --set-name='specific'
    """

rule make_query_combined:
    input:
        'outputs.proportion/combined-pangenome-entire.sig',
        'outputs.proportion/combined-pangenome-specific.sig',
        expand('outputs.proportion/queries/{species}.pangenome-entire.sig.zip',
               species=CORE_NAMES),
        expand('outputs.proportion/queries/{species}.pangenome-specific.sig.zip',
               species=CORE_NAMES),
    output:
        'outputs.proportion/pangenome-many-queries.sig'
    shell: """
        sourmash sig cat {input:q} -o {output:q}
    """

rule proportion_3216:
    input:
        queries='outputs.proportion/pangenome-many-queries.sig',
        against='3216.manifest.csv',
    output:
        'outputs.proportion/combined-pangenomes.x.3216.csv'
    threads: 16
    shell: """
        sourmash scripts manysearch {input.queries} {input.against} -o {output} \
            -k 21 -t 0 -c {threads}
    """
