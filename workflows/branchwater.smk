rule do_branchwater:
    input:
        expand("outputs.branchwater/parquet/gtdb+mags.manysearch.{species}.parquet",
               species=CORE_NAMES),
        expand("outputs.branchwater/parquet/gtdb.manysearch.{species}.parquet",
               species=CORE_NAMES)

rule extract_cds3_sketches:
    input:
        expand("inputs.branchwater/queries/gtdb+mags.{species}.cds3.sig.zip",
               species=CORE_NAMES),
        expand("inputs.branchwater/queries/gtdb.{species}.cds3.sig.zip",
               species=CORE_NAMES)

rule extract_one_sketch:
    input:
        'outputs.cds3/mag+gtdb.cds.singleton.sig.zip',
    output:
        'inputs.branchwater/queries/gtdb+mags.{species}.cds3.sig.zip',
    params:
        exact=lambda w: f"{w.species} " # add space to require exact matching
    shell: """
        sourmash sig grep {params.exact:q} {input:q} -o {output:q}
    """

rule extract_one_gtdb:
    input:
        GTDB_SINGLETON_SIG_ZIP,
    output:
        'inputs.branchwater/queries/gtdb.{species}.cds3.sig.zip',
    params:
        exact=lambda w: f"{w.species} " # add space to require exact matching
    shell: """
        sourmash sig grep {params.exact:q} {input:q} -o {output:q}
    """

rule search_species:
    input:
        q='inputs.branchwater/queries/{label}.{species}.cds3.sig.zip',
        db=BW_ROCKSDB,
    output:
        "outputs.branchwater/csv/{label}.manysearch.{species}.csv",
    threads: 1
    shell: """
        sourmash scripts manysearch -c 1 -t 0 -k 21 -s 1000 \
            {input.q:q} {input.db:q} -o {output:q}
    """

rule add_pvals_branch_mag_gtdb:
    input:
        csv="outputs.branchwater/csv/gtdb+mags.manysearch.{species}.csv",
        queries='outputs.cds3/mag+gtdb.cds.singleton.sig.zip',
        against='wort-703k.mf.csv',
    output:
        pq='outputs.branchwater/parquet/gtdb+mags.manysearch.{species}.parquet',
    shell: """
        scripts/add-pvals-to-manysearch.py {input.csv:q} \
            --queries {input.queries:q} --against {input.against:q} \
            -P {output.pq:q} --remove-missing
    """

rule add_pvals_branch_gtdb:
    input:
        csv="outputs.branchwater/csv/gtdb.manysearch.{species}.csv",
        queries=GTDB_SINGLETON_SIG_ZIP,
        against='wort-703k.mf.csv',
    output:
        pq='outputs.branchwater/parquet/gtdb.manysearch.{species}.parquet',
    shell: """
        scripts/add-pvals-to-manysearch.py {input.csv:q} \
            --queries {input.queries:q} --against {input.against:q} \
            -P {output.pq:q} --remove-missing
    """
