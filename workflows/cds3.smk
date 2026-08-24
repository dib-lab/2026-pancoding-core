import polars as pl

mag_fasta_df = (pl.read_csv(MAG_NAMES_FASTA)
    .with_columns(ident=pl.col('name').str.split(' ').list.get(0))
)
mag_lin_df = pl.read_csv(MAG_LIN_FASTA)

# AtH_MAG_N_ etc.
MAG_NAMES = (mag_fasta_df
    ['ident'].to_list()
    )


rule make_cds3:
    input:
        expand('outputs.cds3/prodigal/{n}.ffn', n=MAG_NAMES)

# retrieve the fasta for wildcards.g, where g is a MAG ident
def _get_mag_fasta_path(w):
    name = w.g
    xx_df = mag_fasta_df.filter(pl.col('ident') == name)
    assert len(xx_df) == 1
    filename = xx_df['genome_filename'].item()
    return filename
    
rule run_prodigal:
    input:
        g=_get_mag_fasta_path
    output:
        ffn="outputs.cds3/prodigal/{g}.ffn",
        faa="outputs.cds3/prodigal/{g}.faa",
    threads: 1
    shell: """
        prodigal -i {input.g} -d {output.ffn} -a {output.faa} -o /dev/null
    """
