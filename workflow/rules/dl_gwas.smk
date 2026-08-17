


rule download_gwas:
    output:
        GWAS_OUT

    params:
        script=FETCH_SCRIPT
    shell:
        """
        python {params.script}
        """

