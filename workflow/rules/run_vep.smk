#rule download_cache:
#    output:
#        cache=CACHE
#    container:
#        VEP_CONTAINER
#    shell:
#        """
#        singularity exec \
#            --bind {output.cache}:/opt/vep/.vep \
#            ensembl-vep.sif \
#            vep_install -a cf -s homo_sapens_refseq -y GRCh38 \
#            --CACHEDIR /opt/vep/.vep --CONVERT --NO_HTSLIB
#        """

rule prepare_hgvs_inpug:
    input:
        var_excel=VAR_EXCEL,
    output:
        hgvs=HGVS
    params:
        script="workflow/scripts/prepare_hgvs.py"
    shell:
        """
        python {params.script} {input.var_excel} {output.hgvs}
        """


rule download_vep_cache:
    output:
        cache_dir=directory(CACHE_PATH),
        cache=CACHE_FILE
    params:
        url=CACHE_URL,
        file_name=CACHE_URL.split('/')[-1],
        outdir=os.path.dirname(CACHE_PATH)
    shell:
        """
        cd resources/
        curl -O {params.url}
        echo {params.file_name}
        tar xzf {params.file_name}
        """


rule run_vep:
    input:
        hgvs=HGVS,
        cache=CACHE_PATH
    output:
        tab=VEP_OUT
    container:
        VEP_CONTAINER
    shell:
        """
        vep --input_file {input.hgvs} --format hgvs \
        --cache --dir_cache resources/ --refseq \
        --assembly GRCh38 --tab --force_overwrite \
        --output_file {output.tab}
        """

