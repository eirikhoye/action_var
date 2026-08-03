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

rule download_vep_cache:
    output:
        cache=CACHE
    params:
        url=CACHE_URL,
        file_name=CACHE_URL.split('/')[-1]
    shell:
        """
        cd resources/
#        curl -O {params.url}
        echo {params.file_name}
        tar xzf {params.file_name}
        """

rule run_vep:
    input:
        hgvs=HGVS,
        cache=CACHE
    output:
        tab=VEP_OUT
    container:
        VEP_CONTAINER
    shell:
        """
        vep --input_file {input.hgvs} --format hgvs \
        --cache --dir_cache {input.cache} --offline --refseq \
        --assembly GRCh38 --tab --force_overwrite \
        --output_file {output.tab}
        """