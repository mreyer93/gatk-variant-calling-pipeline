################################################################################
# Create the index neccessary for BWA alignment.
rule build_ref_index:
    input: REF_FILE
    output: 
        expand(REF_FILE + ".{bwa_ext}",  bwa_ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
    shell: """
        bwa index {input}
    """

################################################################################
#Create the fasta index for alignment
rule build_ref_faidx:
    input: REF_FILE
    output: REF_FILE + ".fai"
    shell: """
        samtools faidx {input}
    """

################################################################################
# Create sequence dictionary for GATK variant calling.
rule create_ref_dict:
    input:  REF_FILE
    output: REF_FILE.split(".")[0] + ".dict"
    shell: """
        gatk CreateSequenceDictionary --REFERENCE {input}
    """
