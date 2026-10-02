#####################################################################################################
### INTERVAL SCATTER ################################################################################
#####################################################################################################
# Mutect2 and HaplotypeCaller run once per shard of the target region, in parallel, and
# their outputs are gathered afterwards. On Apple Silicon this is where the speed comes
# from: Intel's GKL library does not load on arm64, so GATK falls back to its single-
# threaded Java PairHMM and --native-pair-hmm-threads does nothing. On x86 with AVX the
# native PairHMM does use threads, so a couple per shard still help there.
#
# Shards are contiguous runs of whole intervals, balanced by length. Overlapping and
# abutting intervals are merged first, as GATK merges them, so no interval is split
# between shards and each caller sees the same regions it would see unscattered.
#
# Expects the including setup file to define REF_FILE and scatter_targets_bed (a BED file,
# or '' to scatter over every contig in the reference).
import platform

# Fixed rather than derived from --cores: a different shard count changes the jobs, so
# Snakemake would rerun every caller whenever --cores changed.
scatter_count = int(config.get('scatter_count', 8))
gatk_threads = int(config.get('gatk_threads',
                              2 if platform.machine() in ('x86_64', 'AMD64') else 1))
gatk_java_mem_gb = int(config.get('gatk_java_mem_gb', 4))


def _scatter_shards(bed, fai, n):
    """The merged intervals of `bed` (every reference contig when bed is ''), as at most
    n contiguous shards of similar total length. Coordinates are BED (0-based, half-open)."""
    with open(fai) as f:
        contigs = [line.split('\t')[:2] for line in f if line.strip()]
    order = {c: i for i, (c, _) in enumerate(contigs)}
    if bed:
        ivs = []
        with open(bed) as f:
            for line in f:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                c, s, e = line.split('\t')[:3]
                ivs.append((c, int(s), int(e)))
    else:
        ivs = [(c, 0, int(length)) for c, length in contigs]
    ivs.sort(key=lambda iv: (order.get(iv[0], len(order)), iv[0], iv[1], iv[2]))
    merged = []
    for c, s, e in ivs:
        if merged and merged[-1][0] == c and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([c, s, e])
    n = max(1, min(n, len(merged)))
    total = sum(e - s for _, s, e in merged)
    shards, current, done = [], [], 0
    for c, s, e in merged:
        current.append((c, s, e))
        done += e - s
        if len(shards) < n - 1 and done >= total * (len(shards) + 1) / n:
            shards.append(current)
            current = []
    if current:
        shards.append(current)
    return shards


SCATTER = _scatter_shards(scatter_targets_bed, REF_FILE + '.fai', scatter_count)
SHARDS = ['{:04d}'.format(i) for i in range(len(SCATTER))]

wildcard_constraints:
    shard = r'\d{4}'

rule scatter_intervals:
    input:
        [scatter_targets_bed] if scatter_targets_bed else [],
        fai = REF_FILE + '.fai',
    output:
        expand(join(outdir, '00_scatter/{shard}.bed'), shard=SHARDS)
    run:
        for path, shard in zip(output, SCATTER):
            with open(path, 'w') as f:
                for c, s, e in shard:
                    f.write(f'{c}\t{s}\t{e}\n')
