#!/usr/bin/env python3
"""Guard the CADD licensing boundary across the shipped configs.

CADD is licensed for non-commercial use only (https://cadd.gs.washington.edu/download),
so the *_commercial.yaml configs exist to keep it switched off for paid client work. That
is a promise made in the README and in manual/requirements.md, and it is one flipped
boolean away from being silently untrue - which is exactly the kind of thing nobody
notices until a client's legal review does.

This test re-derives, from each config, whether the pipeline would actually request a
CADD output, using the same gating logic the snakefiles use:

  somatic  - call_bam_GATK.snakefile only appends the CADD targets when
             `skip_annotation` is False.
  germline - call_bam_GATK_germline.snakefile appends the combined-annotation target
             when `skip_annotation` is False, and _cadd_input_germline() in
             scripts/variant_annotation_germline.smk then returns [] when `skip_cadd`
             is True, so the CADD rule is never pulled into the DAG.

Run: python3 test/test_config_licensing.py
Needs only pyyaml, not a full snakemake environment.
"""

import sys
from os.path import join, dirname, abspath

import yaml

CONFIG_DIR = join(dirname(dirname(abspath(__file__))), 'call_bam_GATK')

# config filename -> (pipeline, must CADD be able to run?)
EXPECTED = {
    'config_call_bam_GATK.yaml':                    ('somatic',  True),
    'config_call_bam_GATK_local.yaml':              ('somatic',  False),
    'config_call_bam_GATK_commercial.yaml':         ('somatic',  False),
    'config_call_bam_GATK_germline.yaml':           ('germline', True),
    'config_call_bam_GATK_germline_local.yaml':     ('germline', False),
    'config_call_bam_GATK_germline_commercial.yaml':('germline', False),
}


def cadd_would_run(config, pipeline):
    """Mirror the snakefiles' gating. Defaults match setup.smk / setup_germline.smk."""
    if config.get('skip_annotation', False):
        # No annotation targets are requested at all, so the CADD rule is unreachable.
        return False
    if pipeline == 'germline':
        # Granular flag: the combine rule's CADD input becomes [] and the rule is unreachable.
        return not config.get('skip_cadd', False)
    # The somatic pipeline has no granular skip_cadd; annotation on means CADD on.
    return True


def main():
    failures = []
    print('config'.ljust(48), 'pipeline'.ljust(10), 'CADD runs', ' expected')
    print('-' * 84)
    for name, (pipeline, expected) in sorted(EXPECTED.items()):
        path = join(CONFIG_DIR, name)
        try:
            with open(path) as fh:
                config = yaml.safe_load(fh)
        except FileNotFoundError:
            failures.append('{}: config file is missing'.format(name))
            print(name.ljust(48), pipeline.ljust(10), 'MISSING')
            continue

        actual = cadd_would_run(config, pipeline)
        ok = actual == expected
        print(name.ljust(48), pipeline.ljust(10),
              str(actual).ljust(9), str(expected).ljust(8), '' if ok else '  <-- FAIL')
        if not ok:
            failures.append(
                '{}: CADD would{} run, expected it to{}'.format(
                    name, '' if actual else ' not', '' if expected else ' not'))

        # A commercial germline run drops CADD, so it needs AlphaMissense in its place.
        # Without this the run silently produces no deleteriousness annotation at all.
        if 'commercial' in name and pipeline == 'germline':
            if not config.get('alphamissense_file'):
                failures.append(
                    '{}: skips CADD but sets no alphamissense_file, so the run would '
                    'have no deleteriousness annotation at all'.format(name))

    print()
    if failures:
        print('FAILED ({} problem{})'.format(len(failures), '' if len(failures) == 1 else 's'))
        for f in failures:
            print('  - ' + f)
        return 1
    print('PASSED: every commercial and local config keeps CADD switched off, and both '
          'full configs still run it.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
