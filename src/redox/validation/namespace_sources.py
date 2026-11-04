"""Verify namespace-produced records through their two declared path aliases."""
import json
import math
from pathlib import Path, PurePosixPath
import time

from redox.core.hashing import sha256_file
from redox.validation.archived_sources import ArchivedSources


class NamespaceSources(ArchivedSources):
    def __init__(self, spec):
        super().__init__(spec['prior_root'], spec['archive'], relative_to_prior=True)
        self.aliases = [(PurePosixPath(spec['original_campaign']), Path(spec['attempt']).resolve()),
                        (PurePosixPath(spec['original_raw_parent']), Path(spec['raw_parent']).resolve())]
        prior = Path(spec['prior_root'])
        if (not Path(spec['original_campaign']).is_relative_to(prior/'results/polymer_models/campaigns')
                or not Path(spec['original_raw_parent']).is_relative_to(prior/'calcs/polymer_models')):
            raise ValueError('Only the declared polymer attempt and raw-parent aliases are supported')
        for original, target in self.aliases:
            if '..' in original.parts or not original.is_absolute() or not target.is_dir():
                raise ValueError('Explicit existing namespace directories required')
            if target.is_relative_to(self.archive) or target.is_relative_to(Path(str(self.prior))):
                raise ValueError('New attempt aliases cannot write over archived or current-checkout paths')

    def path(self, original):
        path = PurePosixPath(str(original))
        if '..' in path.parts:
            raise ValueError('Evidence paths cannot traverse parents')
        if not path.is_absolute():
            path = self.prior/path
        for source, target in self.aliases:
            try:
                relative = path.relative_to(source)
            except ValueError:
                continue
            return target.joinpath(*relative.parts)
        return super().path(path)


def read_checked_review(spec_path, kind, *, visibility_timeout=30):
    """Wait only for file visibility; changed scientific bytes fail immediately."""
    spec_path = Path(spec_path)
    before = sha256_file(spec_path)
    spec = json.loads(spec_path.read_text())
    if kind not in ('native', 'energy'):
        raise ValueError('Native or electronic review required')
    target = Path(spec['attempt'])/('independent_'+kind+'_review.json')
    deadline = time.monotonic()+visibility_timeout
    while True:
        try:
            report = json.loads(target.read_text())
            break
        except FileNotFoundError:
            if time.monotonic() >= deadline:
                raise
            time.sleep(.5)
    expected = ('alternate_native_state_independently_reconstructed' if kind == 'native'
                else 'alternate_electronic_geometry_and_resource_checked_native_G_pending')
    method_key = 'experimental_method_id' if kind == 'native' else 'portable_method_id'
    if (report['status'] != expected or report['campaign'] != spec['original_campaign']
            or report['task']['key'] != spec['key']
            or report[method_key] != spec['method_id'] or report['primary_candidate_admission'] is not False
            or report['experimental_model']['solvent'] != 'acetonitrile'):
        raise ValueError('Independent review has a different graph, method, solvent or scope')
    if kind == 'energy' and report['native_Hessian_or_G_qualified'] is not False:
        raise ValueError('Electronic review cannot supply native thermochemistry')
    reader = NamespaceSources(spec)
    reader.verify_all(report['source_sha256'])
    plan_path = str(Path(spec['original_campaign'])/'plan.json')
    plan = json.loads(reader.verify(plan_path, report['source_sha256'][plan_path]).read_text())
    if plan['tasks'] != [report['task']] or plan['experimental_model'] != report['experimental_model']:
        raise ValueError('Independent task/model differs from the exact frozen plan')
    raw = Path(spec['original_raw_root'])
    if kind == 'energy':
        raw_path = str(raw/'experimental_energy.json')
        stored = json.loads(reader.verify(raw_path, report['source_sha256'][raw_path]).read_text())
        for name, field in (('E_smd_eV', 'e_smd_eV'), ('E_gas_eV', 'e_gas_eV')):
            if not math.isfinite(report[name]) or abs(report[name]-stored[field]) > 1e-9:
                raise ValueError('Reported electronic energy differs from the bound raw record')
    else:
        raw_path = str(raw/'experimental_partition.json')
        stored = json.loads(reader.verify(raw_path, report['source_sha256'][raw_path]).read_text())
        if set(report['G_by_cutoff_eV']) != {'25', '50', '100'}:
            raise ValueError('Three rotor-cutoff native free energies required')
        for rotor, value in report['G_by_cutoff_eV'].items():
            if not math.isfinite(value) or abs(value-stored['G_by_cutoff_eV'][rotor]) > 1e-8:
                raise ValueError('Reported native free energy differs from its bound raw partition')
    reader.recheck()
    if sha256_file(spec_path) != before:
        raise ValueError('Namespace changed during collection')
    return report, dict(review=str(target), review_sha256=sha256_file(target),
                        namespace=str(spec_path), namespace_sha256=before,
                        verified_bindings=len(reader.receipts), original_source_hashes_preserved=True)
