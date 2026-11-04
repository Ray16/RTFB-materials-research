"""Apply a prospectively tested resource bound to pending matched native tasks.

No numerical setting, namespace, geometry or running reservation is changed.
The reference run must pass independent native and high-water resource checks.
Every new task keeps both its original inner guard and a stricter outer guard.
"""
import json
from pathlib import Path
import tempfile
import shutil

from redox.core.hashing import sha256_file
from redox.workflows import effective_native_namespace as effective
from redox.validation.namespace_sources import read_checked_review


class MatchedNativeProfiles:
    def __init__(self, resource_review):
        self.reference_path=Path(resource_review).resolve()
        resource=json.loads(self.reference_path.read_text())
        if (resource['status']!='effective_native_resource_profile_checked'
                or resource['tested_run_profile_qualified'] is not True
                or resource['experimental_accuracy_inferred'] is not False):
            raise ValueError('An independently checked prospective resource run is required')
        self.namespace=Path(resource['normal_native_binding']['namespace'])
        self.reference,self.binding=read_checked_review(self.namespace,'native')
        self.profile_path=Path(self.reference['effective_runtime_validation']['profile'])
        self.profile,self.spec=effective.verify(self.profile_path)
        if (self.binding!=resource['normal_native_binding']
                or self.profile['native_namespace']!=str(self.namespace)
                or resource['profile_sha256']!=sha256_file(self.profile_path)):
            raise ValueError('Resource reference does not bind its exact native run')
        # Recompute the resource result rather than trusting its passed label.
        with tempfile.TemporaryDirectory(prefix='redox_resource_reference_') as scratch:
            reconstructed=effective.review(self.profile_path,Path(scratch)/'review.json')
        if reconstructed!=resource:
            raise ValueError('Prospective resource qualification does not reconstruct')
        self.bound=resource['effective_bound_GiB']
        self.reference_plan=json.loads((Path(self.spec['attempt'])/'plan.json').read_text())
        files=[self.reference_path,self.namespace,self.profile_path,Path(self.binding['review']),
               Path(self.spec['attempt'])/'plan.json',Path(__file__).resolve()]
        self.source_sha256={str(p):sha256_file(p) for p in files}
        self.profiles={}

    def verify_sources(self):
        if any(sha256_file(p)!=h for p,h in self.source_sha256.items()):
            raise ValueError('Qualified resource reference or prepared profile changed')

    def for_pending(self, namespace):
        """Return a checked profile only for the same constitution/state/basis."""
        namespace=Path(namespace).resolve()
        self.verify_sources()
        spec=effective.verify_spec(namespace)
        target=json.loads((Path(spec['attempt'])/'plan.json').read_text())
        try:
            effective.validate_class(target,self.reference_plan,self.reference,self.bound)
        except ValueError:
            return None
        output=Path(spec['attempt']).parent.parent/'matched_effective_profile'
        path=output/'profile.json'
        if not output.exists():
            effective.prepare(namespace,Path(self.binding['review']),self.binding['review_sha256'],
                              self.bound,'derived',output,reference_namespace=self.namespace)
        profile,prepared=effective.verify(path)
        if (prepared!=spec or profile['native_namespace']!=str(namespace)
                or profile['reference_review_sha256']!=self.binding['review_sha256']
                or profile['effective_bound_GiB']!=self.bound):
            raise ValueError('Existing prospective profile has a different native task/reference')
        self.source_sha256[str(path)]=sha256_file(path)
        self.profiles[str(namespace)]=dict(profile=str(path),profile_sha256=sha256_file(path),
            effective_bound_GiB=self.bound,qualified_reference=str(self.reference_path),
            qualified_reference_sha256=self.source_sha256[str(self.reference_path)])
        return path


def admission_spec(spec, profile):
    """Scheduling copy only; the original immutable native plan stays at 80 GiB."""
    if profile is None:
        return spec
    data,checked=effective.verify(Path(profile))
    if checked!=spec or spec['maximum_host_RSS_GiB']!=80:
        raise ValueError('Scheduling profile belongs to a different immutable namespace')
    return dict(spec,maximum_host_RSS_GiB=data['effective_bound_GiB'])


def first_fitting(pending, resources, profiles, fits):
    """Avoid blocking smaller matched work behind an unrelated larger task."""
    for index,(path,spec) in enumerate(pending):
        profile=profiles.get(str(path))
        if fits(admission_spec(spec,profile),resources):
            return index,profile
    return None,None


def output_room(path, minimum_free_GiB=64):
    """Admission floor on the shared output volume; never interrupt live work."""
    free=shutil.disk_usage(path).free
    return dict(path=str(Path(path).resolve()),free_bytes=free,
                minimum_free_GiB=minimum_free_GiB,admit_new_work=free>=minimum_free_GiB*2**30)
