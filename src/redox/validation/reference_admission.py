"""Current-state gates for experimental residuals and calibration inputs."""
import math

from redox.validation.spin_evidence import state_spin_evidence,spin_comparison_running


def event_admission(gid,event,calculated=None):
    """Return (allowed, reason) without running or promoting any calculation.

    Measured-source quality does not establish computed-state quality. A stale
    serialized potential, unfinished spin decision or unmatched reference stays
    outside calibration, even when an experimental label is tier A.
    """
    from redox.core.common import read_result,free_energy
    from redox.properties.potentials import fc_live,thermal_definition
    if spin_comparison_running(gid):return False,'active spin comparison'
    try:
        ox,red=event.split('->')
        left,right=read_result(gid,ox),read_result(gid,red)
        GO,qO=free_energy(left);GR,qR=free_energy(right);reference=fc_live()
        if (GO is None or GR is None or qO!='ok' or qR!='ok' or reference is None
                or not math.isfinite(reference)):
            return False,'incomplete native endpoints/reference'
        if {thermal_definition(r) for r in (left,right,read_result('ferrocene','ox'))}!={'native_solution_dft_qrrho_spin'}:
            return False,'unmatched native thermal definitions'
        for state in (ox,red):
            evidence=state_spin_evidence(gid,state)
            if not evidence['spin_primary_admissible']:
                return False,evidence['spin_confidence']
        if calculated is not None:
            if isinstance(calculated,bool) or not math.isfinite(float(calculated)):
                return False,'nonfinite serialized potential'
            if abs(float(calculated)-round(GO-GR-reference,3))>1e-6:
                return False,'stale serialized potential'
        return True,'current native endpoints and stated spin model'
    except (ValueError,TypeError,KeyError,OSError):
        return False,'unreadable or changing state evidence'
