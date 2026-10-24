#!/usr/bin/env python
"""Add the existing workflow-reader hash field to a canonical GPU energy marker.

Keep the original marker and name this adapter explicitly. Quantum records,
normal completion receipts, calculation inputs and numerical values are unchanged.
"""
import argparse
import json
from pathlib import Path
import time

from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint


def normalize(campaign):
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    if len(plan['tasks'])!=1:raise ValueError('One exact electronic handoff is required')
    key=plan['tasks'][0]['key'];receipt=campaign/'energy/done'/(key+'.json')
    verify_state_checkpoint(receipt)
    parent=json.loads(receipt.read_text());status=Path(parent['status_path'])
    record=json.loads(status.read_text())
    if record['status']!='state_records_verified' or record.get('physical_QC')!='passed':
        raise ValueError('The normal canonical/physical checks must pass')
    marker=campaign/'workspace_energy_status.json';data=json.loads(marker.read_text())
    if (data['status']!='stationary_electronic_case_checked_native_G_pending'
            or data['source_receipt']!=str(receipt) or data['source_receipt_sha256']!=api.digest(receipt)):
        raise ValueError('The marker points to another or changed canonical receipt')
    expected=api.digest(status)
    if 'energy_receipt_sha256' in data:
        if data['energy_receipt_sha256']!=expected:raise ValueError('Conflicting workflow status hash')
        return
    directory=campaign/'handoff_schema';directory.mkdir(exist_ok=False)
    previous=directory/'original_marker.json';previous.write_bytes(marker.read_bytes())
    data.update(energy_receipt=str(status),energy_receipt_sha256=expected,
        schema_adapter=dict(source=str(Path(__file__).resolve()),source_sha256=api.digest(__file__),
            original_marker=str(previous),original_marker_sha256=api.digest(previous),applied=time.time(),
            quantum_records_changed=False))
    atomic_write_json(marker,data)
    print('Canonical energy marker normalized; quantum records and receipts unchanged')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True)
    a=p.parse_args();normalize(a.campaign.resolve())


if __name__=='__main__':main()
