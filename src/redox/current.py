"""Current public CLI with Nitro-preferred UMA and source-frozen replay modules."""
from redox import __main__ as original


def main(argv=None):
    previous=original.COMMANDS
    original.COMMANDS=dict(previous,
        uma=('redox.qm.uma_current',(),'gas-phase warm starts; prefer Nitro with explicit V100 fallback'),
        dft=('redox.qm.uma_current',('--stage','dft'),'DFT/SMD with current or retained warm-start provenance'),
        launch=('redox.workflows.launch_current',(),'plan current DFT/UMA shards; --execute submits through the reservation gate'))
    try:return original.main(argv)
    finally:original.COMMANDS=previous


if __name__=='__main__':raise SystemExit(main())
