# Chemical TOC artwork

[graphical_abstract.pdf](graphical_abstract.pdf) is the readable vector preview;
[EPS](graphical_abstract.eps) and [600 dpi RGB TIFF](graphical_abstract.tif) are
submission-format exports. [current.json](current.json) identifies the source,
immutable artwork manifest and [export/visual checks](graphical_abstract_check.json).

The exact PMDI precursor and capped graft-site structures come from the checked
candidate snapshot. Atom connectivity and the recorded charge/spin ladder are
checked before drawing. The lower-right open-chain cartoon represents neighbor
assessment; it omits substituent chemistry and stereochemistry and is not an
XYZ geometry, measured loading, association result or electron-transfer rate.
The arrow maps identities and does not establish a successful attachment route.
The graphic retains the requirement that all ten candidate axes pass.

The artwork is 3.25 by 1.75 inches. Color, grayscale, reduced-size and independently
rendered EPS previews have been inspected. The design follows the
[ACS TOC specifications](https://pubsapp.acs.org/paragonplus/submission/toc_abstract_graphics_guidelines.pdf).
The graphic is now inserted on the checked manuscript's final page and included
in its self-contained upload. Author review remains pending.

Regenerate with all six native thread caps set to one:

```bash
PYTHONPATH=src .venv-figures/bin/python scripts/plotting/graphical_abstract.py \
  --out NEW_OUTPUT_STEM
```

The producer refuses to overwrite an existing artwork set. The source snapshot
and every export hash are recorded in the artwork manifest. Preview renders
remain under `calcs/figures/graphical_abstract_2026-10-23`.

The standalone export review predates PDF insertion and is preserved unchanged.
The separate [integration check](graphical_abstract_integration_check.json) binds
the completed insertion to the current PDF and visual review. This avoids a
circular dependency between source evidence and the PDF derived from it.
