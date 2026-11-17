# Portable occupied-sampling replay

From this repository or an extracted source archive, run:

```sh
python examples/occupied_sampling_replay/replay.py
```

The reader uses the shipped exact bytes for all 183 original input bindings,
including twenty XYZ geometries and the historical producer sources. It repeats
graph, chirality, proper-rotation RMSD, duplicate grouping and finite-partition
checks under the unchanged physical contexts. The frozen 240-slot roster,
31 qualified native records and eleven incomplete cells remain explicit.

Expanding the complete occupied cell from ten to twenty starts changes its
finite-set G by -2.399 kJ/mol at the primary cutoff and fails the declared
2-kJ/mol tolerance. No original calculation directory, GPU or network is
required. This example uses the repository's shared readers; copying only this
example folder without the source package is not supported. Required numerical
package versions are listed in `requirements.txt`.

Historical paths remain provenance strings. Every mapped file must match its
recorded hash; missing mappings never fall back to arbitrary original paths.
The currently executing producer and its historical source are kept distinct.
This is a cached sampling replay, not new SCF/Hessian/spin validation, an occupied
equilibrium cycle or evidence of experimental populations.

The maintained exporter is `scripts/release/prepare_sampling_inputs.py`. It can
create a new packet from frozen publication artifacts or repackage an existing
packet using `--existing PATH/portable_inputs.json`. Existing packets are never
overwritten.
