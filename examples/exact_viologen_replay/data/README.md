# Exact viologen cached CPU example

Run `python replay.py replay --package .` after installing the hash-locked NumPy requirement with CPython 3.11/Linux x86_64. This reconstructs two Fc-referenced potentials and a local quaternization from ten cached states. Nine Hessians/mode sets and nuclear partitions reconstruct; six recent states also retain displaced-gradient arrays. Historical Fc/site FD entries remain metadata-only. No SCF, xTB library or original project path is needed.

The example uses a capped local graft, MeCN, 298.15 K and separated chloride. Individual nuclear terms retain the source 1-atm convention; all displayed balanced reactions have zero net solute count, so the 1-M correction cancels. No ion-pairing, counterion-exchange, reflux-temperature or polymer correction is supplied. Spin-only Fc and declared candidate spin states are retained. This is finite sampling and does not admit a complete material or revalidate full historical executable/scratch closure.
