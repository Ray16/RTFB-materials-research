# Two independently checked missing hybrid-energy states

[Summary and bindings](summary.json). The alternating empty-backbone inputs 10
and 19 now have independently checked gas and MeCN electronic energies in the
separate smooth-surface model. The exact reports are
[input 10](alternating__occ_none__conf010.json) and
[input 19](alternating__occ_none__conf019.json).

The archived physical reviewer ran on the original worker host, where the exact
PID/start ownership check uses the correct namespace. It checked the scored
geometry, quantum results, solvent-model/runtime definition, exclusive-residency
closure and sampled resource bounds. Collector verification then checked all
1,520 source bindings per state, retained original source hashes, and compared
reported energies directly with the bound raw electronic records. The observed
peak host RSS values were 80.46 and 80.68 GiB; both are below their predeclared
128 GiB ceilings and above the previous unsuccessful 80 GiB bound.

An initial review attempt from the collector node failed on a remote original
interpreter path. Running the unchanged review on the producing host resolved
the ownership-context problem. A subsequent shared-filesystem visibility delay
was handled by bounded reads and verified reuse of the completed report, without
overwriting it or repeating SCF. The ten collection/host-context tests are in
[tests.xml](tests.xml).

These are accepted **electronic energies only**. Their exact-coordinate native
Hessians and free energies remain separate follow-up. They do not complete the
four missing alternating empty-backbone components or establish any conformer
ensemble, grafting cycle, association free energy or electronic coupling. All
values remain outside the primary-model PMDI cycle.
