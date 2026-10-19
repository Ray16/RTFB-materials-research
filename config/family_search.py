"""Declared first search space; family membership does not prove usable electrons.

These are computational hypotheses. No exact synthesis, solubility or redox
performance is inherited from a parent compound or a family-level measurement.
"""
SCAFFOLD = "Cc1ccc(C[*:1])cc1"
REDUCTIONS = [("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)]
CORES = {
    "viologen": dict(
        template="{anchor}[n+]1ccc(-c2cc[n+]({cap})cc2)cc1",
        states=[("ox2", 2, 1, 0), ("ox1", 1, 2, -1), ("neu", 0, 1, -2)],
        family="pyridine-multi-e"),
    "pmdi": dict(template="{anchor}n1c(=O)c2cc3c(=O)n({cap})c(=O)c3cc2c1=O",
                 states=REDUCTIONS, family="imide (n-type)"),
    "ndi": dict(template="{anchor}N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N({cap})C3=O",
                states=REDUCTIONS, family="imide (n-type)"),
}
CAPS = {
    "methyl": "C", "ethyl": "CC", "npropyl": "CCC", "isopropyl": "C(C)C",
    "nbutyl": "CCCC", "isobutyl": "CC(C)C", "methoxyethyl": "CCOC",
    "methoxyethoxyethyl": "CCOCCOC", "fluoroethyl": "CCF",
}
# Anchor strings end at the bond to the redox-core N. Ether-linked grafts have a
# distinct graph and mass from direct N-benzyl attachment.
ANCHORS = {"direct_N": "[*:1]", "ether_C2": "[*:1]OCC", "ether_C3": "[*:1]OCCC"}
# One ring substitution at a time; symmetry-equivalent products are deduplicated.
RING_SUBSTITUENTS = {"H": None, "F": "F", "CN": "C#N", "OMe": "OC", "Me": "C"}
COMPARISON_PANELS = {
    "initial": (("cap_length", "direct_N", "npropyl", "H"),
                ("ring_F", "direct_N", "ethyl", "F"),
                ("ring_CN", "direct_N", "ethyl", "CN"),
                ("linker_C2", "ether_C2", "ethyl", "H")),
    "donor_rings": (("ring_OMe", "direct_N", "ethyl", "OMe"),
                    ("ring_Me", "direct_N", "ethyl", "Me")),
}
ROUTE_STATUS = ("Exact precursor and graft route unverified. Secondary-alkyl viologen "
                "caps and mono-N-H imide preparation require particular route review.")
