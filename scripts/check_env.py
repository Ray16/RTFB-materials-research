"""Smoke-test the `redox` environment."""
import importlib, sys

MODS = ["torch", "fairchem", "ase", "rdkit", "pyscf", "geometric", "numpy", "pandas",
        "matplotlib", "graphviz"]
# Optional: GPU-only DFT backend. Import needs a GPU, so don't hard-fail on CPU nodes.
OPTIONAL = ["gpu4pyscf.dft"]

def main():
    ok = True
    for m in MODS:
        try:
            mod = importlib.import_module(m)
            ver = getattr(mod, "__version__", "?")
            print(f"  [ok] {m:10s} {ver}")
        except Exception as e:
            ok = False
            print(f"  [--] {m:10s} MISSING: {e}")
    for m in OPTIONAL:
        try:
            importlib.import_module(m)
            print(f"  [ok] {m:14s} (optional GPU DFT backend)")
        except Exception as e:
            print(f"  [..] {m:14s} unavailable (optional): {str(e)[:60]}")
    try:
        import torch
        print(f"  cuda available: {torch.cuda.is_available()}  "
              f"devices: {torch.cuda.device_count()}")
    except Exception:
        pass
    # graphviz's python binding needs the system `dot` executable to render.
    import shutil
    dot = shutil.which("dot")
    print(f"  dot binary: {dot or 'MISSING (pipeline flowchart cannot render)'}")
    print("ENV OK" if ok else "ENV INCOMPLETE")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
