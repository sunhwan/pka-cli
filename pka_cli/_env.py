"""MolGpKa location + import-path setup, isolated from the prediction logic."""
import os
import os.path as osp
import sys
import types


def _default_molgpka_home():
    repo_root = osp.dirname(osp.dirname(osp.abspath(__file__)))
    return osp.join(repo_root, "vendor", "MolGpKa")


def molgpka_home():
    return os.environ.get("MOLGPKA_HOME", _default_molgpka_home())


def _stub_propack():
    """Work around a broken scipy.sparse.linalg._propack Fortran extension
    seen on some macOS builds (dlopen fails on a malformed __thread_bss
    Mach-O section). MolGpKa's GCN inference never calls scipy's SVD
    solvers, so a dummy stub is safe. Only installed if the real import
    actually fails -- a no-op everywhere else.
    """
    class _DummyFn:
        def __call__(self, *a, **k):
            raise NotImplementedError("stubbed PROPACK -- not used by MolGpKa")

    class _PropackStub(types.ModuleType):
        def __getattr__(self, name):
            if name.startswith("__"):
                raise AttributeError(name)
            return _DummyFn()

    pkg = types.ModuleType("scipy.sparse.linalg._propack")
    pkg.__path__ = []
    for name in ["_spropack", "_dpropack", "_cpropack", "_zpropack"]:
        mod = _PropackStub(f"scipy.sparse.linalg._propack.{name}")
        sys.modules[f"scipy.sparse.linalg._propack.{name}"] = mod
        setattr(pkg, name, mod)
    sys.modules["scipy.sparse.linalg._propack"] = pkg


def ensure_importable():
    """Add MolGpKa's src/ to sys.path and chdir into it -- its ionization
    site detector opens utils/smarts_pattern.tsv relative to cwd. Returns
    (molgpka_home, src_dir).
    """
    home = molgpka_home()
    src = osp.join(home, "src")
    if not osp.isdir(src):
        raise FileNotFoundError(
            f"MolGpKa not found at {home}. Run scripts/setup_molgpka.sh "
            "or set MOLGPKA_HOME to an existing checkout."
        )
    try:
        import torch_geometric  # noqa: F401
    except ImportError as e:
        if "_propack" not in str(e):
            raise
        _stub_propack()

    if src not in sys.path:
        sys.path.insert(0, src)
    os.chdir(src)
    return home, src
