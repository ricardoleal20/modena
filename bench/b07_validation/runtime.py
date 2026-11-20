"""The complete environment and pre-import flags for numerical processes."""
PYTHON_FLAGS = ["-I", "-S", "-B", "-X", "pycache_prefix=/dev/null"]
SAFE_ENV = {"PATH": "/usr/bin:/bin:/usr/sbin:/sbin", "LC_ALL": "C", "LANG": "C", "TZ": "UTC"}
THREAD_ENV = {"RAYON_NUM_THREADS":"1", "RUST_TEST_THREADS":"1", "OMP_NUM_THREADS":"1",
              "OPENBLAS_NUM_THREADS":"1", "MKL_NUM_THREADS":"1", "VECLIB_MAXIMUM_THREADS":"1"}


def numerical_env(registered):
    # No ambient variables are inherited. Loader/plugin inputs are absent by construction.
    allowed=set(SAFE_ENV)|set(THREAD_ENV)|{"MODENA_ASE_PYTHON"}
    if set(registered)-allowed:
        raise ValueError("unregistered code-loading environment")
    if any(registered.get(k,v)!=v for k,v in dict(SAFE_ENV,**THREAD_ENV).items()):
        raise ValueError("changed numeric or locale environment")
    return dict(SAFE_ENV,**THREAD_ENV,**{k:v for k,v in registered.items() if k=="MODENA_ASE_PYTHON"})
