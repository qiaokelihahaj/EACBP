import sys
import platform

def main():
    print(f"Python version: {sys.version}")
    print(f"Platform: {platform.platform()}")
    
    modules = [
        ("pydantic", "pydantic"),
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("scipy", "scipy"),
        ("networkx", "networkx"),
        ("yaml", "PyYAML"),
        ("anndata", "anndata"),
        ("scanpy", "scanpy"),
        ("pytest", "pytest"),
        ("eacbp", "eacbp")
    ]
    
    for mod_name, label in modules:
        try:
            mod = __import__(mod_name)
            ver = getattr(mod, "__version__", "unknown")
            print(f"[OK] {label}: version {ver}")
        except ImportError as e:
            print(f"[FAIL] {label}: {e}")

if __name__ == "__main__":
    main()
