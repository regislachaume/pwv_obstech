__version__ = "0.0.1"
    
from pathlib import Path
from importlib import resources

def get_resource(path: str) -> Path:

    root = resources.files(__name__)
    return Path(root) / path

