__version__ = "0.0.1"

def get_resource(path):

    from importlib import resources
    from pathlib import Path

    root = resources.files(__name__)
    return Path(root) / path

