import json
from pathlib import Path
from typing import Union, Optional

PathLike = Union[str, Path]

DEFAULT_CONFIG_FILE = Path('..') / 'config' / 'default_test_config.json'

def read_server_config(config_fpath: PathLike, required_keys: dict = dict()) -> dict:
    used_fpath = Path(config_fpath)
    if not used_fpath.exists():
            raise FileNotFoundError(f"Config file {used_fpath} not found")

    with used_fpath.open('r') as f:
             CONFIG=f.read()

    if isinstance(CONFIG, str):
            CONFIG=json.loads(CONFIG)   # assume JSON string; convert.

    # check CONFIG is a dictionary  
    if not isinstance(CONFIG, dict):
            raise KeyError(f"CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {type(CONFIG)}")
    
    # check that the required keys are in config file.
    missing_keys=required_keys-set(CONFIG.keys())
    if len(missing_keys) > 0:
        raise KeyError(f"Keys: {missing_keys} are required in {config_fpath} but were not found.")

    return CONFIG

