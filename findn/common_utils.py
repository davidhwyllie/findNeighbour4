import json
from pathlib import Path
from typing import Union, Optional

PathLike = Union[str, Path]
ConfigLike = Union[str, dict]

DEFAULT_CONFIG_FILE = Path('..') / 'config' / 'default_test_config.json'

def validate_server_config(config_like: ConfigLike, required_keys: dict = dict()) -> None:
    CONFIG = _enforce_config_object_type(config_like)
    _enforce_key_presence(CONFIG, required_keys)

def read_server_config(config_fpath: PathLike, required_keys: dict = dict()) -> dict:
    used_fpath = Path(config_fpath)
    if not used_fpath.exists():
            raise FileNotFoundError(f"Config file {used_fpath} not found")

    with used_fpath.open('r') as f:
             CONFIG=f.read()

    CONFIG = _enforce_config_object_type(CONFIG)
    _enforce_key_presence(CONFIG, required_keys)

    return CONFIG


def _enforce_key_presence(key_dict: dict, required_keys: dict) -> None:
    """check that the required keys are in config file"""
    missing_keys=required_keys-set(key_dict.keys())
    if len(missing_keys) > 0:
        raise KeyError(f"Keys: {missing_keys} are required but were not found.")

def _enforce_config_object_type(config_like: ConfigLike) -> dict:
    CONFIG = config_like
    if isinstance(CONFIG, str):
        CONFIG=json.loads(CONFIG)   # assume JSON string; convert.

    if not isinstance(CONFIG, dict):
            raise KeyError(f"Configuration object must be either a dictionary or a JSON string encoding a dictionary, but is {type(CONFIG)}")
    return CONFIG
