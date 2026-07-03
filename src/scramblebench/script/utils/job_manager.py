"""This file handles checkpoint manager for downstream analysis"""
import json

from enum import Enum
from typing import Optional
from pathlib import Path
from typing_extensions import Self

class CheckpointStatus(Enum):
    """enumeration of checkpoint status

    Args:
        Enum (Enum): Enum class
    """
    COMPLETED = "COMPLETED"
    PENDING = "PENDING"
    RUNNING = 'RUNNING'
    FAILED = "FAILED"


class CheckpointManager():
    """class for genbench3d checkpoint. This class is made so that we can track the progress of genbench3d.
    Although it is possible to use the filename directly to check completion, I am using this class to also
    submit only necessary files, seen in the other functions.

    P.S. I also like being fancy trying to do this once in a while :p
    """

    def __init__(self):
        """initialize class
        """
        self.checkpoint_fname = 'checkpoint.json'
        self.state = None


    def from_dict(self, state_dict: dict[str, dict]) -> Self:
        """create a save state from a provided dictionary

        Args:
            state_dict (dict[str, dict]): dictionary containing model names and checkpoint enums

        Returns:
            Self: CheckpointManager
        """

        assert isinstance(state_dict, dict)
        self.state = state_dict

        return self


    def from_json(self, json_fname: str) -> Self:
        """create a save state from provided json file, usually load after a json file has been created.

        Args:
            json_fname (str): json file containing job state

        Returns:
            Self: CheckpointManager
        """

        assert Path(json_fname).is_file() and Path(json_fname).suffix == '.json'
        self.checkpoint_fname = json_fname
        self.state = self.load_state()

        return self


    def load_state(self) -> dict[str, dict]:
        """return a dictionary from the checkpoint filename 

        Returns:
            dict[str, dict]: state as dictionary
        """
        with open(self.checkpoint_fname, 'r', encoding='utf-8') as checkpoint_f:
            return json.load(checkpoint_f)


    def save_state(self, checkpoint_fname: Optional[str]=None) -> Self:
        """write the job progress. Will return Self for chaining

        Returns:
            Self: CheckpointManager
        """
        if checkpoint_fname:
            assert Path(checkpoint_fname).suffix == '.json'
            checkpoint_destination = checkpoint_fname
        else:
            checkpoint_destination = self.checkpoint_fname

        with open(checkpoint_destination, 'w', encoding='utf-8') as checkpoint_f:
            json.dump(self.state, checkpoint_f, indent=4)

        return self
