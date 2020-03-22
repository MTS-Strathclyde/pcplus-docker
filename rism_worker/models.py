from typing import List
from dataclasses import dataclass


@dataclass()
class Results:
    PC_plus_exchem: float
    closure_exchem: float
    PC_exchem: float
    PMV: float
    solute_solvent_potential_energy: float


class CalculationResponse:
    def __init__(self, converged, results=None):
        self.converged = converged
        if results is not None:
            self.results = Results(**results)
        else:
            self.results = None


@dataclass()
class CalculationRequest:
    SMILES: str
    solvent: str
    closure: str
    tolerance: float = 1.0e-5
