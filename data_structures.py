from dataclasses import dataclass, field

@dataclass
class ProtChain:
    id: str = ""
    sequence: str = ""
    residues: list[str] = field(default_factory=list)
    aligned_sequence: str = ""
    score: float = 0.0
    aggreprot_scores: list[float] = field(default_factory=list)
    
    
@dataclass
class ProtChainPart:
    par1: ProtChain
    par2: ProtChain
    id: str = ""
    sequence: str = ""
    residues: list[str] = field(default_factory=list)
    aggreprot_score: float = 0.0

@dataclass
class Point:
    x: float
    y: float
    z: float


@dataclass
class Atom:
    pdb_index: int = 0
    name: str = ""
    element: str = ""
    coords: Point = field(default_factory=lambda: Point(0, 0, 0))

@dataclass
class Residue:
    seq_index: int = 0  #index in sequence
    pdb_index: int = 0  #index parsed from pdb file (initial index could be != 0)
    pdb_insertion_code: str = "" #for example 42A, 42B
    res_name: str = ""
    res_char: str = ""
    atoms: list[Atom] = field(default_factory=list)

@dataclass
class PdbProtChain:
    id: str = ""
    sequence: str = ""
    residues: list[Residue] = field(default_factory=list)
    aligned_sequence: str = ""
