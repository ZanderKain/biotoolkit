from bio_structs import *

class bio_seq:
    """DNA sequence class. Default value :ATCG, DNA, No label"""
    
    def __init__(self, seq="ATCG",seq_type="DNA",label="No label") -> None:
        self.seq = seq.upper()
        self.seq_type = seq_type
        self.label = label
        self.is_valid = self.validate()
        
        assert self.is_valid, f"The sequence inputted seems not to be valid {self.seq_type}. Please try a new one"

            
    def validate(self):
       """Checks to see if passed string is a valid DNA string"""
       return set(nucleotides).issuperset(self.seq) 