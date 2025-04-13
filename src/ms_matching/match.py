#match.py

#function that calculates_b_y
def calculate_b_y(*args, **kwargs):
  pass
#function for testing below
def func(a=1, b=2):

def match_ms2(*args, ms2_file=None, fasta_db=None, **kwargs):
  print('welcome to match_ms2')
  print(f'args: {args}')
  print(f'kwargs: {kwargs}')
  print(f"ms2_file: {ms2_file}")
  

  pass
# ============================
#NEW
#1. Function: make_fragments
# ===========================================
def make_fragments(peptide, types=('b', 'y'), maxcharge=1):
    """
    Generate all possible m/z for fragments of types
    `types` and of charges from 1 to `maxcharge`.
    """
    for i in range(1, len(peptide)):
        for ion_type in types:
            for charge in range(1, maxcharge + 1):
                if ion_type[0] in 'abc':
                    yield mass.fast_mass(peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    yield mass.fast_mass(peptide[i:], ion_type=ion_type, charge=charge)

# ===========================================
# 2. Function: do_match 
# ===========================================
def do_match(peptide_sequence, raw_mz, error_tol=20, error_tol_units="ppm", verbose=True, **kwargs):
    """
    Match a peptide sequence to raw m/z data.
    Parameters:
    - peptide_sequence: str
    - raw_mz: dict with "m/z array" and "intensity array"
    - error_tol: float
    - error_tol_units: 'ppm' or 'Da'
    """
    if error_tol_units not in ("ppm", "Da"):
        raise ValueError("error_tol_units must be one of ('ppm', 'Da')")

    mz_array = raw_mz.get("m/z array", [])
    if error_tol_units == "ppm":
        error_tol = error_tol * 1e-6 * sum(mz_array)/len(mz_array)

    theo_fragments = list(make_fragments(peptide_sequence))
    matched = []

    for theo_mz in theo_fragments:
        for obs_mz in mz_array:
            if abs(theo_mz - obs_mz) <= error_tol:
                matched.append((theo_mz, obs_mz))
                if verbose:
                    print(f"Matched theo {theo_mz:.4f} with obs {obs_mz:.4f}")
    return matched
