#pyteomics_demo.py
from pyteomics import parser, mass
import matplotlib.pyplot as plt
#I added the line below as chaptgpt said, i need it to define aa_masses
from pyteomics import mass

aa_masses = mass.std_aa_mass  # dictionary of standard amino acid masses

# ============================
# 1. Define protein sequence
# ============================
protein_seq = "RELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPI"

# ============================
# 2. In-silico trypsin digestion
# ============================
peptides = parser.cleave(protein_seq, parser.expasy_rules['trypsin'])

# ============================
# 3. Calculate peptide masses
# ============================
mass_range = (200, 10000)  # Da
length_range = (7, 50)
filtered_peptides = []

for pep in peptides:
    pep_mass = mass.calculate_mass(sequence=pep)
    pep_length = len(pep)
    if (mass_range[0] <= pep_mass <= mass_range[1]) and (length_range[0] <= pep_length <= length_range[1]):
        filtered_peptides.append((pep, pep_mass))

if not filtered_peptides:
    print("\nNo peptides in mass range!")
    exit()

# ============================
# 4. Generate b- and y-ions for first filtered peptide
# ============================
#All the code below I added # to 'remove' it and replaced it with the code right below , it will go through the possible peptides together in a loop instead of one by one
for idx, (peptide, peptide_mass) in enumerate(filtered_peptides):
    print(f"\nPeptide {idx+1}: {peptide}")
    print(f"Monoisotopic Mass: {peptide_mass:.4f} Da")

    # === Fragment ion calculation ===
    b_ions = []
    y_ions = []
    for i in range(1, len(peptide)):
        b_mass = sum(aa_masses[aa] for aa in peptide[:i]) + mass.calculate_mass(formula='H')
        b_ions.append((f'b{i}', b_mass))

        y_mass = sum(aa_masses[aa] for aa in peptide[-i:]) + mass.calculate_mass(formula='H2O') + mass.calculate_mass(formula='H')
        y_ions.append((f'y{i}', y_mass))

    # === Plot each spectrum ===
    b_labels = [ion for ion, _ in b_ions]
    b_mzs = [mz for _, mz in b_ions]
    y_labels = [ion for ion, _ in y_ions]
    y_mzs = [mz for _, mz in y_ions]

    plt.figure(figsize=(12, 5))
    plt.stem(b_mzs, [1]*len(b_mzs), basefmt=" ", linefmt="red")
    plt.stem(y_mzs, [1]*len(y_mzs), basefmt=" ", linefmt="blue")

    for mz, label in zip(b_mzs, b_labels):
        plt.text(mz, 1.02, label, rotation=90, ha='center', va='bottom', fontsize=9, color='red')
    for mz, label in zip(y_mzs, y_labels):
        plt.text(mz, 1.02, label, rotation=90, ha='center', va='bottom', fontsize=9, color='blue')

    plt.xlabel('m/z')
    plt.ylabel('Relative Intensity')
    plt.title(f"b- and y-ions for peptide {idx+1}: {peptide}")
    plt.ylim(0, 1.2)
    plt.tight_layout()
    plt.show()
    #To save the plots
    plt.savefig(f"spectrum_peptide_{idx+1}.png", dpi=300)

#peptide = filtered_peptides[0][0]
#I previously had the line above but got an error so I removed it and added the line below
#I changed the line 
#peptide, peptide_mass = filtered_peptides[0]
#aa_masses = mass.std_aa_mass
#b_ions = []
#y_ions = []

# this is manual calculation
#for i in range(1, len(peptide)):
 #   b_mass = sum(aa_masses[aa] for aa in peptide[:i]) + mass.calculate_mass(formula='H')
  #  b_ions.append((f'b{i}', b_mass))

#for i in range(1, len(peptide)):
    #y_mass = sum(aa_masses[aa] for aa in peptide[-i:]) + mass.calculate_mass(formula='H2O') + mass.calculate_mass(formula='H')
   # y_ions.append((f'y{i}', y_mass))

# we can also use pyteomics.mass.fast_mass to perform all of these calculations for us
#def fragments(peptide, types=('b', 'y'), maxcharge=1): # from https://pyteomics.readthedocs.io/en/latest/examples/example_msms.html
    #"""
    #The function generates all possible m/z for fragments of types
    #`types` and of charges from 1 to `maxharge`.
    #"""
    #for i in range(1, len(peptide)):
        #for ion_type in types:
            #for charge in range(1, maxcharge+1):
                #if ion_type[0] in 'abc':
                  #  yield mass.fast_mass(
                 #           peptide[:i], ion_type=ion_type, charge=charge)
                #else:
               #     yield mass.fast_mass(
               #             peptide[i:], ion_type=ion_type, charge=charge)
                            
                            
#from chatgpt
#print(f"\nPeptide: {peptide}")
#print(f"Monoisotopic Mass: {peptide_mass:.4f} Da")
    
# ============================
# 5. Plot spectrum nicely
# ============================
# Separate b and y
#b_labels = [ion for ion, _ in b_ions]
#b_mzs = [mz for _, mz in b_ions]
# b_mzs_2 = fragments(peptide, 'b') # this should match b_mzs
#y_labels = [ion for ion, _ in y_ions]
#y_mzs = [mz for _, mz in y_ions]

#plt.figure(figsize=(12, 5))

# Plot stems
#markerline, stemlines, baseline = plt.stem(b_mzs, [1]*len(b_mzs), basefmt=" ", linefmt="red")
#plt.setp(markerline, visible=False)  # Hide dots
#markerline, stemlines, baseline = plt.stem(y_mzs, [1]*len(y_mzs), basefmt=" ", linefmt="blue")
#plt.setp(markerline, visible=False)  # Hide dots

# Annotate ions above each stem
# Annotate
#for mz, label in zip(b_mzs, b_labels):
 #   plt.text(mz, 1.02, label, rotation=90, ha='center', va='bottom', fontsize=9, color='red')

#for mz, label in zip(y_mzs, y_labels):
 #   plt.text(mz, 1.02, label, rotation=90, ha='center', va='bottom', fontsize=9, color='blue')



#plt.xlabel('m/z')
#plt.ylabel('Relative Intensity')
#plt.title(f"Theoretical b- and y-ions for peptide: {peptide}")
#plt.ylim(0, 1.2)
#plt.tight_layout()
#plt.show()

#to save the plot after its generated
#plt.savefig("spectrum.png", dpi=300)

# ============================
# NEW LOOP: Process all filtered peptides
# ============================
#for i, (peptide, peptide_mass) in enumerate(filtered_peptides):
 #   print(f"\nPeptide {i+1}: {peptide}")
  #  print(f"Monoisotopic Mass: {peptide_mass:.4f} Da")

# ============================
#NEW
#5. Function: make_fragments
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
# 6. Function: do_match 
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
