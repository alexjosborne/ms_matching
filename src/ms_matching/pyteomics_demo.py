#pyteomics_demo.py
from pyteomics import parser, mass
import matplotlib.pyplot as plt

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
#peptide = filtered_peptides[0][0]
#I previously had the line above but got an error so I removed it and added the line below
peptide, peptide_mass = filtered_peptides[0]
aa_masses = mass.std_aa_mass
b_ions = []
y_ions = []

# this is manual calculation
for i in range(1, len(peptide)):
    b_mass = sum(aa_masses[aa] for aa in peptide[:i]) + mass.calculate_mass(formula='H')
    b_ions.append((f'b{i}', b_mass))

for i in range(1, len(peptide)):
    y_mass = sum(aa_masses[aa] for aa in peptide[-i:]) + mass.calculate_mass(formula='H2O') + mass.calculate_mass(formula='H')
    y_ions.append((f'y{i}', y_mass))

# we can also use pyteomics.mass.fast_mass to perform all of these calculations for us
def fragments(peptide, types=('b', 'y'), maxcharge=1): # from https://pyteomics.readthedocs.io/en/latest/examples/example_msms.html
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxharge`.
    """
    for i in range(1, len(peptide)):
        for ion_type in types:
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.fast_mass(
                            peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    yield mass.fast_mass(
                            peptide[i:], ion_type=ion_type, charge=charge)
#from chatgpt
print(f"\nPeptide: {peptide}")
print(f"Monoisotopic Mass: {peptide_mass:.4f} Da")
    
# ============================
# 5. Plot spectrum nicely
# ============================
# Separate b and y
b_labels = [ion for ion, _ in b_ions]
b_mzs = [mz for _, mz in b_ions]
# b_mzs_2 = fragments(peptide, 'b') # this should match b_mzs
y_labels = [ion for ion, _ in y_ions]
y_mzs = [mz for _, mz in y_ions]

plt.figure(figsize=(12, 5))

# Plot stems
markerline, stemlines, baseline = plt.stem(b_mzs, [1]*len(b_mzs), basefmt=" ", linefmt="red")
plt.setp(markerline, visible=False)  # Hide dots
markerline, stemlines, baseline = plt.stem(y_mzs, [1]*len(y_mzs), basefmt=" ", linefmt="blue")
plt.setp(markerline, visible=False)  # Hide dots

# Annotate ions above each stem
# Annotate
for mz, label in zip(b_mzs, b_labels):
    plt.text(mz, 1.02, label, rotation=90, ha='center', va='bottom', fontsize=9, color='red')

for mz, label in zip(y_mzs, y_labels):
    plt.text(mz, 1.02, label, rotation=90, ha='center', va='bottom', fontsize=9, color='blue')



plt.xlabel('m/z')
plt.ylabel('Relative Intensity')
plt.title(f"Theoretical b- and y-ions for peptide: {peptide}")
plt.ylim(0, 1.2)
plt.tight_layout()
plt.show()

#to save the plot after its generated
plt.savefig("spectrum.png", dpi=300)

