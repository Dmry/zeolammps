from maze import Zeolite
from ase.geometry.analysis import Analysis
from ase.lattice.orthorhombic import SimpleOrthorhombic
from ase.spacegroup import crystal
from ase.build import add_vacuum

bea_zeolite = Zeolite.make('BEA') 
bea_zeolite = bea_zeolite.cap_atoms()

for key, idx in bea_zeolite.get_atom_types().items():
    if key == "framework-Si":
        for idx, atom_i in enumerate(idx):
            bea_zeolite[atom_i].charge = 2.1
            bea_zeolite[atom_i].tag = 1
    elif key == "framework-O":
        for idx, atom_i in enumerate(idx):
            bea_zeolite[atom_i].charge = -1.05
            bea_zeolite[atom_i].tag = 2

bea_zeolite = crystal(symbols=bea_zeolite, size=(3,3,3), cell=bea_zeolite.cell)
bea_zeolite.center(vacuum=5, axis=(1,1))

def write_header(f):
    f.write("""
    write_once('In Init'){
        units real                              # angstroms, kCal/mole, Daltons, Kelvin
        atom_style full                         # select column format for Atoms section
        pair_style lj/cut/coul/long 12.0        # params needed: epsilon sigma
        bond_style harmonic                     # parameters needed: k_bond, r0
        angle_style harmonic                    # parameters needed: k_theta, theta0
        kspace_style pppm 0.0001                # long-range electrostatics sum method
    }\n"""
    )

def write_atoms(f, zeolite, indentation):
    f.write(indentation)
    f.write("write('Data Atoms') {\n")
    #atomID molID atomType charge coordX coordY coordZ
    for i, atom in enumerate(zeolite):
        x, y, z = atom.position
        f.write(f'{indentation}\t$atom:{i}\t$mol:.\t@atom:{atom.symbol}\t{atom.charge}\t{x}\t{y}\t{z}\n')
    f.write(indentation)
    f.write("}\n")

def write_masses(f, zeolite, indentation):
    element_indices = zeolite.symbols.indices()
    element_mass = {element:zeolite[index[0]].mass for (element,index) in element_indices.items()}

    f.write(indentation)
    f.write("write_once('Data Masses') {\n")
    for element, mass in element_mass.items():
        f.write(f"{indentation}\t@atom:{element}\t{mass}\n")
    f.write(indentation)
    f.write("}\n")

def write_angles(f, zeolite, indentation, el1, el2):
    f.write(indentation)
    f.write("write('Data Angles') {\n")

    tag = f"{el2}{el1}{el2}"

    ana = Analysis(zeolite)
    angles = ana.get_angles(el2, el1, el2, unique=True)
    for idx, angle in enumerate(angles[0]):
        f.write(f"{indentation}\t$angle:{tag}{idx}\t@angle:{tag}\t$atom:{angle[0]}\t$atom:{angle[1]}\t$atom:{angle[2]}\n")

    f.write(indentation)
    f.write("}\n")

def write_bonds(f, zeolite, indentation, el1, el2):
    f.write(indentation)
    f.write("write('Data Bonds') {\n")

    # # bondID bondType atomID1 atomID2
    tag = f"{el1}{el2}"

    ana = Analysis(zeolite)
    bonds = ana.get_bonds(el1, el2, unique=True)
    for idx, bond in enumerate(bonds[0]):
        f.write(f"{indentation}\t$bond:{tag}{idx}\t@bond:{tag}\t$atom:{el1}{bond[0]}\t$atom:{el2}{bond[1]}\n")

    f.write(indentation)
    f.write("}\n")

def write_settings(f, zeolite, indentation):
    f.write(indentation)
    f.write("write_once('In Settings') {\n")

    f.write(f"{indentation}\tpair_coeff @atom:Si @atom:O  8.8401e-4 4.36375\n")
    f.write(f"{indentation}\tpair_coeff @atom:Si @atom:Si  5.0298e-6 5.562\n")
    f.write(f"{indentation}\tpair_coeff @atom:O @atom:O  0.1554 3.1655\n")
    f.write(f"{indentation}\tangle_coeff @angle:OSiO  109.47 0\n")

    f.write(indentation)
    f.write("}\n")

#### System.lt

def write_imports(f, indentation):
    f.write("import test.lt\n")

def write_boundaries(f, zeolite, indentation):
    cell = zeolite.cell
    f.write("write_once('Data Boundary') {\n")
    f.write(f"{indentation}0\t{cell[0][0]}\txlo\txhi\n")
    f.write(f"{indentation}0\t{cell[1][1]}\tylo\tyhi\n")
    f.write(f"{indentation}0\t{cell[2][2]}\tzlo\tzhi\n")
    f.write("}\n")

def write_zeolite(f, name, indentation):
    f.write(f"bea_zeolite = new {name}[1]")

name = "Bea_zeolite"
indentation = 1*'\t'

with open('test.lt', 'w') as f:
    f.write(name + " {")
    write_header(f)
    write_masses(f, bea_zeolite, indentation)
    write_atoms(f, bea_zeolite, indentation)
    #write_angles(f, bea_zeolite, indentation, "Si", "O")
    #write_bonds(f, bea_zeolite, indentation, "Si", "O")
    write_settings(f, bea_zeolite, indentation)
    f.write("}")

with open('system.lt', 'w') as f:
    write_imports(f, indentation)
    write_boundaries(f, bea_zeolite, indentation)
    write_zeolite(f, name, indentation)
