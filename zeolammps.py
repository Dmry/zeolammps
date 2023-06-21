from maze import Zeolite
from ase.geometry.analysis import Analysis
from ase.lattice.orthorhombic import SimpleOrthorhombic
from ase.spacegroup import crystal
from ase.build import add_vacuum

class Block():
    def __init__(self, header, *content):
        self.content = [line for line in content]
        self.header = header

    def get(self):
        out = [self.header + ' {']

        for content in self.content:
            if isinstance(content, Block):
                content = content.get()

            out += ['\t' + line for line in content]

        out.append('}')
        return out
    
    def write(self, location):
        with open(location, 'w') as f:
            f.write('\n'.join(self.get()) + '\n')

class FunctionBlock(Block):
    def __init__(self, func, tag, *content):
        definition = f"{func}('{tag}')"
        super().__init__(definition, *content)

    def get(self):
        return super().get()
    
    def write(self, location):
        super().write(location)

def header_block():
    content = [
        "units real                              # angstroms, kCal/mole, Daltons, Kelvin",
        "atom_style full                         # select column format for Atoms section",
        "pair_style lj/cut/coul/long 12.0        # params needed: epsilon sigma",
        "bond_style harmonic                     # parameters needed: k_bond, r0",
        "angle_style harmonic                    # parameters needed: k_theta, theta0",
        "kspace_style pppm 0.0001                # long-range electrostatics sum method"
    ]

    return FunctionBlock('write_once', 'In Init', content)


def atoms_block(zeolite):
    #atomID molID atomType charge coordX coordY coordZ
    content = []
    for i, atom in enumerate(zeolite):
        x, y, z = atom.position
        content.append(f'$atom:{i}\t$mol:.\t@atom:{atom.symbol}\t{atom.charge}\t{x}\t{y}\t{z}')

    return FunctionBlock('write', 'Data Atoms', content)

def mass_block(zeolite):
    element_indices = zeolite.symbols.indices()
    element_mass = {element:zeolite[index[0]].mass for (element,index) in element_indices.items()}
    
    content = [f"@atom:{element}\t{mass}" for element, mass in element_mass.items()]

    return FunctionBlock('write_once', 'Data Masses', content)

def angles_block(zeolite, el1, el2):
    tag = f"{el2}{el1}{el2}"

    ana = Analysis(zeolite)
    angles = ana.get_angles(el2, el1, el2, unique=True)
    content = [f"$angle:{tag}{idx}\t@angle:{tag}\t$atom:{angle[0]}\t$atom:{angle[1]}\t$atom:{angle[2]}" for idx, angle in enumerate(angles[0])]
    return FunctionBlock('write', 'Data Angles', content)

def bonds_block(zeolite, el1, el2):
    # # bondID bondType atomID1 atomID2
    tag = f"{el1}{el2}"

    ana = Analysis(zeolite)
    bonds = ana.get_bonds(el1, el2, unique=True)
    content = [f"$bond:{tag}{idx}\t@bond:{tag}\t$atom:{el1}{bond[0]}\t$atom:{el2}{bond[1]}content" for idx, bond in enumerate(bonds[0])]

    return FunctionBlock('write', 'Data Bonds', content)

def coeff_block():
    content = [
        "pair_coeff @atom:Si @atom:Si 0.0000018402 3.3019566252",
        "pair_coeff @atom:O  @atom:O 0.1554164124 3.1655200879",
        "pair_modify mix arithmetic"
    ]

    return FunctionBlock("write_once", "In Settings", content)

#### System.lt

def boundaries_block(zeolite):
    cell = zeolite.cell

    content = [
        f"0\t{cell[0][0]}\txlo\txhi",
        f"0\t{cell[1][1]}\tylo\tyhi",
        f"0\t{cell[2][2]}\tzlo\tzhi"
    ]
    
    return FunctionBlock("write_once", 'Data Boundary', content)

def set_charge_and_tag(zeolite):
    for key, idx in zeolite.get_atom_types().items():
        if key == "framework-Si":
            for idx, atom_i in enumerate(idx):
                zeolite[atom_i].charge = 2.1
                zeolite[atom_i].tag = 1
        elif key == "framework-O":
            for idx, atom_i in enumerate(idx):
                zeolite[atom_i].charge = -1.05
                zeolite[atom_i].tag = 2

    return zeolite

def main():
    bea_zeolite = Zeolite.make('BEA') 
    bea_zeolite = set_charge_and_tag(bea_zeolite)

    bea_zeolite = crystal(symbols=bea_zeolite, size=(3,3,3), cell=bea_zeolite.cell)
    bea_zeolite.center(vacuum=5, axis=(1,1))

    name = "Bea_zeolite"

    header = header_block()
    masses = mass_block(bea_zeolite)
    atoms = atoms_block(bea_zeolite)
    #angles = angles_block(bea_zeolite, "Si", "O")
    #bonds = bonds_block(bea_zeolite, "Si", "O")
    coeffs = coeff_block()

    total = Block(name, header, masses, atoms, coeffs)

    total.write('test.lt')

    with open('system.lt', 'w') as f:
        cell = bea_zeolite.cell
        comment = (f"# Unit cell vectors",
                   f"# {cell[0]}",
                   f"# {cell[1]}",
                   f"# {cell[2]}")
        f.write('\n'.join(comment))
        f.write('\n\n')
        f.write("import test.lt\n")
        boundaries = boundaries_block(bea_zeolite)
        f.write('\n'.join(boundaries.get()) + '\n')
        f.write(f"bea_zeolite = new {name}[1]")


if __name__ == "__main__":
    main()
