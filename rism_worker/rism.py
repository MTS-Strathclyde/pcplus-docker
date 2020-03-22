import tempfile
import os
import shutil
import time
import datetime
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import pathlib

CUR_DIR = pathlib.Path(__file__).parent.absolute()
SOLVENT_XVV_MAP = {
    'acetonitrile': 'acetonitrile.xvv',
    'benzene': 'benzene.xvv',
    'bromobenzene': 'bromobenzene.xvv',
    'carbtet': 'carbtet.xvv',
    'chloroform': 'chloroform.xvv',
    'CS2': 'cs2.xvv',
    'water': 'cSPCE_298.15.xvv',
    'cyclohexane': 'cyclohexane.xvv',
    'decane': 'decane_chain.xvv',
    'dichloroethane': 'dichloroethane.xvv',
    'diethylether': 'diethylether.xvv',
    'dmso': 'dmso.xvv',
    'ethylacetate': 'ethylacetate.xvv',
    'heptane': 'heptane_chain.xvv',
    'isooctane': 'isooctane_chain.xvv',
    '0.1M NaCl solution': 'nacl_0.1_298.15.xvv',
    '0.2M NaCl solution': 'nacl_0.2_298.15.xvv',
    '0.3M NaCl solution': 'nacl_0.3_298.15.xvv',
    '0.4M NaCl solution': 'nacl_0.4_298.15.xvv',
    '0.5M NaCl solution': 'nacl_0.5_298.15.xvv',
    'octanol': 'octanol.xvv',
    'olive_oil': 'olive_oil.xvv',
    'toluene': 'toluene.xvv',
    'xylene': 'xylene.xvv'
}

K_B = 1.9872041E-3  # boltzmann const in kcal/mol/K
N_A = 6.022141e23  # avogadro's constant

RUNLEAP = """source leaprc.gaff2
mol = loadmol2 "{name}.mol2"
check mol
loadamberparams "{name}.frcmod"
SaveAmberParm mol "{name}.prmtop" "{name}.incrd"
SavePdb mol "{name}.pdb"
quit
"""

MIN_SCRIPT = """Normal minimization
   &cntrl
        imin=1,      ! perform minimization
        maxcyc=200,  ! The maximum number of cycles of minimization
        drms=1e-3,   ! RMS force
        ntmin=3,     ! xmin algorithm
        ntb=0,       ! no periodic boundary
        cut=999.,    ! non-bonded cutoff
        ntpr=5       ! printing frequency
   /
"""


class Xvv(object):
    """ Wrapper around xvvfile used to compute 3d-rism pressure """
    def __init__(self, fname):
        """ Read xvvfile and set instance attributes

        Parameters
        ----------

        fname : string
            Path to a valid xvv file
        """
        self.fname = fname
        self.ngrid = None
        self.nsites = None
        self.nspecies = None
        self.temperature = None
        self.dr = None
        self.atom_names = None
        self.densities = None
        self.xvv_data = None
        self.multiplicities = None
        self.unique_sites_per_species = None
        self.total_sites_per_species = None
        self.species_densities = None
        self.normalized_densities = None
        self._read_xvvfile()
        self._compute_species_properties()

    def _read_xvvfile(self):
        with open(self.fname) as f:
            lines = f.readlines()
        tot_lines = len(lines)
        for i, line in enumerate(lines):
            line = line.split()
            if len(line) <= 1:
                continue
            if line[1] == 'POINTERS':
                data = list(map(int, lines[i + 2].split()))
                self.ngrid, self.nsites, self.nspecies = data
            if line[1] == 'MTV':
                self.multiplicities = list(map(int, lines[i + 2].split()))
            if line[1] == 'NVSP':
                self.unique_sites_per_species = list(map(int, lines[i + 2].split()))
            if line[1] == 'THERMO':
                data = lines[i + 2].split()
                self.temperature = float(data[0])  # K
                self.dr = float(data[4])  # Angstrom
            if line[1] == 'ATOM_NAME':
                data = lines[i + 2].strip()
                # split into groups of 4
                self.atom_names = [data[i:i + 4].strip() for i in range(0, len(data), 4)]
            if line[1] == 'RHOV' and len(line) == 2:
                self.densities = list(map(float, lines[i + 2].split()))
                # are there more lines with density?
                counter = 3
                while lines[i + counter].startswith(' '):
                    self.densities.extend(list(map(float, lines[i + counter].split())))
                    counter += 1
                try:
                    assert len(self.densities) == len(self.atom_names)
                except AssertionError:
                    print('Inconsistent number of densities and atom names')
                    print(self.densities)
                    print(self.atom_names)
                    raise ValueError
            if line[1] == 'XVV' and len(line) == 2:
                self.xvv_data = []
                xvv_ind = i + 2
                while xvv_ind < tot_lines and not lines[xvv_ind].startswith('%'):
                    self.xvv_data.extend(lines[xvv_ind].split())
                    xvv_ind += 1
                break
        assert len(self.xvv_data) == self.ngrid * self.nsites * self.nsites
        self.xvv_data = np.array(self.xvv_data, dtype=float)
        self.xvv_data = np.reshape(self.xvv_data, (self.ngrid, self.nsites, self.nsites), order='F')

    def _compute_species_properties(self):
        self.normalized_densities = []
        for density, multiplicity in zip(self.densities, self.multiplicities):
            self.normalized_densities.append(density / multiplicity)
        self.species_densities = []
        self.total_sites_per_species = []
        pointer = 0
        for sp_sites in self.unique_sites_per_species:
            pointer += sp_sites
            total_sites = sum(self.multiplicities[pointer - sp_sites:pointer])
            self.total_sites_per_species.append(total_sites)
            self.species_densities.append(self.normalized_densities[pointer - 1])
        assert len(self.species_densities) == self.nspecies

    def compute_3drism_pressures(self, k=0):
        """ Compute 3drism pressure using loaded xvv file.
        Uses equation 20 from the article by Sergiievskyi et al.
        (http://dx.doi.org/10.1063/1.4935065).

        Parameters
        ----------
        k : int
            Which k value to use to compute pressure. The pressure can be pretty
            sensitive to it. It is recommended to experiment with a couple of
            k values or better, plot dependency of pressure on it to see
            which value works best.

        Return
        ------
        pressures : tuple of floats
            Tuple containeing two pressures.
            First element is 3D-RISM pressure (used in PC), second element is
            3D-RISM pressure minus ideal gas pressure (used in PC+).
            Both have units of kcal/mol/A^3.
        """
        xvv_k = self.xvv_data[k, :, :]
        density_vec = np.array(self.normalized_densities)
        mult_vec = np.array(self.multiplicities)
        # Z_k from sergievskyi's article
        z_k = mult_vec / density_vec * (np.identity(self.nsites) - np.linalg.inv(xvv_k))
        z_k_sum_densities2 = np.sum(density_vec * z_k * density_vec.T)
        densities_times_sites = [
            sites * dens for sites, dens in zip(self.total_sites_per_species, self.species_densities)
        ]
        pressure = sum(densities_times_sites) - .5 * z_k_sum_densities2
        pressure = pressure * self.temperature * K_B
        ideal_pressure = sum(self.species_densities) * K_B * self.temperature
        return pressure, pressure - ideal_pressure


def smi_to_mol(smi):
    """Create a reasonable 3d geometry"""
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def prepare_calc_directory(smi):
    """Copy pdb file into the directory with the same name. If such directory
    doesn't exist it will try to create it.

    Returns
    -------
    dir_name: string
        Name of the calculation directory
    name: string
        Full path to input file without extension
    ext: string
        Extension of the file (everything after last '.').
    """

    mol = smi_to_mol(smi)
    temp_dir = tempfile.mkdtemp(prefix='rism_')
    mol_p = os.path.join(temp_dir, 'molecule.pdb')
    tmp_p = os.path.join(temp_dir, 'tmp.pdb')
    Chem.MolToPDBFile(mol, tmp_p)
    # remove any connect records
    with open(tmp_p) as f:
        lines = f.readlines()
    lines = [l for l in lines if not l.startswith('CONECT')]
    with open(mol_p, 'w') as f:
        f.writelines(lines)
    return temp_dir, mol_p[:-4], 'pdb'


def prepare_logfile(name):
    """Create a logfile which will be used throught the calculation.

    Parameters
    ----------
    name : string
        Full path to calculation directory/name

    argv : list
        Command used to start script

    Returns
    -------
    out: file object
        Returns logfile.
    """
    p, _ = os.path.split(name)
    if p == '':
        p = '.'
    log_name = '{}.log'.format(name)
    logfile = open(log_name, 'w')
    logfile.write(str(datetime.datetime.now()) + '\n')  # timestamp
    logfile.flush()
    return logfile


def generate_prmtop(name, ext, logfile, molcharge=0, multiplicity=1):
    """Generate topology file using GAFF and AM1-BCC charges, scaled by the supplied
    factor.

    Parameters
    ----------
    name : string
        Full path to molecule structure file

    ext : string
        Structre file extension (indicating structure type)

    logfile : A writable file object
        A file to which calculation std. output will be written

    molcharge : int | filename (string)
        Charge of the solute

    multiplicity : int
        Multiplicty of the solute

    Returns
    -------
    out: string
        Returns the name of prepared topology file
    """
    p, no_p_name = os.path.split(name)
    if p == '':
        p = '.'
    #Firstly we use antechamber to recognize atom and bonding types, and
    #generate topology
    chg = [
        '-c',
        'bcc',  #charge method  (AM1-BCC)
        '-nc',
        str(molcharge),  #Net molecule charge
        '-m',
        str(multiplicity)  #Multiplicity
    ]
    ante_out = subprocess.check_output(
        [
            'antechamber',
            '-i',
            '{}.{}'.format(no_p_name, ext),
            '-fi',
            ext,
            '-o',
            '{}.mol2'.format(no_p_name),  # output file
            '-fo',
            'mol2',  # output format describing each residue
            '-at',
            'amber',  # atom types (gaff2)
            '-s',
            '2'  # status info ; 2 means verbose
        ] + chg,
        cwd=p)
    logfile.write(ante_out.decode('utf8'))
    #Run parmchk to generate missing gaff force field parameters
    parm_out = subprocess.check_output(
        ['parmchk2', '-i', '{}.mol2'.format(no_p_name), '-f', 'mol2', '-o', '{}.frcmod'.format(no_p_name)
         ],  #file with missing FF params
        cwd=p)
    logfile.write(parm_out.decode('utf8'))
    logfile.flush()
    #Run tleap to generate topology and coordinates for the molecule
    leap_input_name = os.path.join(p, 'runleap.in')
    with open(leap_input_name, 'w') as f:
        f.write(RUNLEAP.format(name=no_p_name))
    leap_out = subprocess.check_output(['tleap', '-f', 'runleap.in'], cwd=p)
    logfile.write(leap_out.decode('utf8'))
    logfile.flush()
    prmtop_name = '{}.prmtop'.format(no_p_name)
    return prmtop_name


def prepare_prmtop(name, ext, dir_name, logfile):
    """ Places appropriate prmtop file into the calculation folder and scales
    it's charges if necessary.

    Parameters
    ----------
    name : string
        Calculation name
    ext : string
        Structure type (extension).
    dir_name : string
        Name of calculation directory
    logfile : File_object
        Calculation log

    Returns
    -------

    out : string
        Path to prmtop file.

    """
    # ASSUME chg=0, multiplicity=1
    chg = 0
    multiplicity = 1
    return generate_prmtop(name, ext, logfile, chg, multiplicity)


class RISM3D_Singlpnt(object):
    """ A class used to assist setting up 3D-RISM calculation.

    Init is used to specify temperature as well as non-standard names for
    topology (prmtop) or water susceptibility (xvv) files.

    The calculation details like closure or tolerance are defined in
    setup_calculation method.
    """
    def __init__(self, name, T, logfile, prmtop_name=None, xvv=None):
        """ Create a class for running rism3d.snglpnt.

        Parameters
        ----------
        name : string
            Full path to pdb file without extension

        T : float
            A calculation temperature

        logfile : A writable file object
            A file to which calculation std. output will be written

        prmtop_name : string, default None
            A name of topology file for calculation. If
            it is not specified defaults to name.prmtop. Path should be given
            relative to the directory in which pdb file is located.

        xvv : string, default None
            A name of susceptibility file for this calculation. Defaults
            to water_T.xvv, where T is calculation temperature rounded to
            two digits. Path should be given relative to the directory in
            which pdb file is located.

        """
        self.name = name
        self.T = T
        self.p, self.no_p_name = os.path.split(name)
        self.logfile = logfile
        if prmtop_name:
            self.prmtop_name = prmtop_name
        else:
            self.prmtop_name = '{}.prmtop'.format(self.no_p_name)
        self.xvv_name = xvv
        self.run_flags_list = None

    def setup_calculation(self,
                          closure='pse3',
                          no_asymp_cutoff=True,
                          write_g=False,
                          write_h=False,
                          write_c=False,
                          write_u=False,
                          write_asymp=False,
                          noasympcorr=False,
                          buffer_distance=20,
                          solvbox=False,
                          grdspc=(0.5, 0.5, 0.5),
                          tolerance=1e-5,
                          polar_decomp=False,
                          verbose=0,
                          maxstep=500,
                          rism3d_path='rism3d.snglpnt'):
        """ Setup calculation rism3d.snglpnt. calculation.

        More details on each of the parameter can be found in AmberTools
        manual RISM section.

        Parameters
        ----------
        no_asymp_cutoff: bool, default True
            Do not set asympKSpaceTolerance parameter. With it water should
            be more accurate, but it breaks everything else.

        closure : string, default hnc
            Allowed closure values are kh, hnc, pseN. Here N is an
            integer.

        write_g : boolean, default False
            Specifies whether program will write radial distribution
            functions.

        write_h : boolean, default False
            Specifies whether program will write total correlation
            functions in k space.

        write_c : boolean, default False
            Specifies wheter program will write direct correlation
            functions.

        write_u : boolean, default False
            Specifies wheter program will write potential energy
            grid.

        write_asymp : boolean, default False
            Write asymptotics of total and direct correlation fuctions in
            real space.

        noasympcorr : boolean, default False
            Don't use long range corrections to compute thermodynamics.

        buffer_distance : float, default 25.0
            Minimum distance between the solute and the edge of
            the solvent box in A.

        solvbox : array-like (should contain 3 floats)
            Size of the box in x y and z directions. Overrides buffer_distance.

        grdsp: array-like (should contain 3 floats), default (0.5, 0.5, 0.5)
            Comma separated linear grid spacings for x, y and z dimensions.

        tolerance: float, default 1e-10
            Maximum residual values for solution convergence.

        polar_decomp: boolean, default False
            Decomposes solvation free energy into polar and non-polar
            components

        verbose: int, default 0
            Either 0, 1 or 2. Determines verbosity of caluclation.

        maxstep: int, default 1000
            Number of iterations in 3D-RISM calculation.

        rism3d_path : str, default rism3d.snglpnt
            Absolute path or exact name of rism3d.snglpnt program
        """
        grdspc = ','.join(map(str, grdspc))
        if solvbox:
            solvbox = ','.join(map(str, solvbox))
        self.run_flags_list = [
            rism3d_path,
            '--pdb',
            '{}.pdb'.format(self.no_p_name),
            '--prmtop',
            self.prmtop_name,
            '--rst',
            '{}.incrd'.format(self.no_p_name),
            '--xvv',
            self.xvv_name,
            '--grdspc',
            grdspc,
        ]
        self.run_flags_list.extend(['--tolerance', str(tolerance)])
        self.run_flags_list.extend(['--closure', closure])
        if no_asymp_cutoff:
            self.run_flags_list.extend(['--asympKSpaceTolerance', '0'])
        if solvbox:
            self.run_flags_list.extend(['--solvbox', solvbox])
        else:
            self.run_flags_list.extend(['--buffer', str(buffer_distance)])
        if write_g:
            self.run_flags_list.extend(['--guv', 'g_{}'.format(self.no_p_name)])
        if write_h:
            self.run_flags_list.extend(['--huv', 'h_{}'.format(self.no_p_name)])
        if write_c:
            self.run_flags_list.extend(['--cuv', 'c_{}'.format(self.no_p_name)])
        if write_u:
            self.run_flags_list.extend(['--uuv', 'u_{}'.format(self.no_p_name)])
        if write_asymp:
            self.run_flags_list.extend(['--asymp', 'a_{}'.format(self.no_p_name)])
        if noasympcorr:
            self.run_flags_list.extend(['--noasympcorr'])
        if polar_decomp:
            self.run_flags_list.extend(['--polarDecomp'])
        if verbose:
            self.run_flags_list.extend(['--verbose'])
            self.run_flags_list.extend(['{}'.format(verbose)])
        if maxstep:
            self.run_flags_list.extend(['--maxstep'])
            self.run_flags_list.extend(['{}'.format(maxstep)])

    def run_calculation_and_log(self):
        """Run 3D-RISM single point calculation and log.
        """
        start_time = time.time()
        #print(self.run_flags_list)
        self.logfile.write('3D-RISM command: {}\n'.format(self.run_flags_list))
        rism_out = subprocess.check_output(self.run_flags_list, cwd=self.p)
        self.logfile.write(rism_out.decode('utf8'))
        self.logfile.flush()
        #write timestamp and close
        end_time = time.time()
        self.logfile.write(str(datetime.datetime.now()) + '\n')
        runtime = end_time - start_time
        self.logfile.write('3D-RISM runtime: {:.0f}'.format(runtime))
        self.logfile.flush()
        self.logfile.close()


def parse_results(name, xvv_obj):
    """ Parses log file and writes free energies and corrections to
    results.txt.

    Parameters
    ----------
    name : string
        Full path to pdb file without extension

    xvv_obj : Xvv class instance
        Wrapper around xvv file used for calculation
    """
    p, _ = os.path.split(name)
    log_name = '{}.log'.format(name)
    exchem = None
    pmv = None
    uv = None
    with open(log_name, 'r') as f:
        for line in f:
            if line.startswith('rism_excessChemicalPotential'):
                exchem = float(line.split()[1])
            if line.startswith("rism_partialMolarVolume"):
                pmv = float(line.split()[1])
            if line.startswith("rism_solventPotentialEnergy"):
                uv = float(line.split()[1])
    if not pmv:
        raise ValueError("Cannot find pmv value in log file. Most likely calculation didn't converge.")
    # compute PC
    pres, pres_plus = xvv_obj.compute_3drism_pressures()  # [kcal/mol/A^3]
    PC = exchem - pres * pmv  # pressure correction [kcal/mol]
    PC_plus = exchem - pres_plus * pmv  # pressure correction plus [kcal/mol]
    return {
        'PC_plus_exchem': PC_plus,
        'closure_exchem': exchem,
        'PC_exchem': PC,
        'PMV': pmv,
        'solute_solvent_potential_energy': uv,
    }


def run_rism(calc_params):
    logfile = None
    dir_name = None
    xvv_name = SOLVENT_XVV_MAP[calc_params.solvent]
    xvv = os.path.join(CUR_DIR, 'solvent_files', xvv_name)
    try:
        dir_name, name, ext = prepare_calc_directory(calc_params.SMILES)
        logfile = prepare_logfile(name)
        prmtop_name = prepare_prmtop(name, ext, dir_name, logfile)
        xvv_obj = Xvv(xvv)
        rism_calc = RISM3D_Singlpnt(name, xvv_obj.temperature, logfile, prmtop_name=prmtop_name, xvv=xvv)
        rism_calc.setup_calculation(closure=calc_params.closure, tolerance=calc_params.tolerance)
        rism_calc.run_calculation_and_log()
        results = parse_results(name, xvv_obj)
        logfile.close()
        return results
    finally:
        if logfile is not None:
            logfile.close()
        if dir_name is not None:
            shutil.rmtree(dir_name)


if __name__ == '__main__':
    from models import CalculationRequest
    # calc = CalculationRequest('CC', 'water', 'pse3', 1.0e-5)
    # print(run_rism(calc))

    # calc = CalculationRequest('O=C(C)Oc1ccccc1C(=O)O', 'water', 'pse3', 1.0e-5)
    # print(run_rism(calc))

    for solv in SOLVENT_XVV_MAP.keys():
        print(solv)
        try:
            calc = CalculationRequest('O=C(C)Oc1ccccc1C(=O)O', solv, 'pse3', 1.0e-5)
            print(run_rism(calc))
        except Exception as e:
            print(e)
            print('failed')
