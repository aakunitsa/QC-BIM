#!/sw/bw/bwpy/0.3.0/python-single/usr/bin/python
import sys
import numpy as np
import psi4
import os

ANG2BOHR = 1.889725989
BOHR2ANG = 1.0/ANG2BOHR

def parse(inp):
    geom_str = ''
    bq_list = []
    options = {}

    for line in inp:
        if 'GEOM' in line:
            print 'found geom'
            for line2 in inp:
                if 'END' in line2:
                    break
                else:
                    geom_str += line2 # GEOM IN ANGSTROM
        if 'BQ_CHARGES' in line:
            for line2 in inp:
                if 'END' in line2:
                    break
                else:
                    bq = np.fromstring(line2, sep=' ') # BQ POS ALSO ANGSTROM
                    bq_list.append(bq)
        dat = line.split()
        if len(dat) > 1:
            key = dat[0].lower()
            value = dat[-1].lower()
            options[key] = value

    geom_str += "\nsymmetry c1\nno_reorient\nno_com\n"
    print "GEOM"
    print geom_str
    print "CALCTYPE", options['calc_type']
    mol = psi4.geometry(geom_str)
    mol.update_geometry()

    calc_type = 'energy' if 'calc_type' not in options else options['calc_type']
    mult = 1 if 'mult' not in options else int(options['mult'])
    charge = 0 if 'charge' not in options else int(options['charge'])

    if mult != 1:
        mol.set_multiplicity(mult)
    if charge != 0:
        mol.set_molecular_charge(charge)

    options['bq_list'] = bq_list
    bqfield = psi4.QMMM() # Create external potential in angstrom
    for bq in bq_list:
        bqfield.extern.addCharge(bq[0], bq[1], bq[2], bq[3])

    return calc_type, mol, bqfield, options


def run_energy(mol, bqfield, options):
    basis = 'aug-cc-pvdz' if 'basis' not in options else options['basis']
    theory = 'mp2' if 'theory' not in options else options['theory']
    psi4.core.set_global_option('freeze_core', 'True')
    method_str = theory + '/' + basis
    psi4.core.set_global_option_python('EXTERN', bqfield.extern)
    psi4.core.set_global_option('BASIS', basis)
    E = psi4.energy(method_str)
    with open('en_grad.dat', 'w') as fp:
        fp.write("%24.16f\n" % E)

def run_gradient(mol, bqfield, options):
    basis = 'aug-cc-pvdz' if 'basis' not in options else options['basis']
    theory = 'mp2' if 'theory' not in options else options['theory']
    psi4.core.set_global_option('freeze_core', 'True')
    method_str = theory + '/' + basis
    psi4.core.set_global_option_python('EXTERN', bqfield.extern)
    psi4.core.set_global_option('BASIS', basis)
    grad, wfn = psi4.gradient(method_str, return_wfn=True)
    
    E = psi4.core.get_variable('CURRENT ENERGY')
    gradarr = psi4.p4util.mat2arr(grad)
    with open('en_grad.dat', 'w') as fp:
        fp.write("%24.16f\n" % E)
        for g in gradarr:
            fp.write(("%24.16f"*3+"\n") % (g[0], g[1], g[2]))

    bq_list = options['bq_list'] # feed BQ positions into Psi4 as angstrom
    with open('grid.dat', 'w') as fp:
        for bq in bq_list:
            fp.write('%16.10f%16.10f%16.10f\n' % (bq[1], bq[2], bq[3]))
    psi4.oeprop(wfn, 'GRID_FIELD')

def run_refgradient(mol, bqfield, options):
    """The function calculates SCF gradient in addition
        to MP2 or CCSD gradients and saves it the file en_grad0.dat
    """
    
    basis = 'aug-cc-pvdz' if 'basis' not in options else options['basis']
    theory = 'mp2' if 'theory' not in options else options['theory']
    psi4.core.set_global_option('freeze_core', 'True')
    method_str = theory + '/' + basis
    ref_str = 'scf/' + basis
    psi4.core.set_global_option_python('EXTERN', bqfield.extern)
    psi4.core.set_global_option('BASIS', basis)

    grad0, wfn0 = pri4.gradient(ref_str, return_wfn=True)
    E0 = psi4.core.get_variable('CURRENT ENERGY')
    grad, wfn = psi4.gradient(method_str, return_wfn=True, ref_wfn=wfn0)
    E = psi4.core.get_variable('CURRENT ENERGY')
    
    gradarr0 = psi4.p4util.mat2arr(grad0)
    gradarr = psi4.p4util.mat2arr(grad)

    with open('en_grad_ref.dat', 'w') as fp:
        fp.write("%24.16f\n" % E0)
        for g in gradarr0:
            fp.write(("%24.16f"*3+"\n") % (g0[0], g0[1], g0[2]))
    with open('en_grad_correl.dat', 'w') as fp:
        fp.write("%24.16f\n" % E)
        for g, g0 in zip(gradarr, gradarr0):
            fp.write(("%24.16f"*3+"\n") % (g[0], g[1], g[2]))

    bq_list = options['bq_list'] # feed BQ positions into Psi4 as angstrom
    with open('grid.dat', 'w') as fp:
        for bq in bq_list:
            fp.write('%16.10f%16.10f%16.10f\n' % (bq[1], bq[2], bq[3]))
    psi4.oeprop(wfn, 'GRID_FIELD')

def run_bqgradient(mol, bqfield, options):
    # This subroutine doesn't need to be used
    # Easier for Fortran to simply read grid_field.dat
    raise RuntimeError("We are not using this")
    bq_field = np.loadtxt('grid_field.dat')
    bq_list = options['bq_list']
    assert bq_field.shape == (len(bq_list), 3)
    print "@QCBIM BQFORCE"
    for chg, field_vec in zip(bq_list, bq_field):
        force = -1.0 * chg[0] * field_vec 
        print ("%24.16f"*3) % (force[0], force[1], force[2])

def run_hessian(mol, bqfield, options):
    basis = 'aug-cc-pvdz' if 'basis' not in options else options['basis']
    theory = 'mp2' if 'theory' not in options else options['theory']
    psi4.core.set_global_option('freeze_core', 'True')
    method_str = theory + '/' + basis
    psi4.core.set_global_option_python('EXTERN', bqfield.extern)
    psi4.core.set_global_option('BASIS', basis)
    hess = psi4.hessian(theory)
    print "@QCBIM HESSIAN WRITTEN TO hess.dat"
    hessarr = np.array(psi4.p4util.mat2arr(hess))
    M, N = hessarr.shape
    i = 0
    with open('hess.dat', 'w') as fp:
        for row in range(M):
            for col in range(row+1):
                fp.write('%24.16E\n' % hessarr[row, col])
                i += 1
    natom = mol.natom()
    assert i == 3*natom + (3*natom)*(3*natom-1)/2

def main(path):
    os.chdir(os.path.dirname(path))
    fname = os.path.basename(path)
    psi4.set_memory(int(3.8e9))
    with open(fname) as fp:
        calc_type, mol, bqfield, options = parse(fp) # ESP, EN, GRAD, HESS

    method = getattr(sys.modules['__main__'], 'run_%s' % calc_type)
    result = method(mol, bqfield, options)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise RuntimeError("Need input path")
    path = sys.argv[1]
    if not os.path.exists(path):
        raise RuntimeError("%s is not a valid input file" % path)
    main(path)
