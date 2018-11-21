#!/usr/bin/env python
#Copyright (c) 2018 Yuichi Motoyama

from __future__ import print_function
import os
import os.path
import sys
import shutil
import struct
import xml.etree.ElementTree as ET
import numpy as np

class Iotk_dat():
    def __init__(self, filename, endian=sys.byteorder):
        self.f = open(filename, 'rb')
        self.integer_nptype = {2 : np.int16, 4 : np.int32, 8 : np.int64}
        self.integer_fmt = { 2 : 'h', 4 : 'i', 8 : 'q' }
        self.float_nptype = {4 : np.float32, 8 : np.float64, 16 : np.float128}
        self.float_fmt = { 4 : 'f', 8 : 'd' }
        self.complex_nptype = {4 : np.complex64, 8 : np.complex128, 16 : np.complex256}

        if endian == 'little':
            self.endian_fmt = '<'
        elif endian == 'big':
            self.endian_fmt = '>'
        else:
            raise RuntimeError('unknown endian, {0}'.format(endian))

    def __del__(self):
        self.f.close()

    def __enter__(self):
        return self

    def __exit__(self, ex_type, ex_value, trace):
        self.f.close()

    def getattrs(self, line):
        ret = {}
        for word in line.split():
            if b'=' in word:
                w = word.split(b'=')
                ret[w[0]] = w[1].split(b'"')[1]
        return ret

    def load(self, name, rewind=False, raw=False):
        if type(name) is str:
            name = name.encode()
        if rewind:
            self.f.seek(0)
        s = self.f.readline().strip()
        while not s.startswith(b'<'+name):
            s = self.f.readline().strip()
        attrs = self.getattrs(s)
        nc = int(attrs.get(b'columns', '1'))
        nr = int(attrs[b'size'])//nc
        kind = int(attrs[b'kind'])

        self.f.read(4) # size of the last line (write)

        if raw:
            b = struct.unpack(self.endian_fmt+'i', self.f.read(4))[0] - 4
            self.f.read(4) # some flag (complex or not???)
            return b, self.f.read(b)
        else:
            if attrs[b'type'] == b'integer':
                nptype = self.integer_nptype[kind]
                strfmt = self.endian_fmt + self.integer_fmt[kind]
                bcomplex = False
            elif attrs[b'type'] == b'real':
                nptype = self.float_nptype[kind]
                strfmt = self.endian_fmt + self.float_fmt[kind]
                bcomplex = False
            elif attrs[b'type'] == b'complex':
                nptype = self.complex_nptype[kind]
                strfmt = self.endian_fmt + self.float_fmt[kind]
                bcomplex = True

            self.f.read(8) # data length and some flag (complex or not???)
            ret = np.zeros((nr,nc),nptype)
            for r in range(nr):
                for c in range(nc):
                    x = struct.unpack(strfmt, self.f.read(kind))[0]
                    if bcomplex:
                        y = struct.unpack(strfmt, self.f.read(kind))[0]
                        z = complex(x,y)
                        ret[r,c] = z
                    else:
                        ret[r,c] = x
            return ret


def qe2respack(dirname, endian=sys.byteorder):
    if endian == 'little':
        endian_fmt = '<'
    elif endian == 'big':
        endian_fmt = '>'
    else:
        raise RuntimeError('unknown endian, {0}'.format(endian))

    if not os.path.exists('dir-wfn'):
        os.mkdir('dir-wfn')

    xmlfile = os.path.join(dirname, 'data-file.xml')
    print('loading {0}'.format(xmlfile))
    tree = ET.parse(xmlfile)
    root = tree.getroot()

    child = root.find('BAND_STRUCTURE_INFO')
    num_k = int(child.find('NUMBER_OF_K-POINTS').text)
    print('num_k = {0}'.format(num_k))
    num_b = int(child.find('NUMBER_OF_BANDS').text)
    print('num_b = {0}'.format(num_b))
    eFermi = float(child.find('FERMI_ENERGY').text)
    print('eFermi = {0}'.format(eFermi))

    k_vec = np.zeros((3,num_k))

    child = root.find('EIGENVALUES')
    for i in range(num_k):
        ev = child.find('K-POINT.{0}'.format(i+1))
        k_vec[:,i] = [float(x) for x in ev.find('K-POINT_COORDS').text.strip().split()]

    child = root.find('EIGENVECTORS')
    num_Gk = [int(child.find('K-POINT.{0}'.format(i+1))
                       .find('NUMBER_OF_GK-VECTORS').text) for i in range(num_k)]

    A = np.zeros((3,3))
    child = root.find('CELL')
    celldm = float(child.find('LATTICE_PARAMETER').text)
    cc = child.find('DIRECT_LATTICE_VECTORS')
    for i in range(3):
        A[i,:] = [float(x) for x in cc.find('a{0}'.format(i+1)).text.split()]
    Ainv = np.linalg.inv(A.transpose())

    child = root.find('IONS')
    n_atoms = int(child.find('NUMBER_OF_ATOMS').text)
    atom_symbs = ['' for i in range(n_atoms)]
    atom_positions = np.zeros((3,n_atoms))
    
    for i in range(n_atoms):
        atom = child.find('ATOM.{0}'.format(i+1)).attrib
        atom_symbs[i] = atom['SPECIES'].strip()
        atom_positions[:,i] = [float(x) for x in atom['tau'].split()]

    child = root.find('PLANE_WAVES')
    Ecut_for_psi = float(child.find('WFC_CUTOFF').text)

    child = root.find('SYMMETRIES')
    n_sym = int(child.find('NUMBER_OF_SYMMETRIES').text)
    print('n_sym = {0}'.format(n_sym))
    ftau = np.zeros((3,n_sym))
    mat_sym = np.zeros((3, 3, n_sym), np.int)
    for i in range(n_sym):
        sym = child.find('SYMM.{0}'.format(i+1))
        rot = sym.find('ROTATION').text.strip().split('\n')
        for j in range(3):
            mat_sym[j,:,i] = [int(x) for x in rot[j].split()]
        ftau[:,i] = [float(x) for x in sym.find('FRACTIONAL_TRANSLATION').text.split()]

    ## end of read XML file

    k_vec = np.dot(A,k_vec) / celldm

    print('generating dir-wfn/dat.sample-k')
    with open('./dir-wfn/dat.sample-k', 'w') as f:
        f.write('{0}\n'.format(num_k))
        for i in range(num_k):
            f.write('{0} {1} {2}\n'.format(k_vec[0,i],k_vec[1,i],k_vec[2,i]))

    print('generating dir-wfn/dat.lattice')
    with open('./dir-wfn/dat.lattice', 'w') as f:
        for i in range(3):
            f.write('{0} {1} {2}\n'.format(A[i,0], A[i,1], A[i,2]))
    
    print('generating dir-wfn/dat.bandcalc')
    with open('./dir-wfn/dat.bandcalc', 'w') as f:
        f.write('{0}\n'.format(Ecut_for_psi*2.0)) # Hartree => Rydberg
        f.write('{0}\n'.format(eFermi))
        f.write('0.0\n')

    print('generating dir-wfn/dat.nkm')
    with open('./dir-wfn/dat.nkm', 'w') as f:
        for i in range(num_k):
            f.write('{0}\n'.format(num_Gk[i]))

    print('generating dir-wfn/dat.kg')
    with open('./dir-wfn/dat.kg', 'w') as f:
        for k in range(num_k):
            f.write('{0}\n'.format(num_Gk[k]))
            with Iotk_dat(os.path.join(dirname, 'K{0:0>5}/gkvectors.dat'.format(k+1)), endian=endian) as inp:
                dat = inp.load('GRID')
                nr,nc = dat.shape
                for r in range(nr):
                    for c in range(nc):
                        f.write('{0} '.format(dat[r,c]))
                    f.write('\n')

    print('generating dir-wfn/dat.wfn')
    with open('./dir-wfn/dat.wfn', 'wb') as f:
        f.write(struct.pack(endian_fmt+'i',4))
        f.write(struct.pack(endian_fmt+'i',1))
        f.write(struct.pack(endian_fmt+'i',4))
        for k in range(num_k):
            with Iotk_dat(os.path.join(dirname, 'K{0:0>5}/evc.dat'.format(k+1)), endian=endian) as inp:
                for ib in range(num_b):
                    size,dat = inp.load('evc.{0}'.format(ib+1), raw=True)
                    f.write(struct.pack(endian_fmt+'i', size))
                    f.write(dat)
                    f.write(struct.pack(endian_fmt+'i', size))

    print('generating dir-wfn/dat.eigenvalue')
    with open('./dir-wfn/dat.eigenvalue', 'w') as f:
        f.write('{0}\n'.format(num_b))
        for k in range(num_k):
            tree = ET.parse(os.path.join(dirname, 'K{0:0>5}/eigenval.xml'.format(k+1)))
            root = tree.getroot()
            child = root.find('EIGENVALUES')
            for eigval in child.text.split():
                f.write(eigval)
                f.write('\n')

    print('generating dir-wfn/dat.atom_position')
    with open('./dir-wfn/dat.atom_position', 'w') as f:
        f.write('{0}\n'.format(n_atoms))
        X = np.dot(Ainv, atom_positions)
        for i in range(n_atoms):
            f.write('{0} {1} {2} {3}\n'.format(atom_symbs[i], X[0,i], X[1,i], X[2,i]))

    print('generating dir-wfn/dat.symmetry')
    with open('./dir-wfn/dat.symmetry', 'w') as f:
        for i in range(1,1001):
            tmp = (ftau*i) - (ftau*i).round()
            if (abs(tmp)<=1.0e-7).all():
                f.write('{0}\n'.format(n_sym))
                f.write('{0}\n'.format(i))
                tau = -((ftau*i).round())
                break
        for i in range(n_sym):
            mat = np.zeros((3,3))
            mat[:,:] = mat_sym[:,:,i]
            mat[:,:] = np.linalg.inv(mat)

            for j in range(3):
                for k in range(3):
                    f.write('{0} '.format(int(mat[j,k])))
            f.write('\n')
            f.write('{0} {1} {2}\n'.format(int(tau[0,i]), int(tau[1,i]), int(tau[2,i])))

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print('Please specify the output directory where data-file.xml exists')
        print('Usage: {0} {1} output_directory_of_QE'.format(sys.executable, sys.argv[0]))
        sys.exit(1)

    datadir = sys.argv[1]

    if os.path.exists('dir-wfn'):
        print('dir-wfn exists. Original directory is saved as dir-wfn_original.')
        print('If it is not necessary, please remove it (rm -r dir-wfn_original).')
        shutil.move('dir-wfn', 'dir-wfn_original')

    qe2respack(datadir)
