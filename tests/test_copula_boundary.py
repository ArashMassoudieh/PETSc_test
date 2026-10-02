"""Exercise the production copula matrix/RHS assembly without a PETSc install.

A recording matrix replaces PETSc storage; the assembly statements are extracted
verbatim from grid.cpp. This tests the discretization, not the PETSc backend.
Run: python3 -m unittest discover -s tests -p 'test_copula_boundary.py'
"""
import pathlib
import subprocess
import tempfile
import unittest
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]

class CopulaBoundary(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = (ROOT / 'grid.cpp').read_text().split('void Grid2D::SolveCopulaTransport(', 1)[1]
        assembly = source.split('PETScMatrix* Amat =', 1)[1].split('// 7.  BTC bookkeeping', 1)[0]
        assembly = 'PETScMatrix* Amat =' + assembly
        rhs = source.split('double rhs = static_cast<double>(c_arr[Irow]) * inv_dt;', 1)[1].split('b_vec.setValue', 1)[0]
        cls.tmp = tempfile.TemporaryDirectory()
        path = pathlib.Path(cls.tmp.name)
        cpp = r'''
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
using PetscInt = int;
const int ADD_VALUES = 0;
int nnzPerRow(int n) { return n + 4; }
struct PETScMatrix {
 int n; std::vector<double> a;
 PETScMatrix(int n_, int, int):n(n_),a(n_*n_) {}
 void ownershipRange(int& a, int& b) { a=0;b=n; }
 void setValue(int i,int j,double v,int) { a[i*n+j]+=v; }
 void assemble() {}
};
int main(int argc,char** argv) {
 const int NX=std::atoi(argv[1]),nu=2,N=NX*nu;
 const double DX=1.0/NX, inv_dx=1/DX,inv_dx2=inv_dx*inv_dx;
 const double D=std::atof(argv[2]),inv_dt=20.0,c_left_=1.0;
 const std::vector<double> v={1.0,2.0},K_diag={0.3,0.3};
 auto K=[](int j,int k){return j==k?0.0:0.3;};
 auto gIdxXfast=[NX](int i,int j){return j*NX+i;};
'''+assembly+r'''
 std::cout<<std::setprecision(17);
 for(double value:Amat->a) std::cout<<value<<" ";
 std::cout<<"\n";
 for(int Irow=0;Irow<N;++Irow) {
   int i_x=Irow%NX,j_u=Irow/NX;double vj=v[j_u],rhs=inv_dt;
'''+rhs+r'''
   std::cout<<rhs<<" ";
 }
 delete Amat;
}
'''
        (path/'test.cpp').write_text(cpp)
        cls.exe = path/'test'
        subprocess.run(['g++','-std=c++17',str(path/'test.cpp'),'-o',str(cls.exe)],check=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_uniform_inlet_is_stationary_and_bounded(self):
        for nx in (1,2,8):
            for diffusion in (0.,0.1,1.):
                with self.subTest(nx=nx,diffusion=diffusion):
                    output=subprocess.check_output([str(self.exe),str(nx),str(diffusion)],text=True).splitlines()
                    a=np.fromstring(output[0],sep=' ').reshape(2*nx,2*nx)
                    b=np.fromstring(output[1],sep=' ')
                    np.testing.assert_allclose(a @ np.ones(2*nx),b,rtol=1e-12,atol=1e-12)
                    c=np.zeros(2*nx)
                    for _ in range(100):
                        c=np.linalg.solve(a,20*c+b-20)
                        self.assertGreaterEqual(c.min(),-1e-12)
                        self.assertLessEqual(c.max(),1+1e-12)

if __name__=='__main__': unittest.main()
