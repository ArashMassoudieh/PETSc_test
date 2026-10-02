"""Integration checks; set PETSC_TEST_EXECUTABLE to a built script2 executable.

MPI must be allowed to initialize. Requires NumPy, mpirun, and two process slots.
Run: PETSC_TEST_EXECUTABLE=/path/to/PETSc_test python3 -m unittest discover -s tests
"""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np

EXE=os.environ.get('PETSC_TEST_EXECUTABLE')

@unittest.skipUnless(EXE, 'Set PETSC_TEST_EXECUTABLE for PETSc integration checks')
class ExecutableRegressions(unittest.TestCase):
    def run_script(self, folder, script, processes=1):
        path=folder/'test.uscript2';path.write_text(script)
        command=[EXE, f'--script={path}']
        if processes>1: command=['mpirun','-np',str(processes),*command]
        return subprocess.run(command,cwd=folder,capture_output=True,text=True,timeout=30)

    def test_uniform_diffusion_does_not_create_solute(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            (d/'copula.csv').write_text('1,0\n0,1\n')
            (d/'cdf.csv').write_text('0,1\n1,1\n')
            script=f'''set locations = 0.1,0.5,0.9
grid g = nx:8, ny:2, Lx:1.0, Ly:1.0
g.load_copula name=adv, file={d}/copula.csv
g.load_copula name=diff, file={d}/copula.csv
g.load_cdf name=v, file={d}/cdf.csv, u_cap=1.0
g.solve_copula_transport theta_adv=adv, theta_diff=diff, v_of_u=v, lambda_a=0.2, lambda_d=0.1, diffusion=0.1, t_end=5.0, dt=0.01, c_left=1.0, btc_x=$locations, save_btc={d}/btc.csv
'''
            result=self.run_script(d,script)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            # save_btc contains interval derivatives at midpoint times.
            # Their sum times dt reconstructs concentration; trapezoidal
            # integration would incorrectly omit half of the end intervals.
            pdf=np.loadtxt(d/'btc.csv',delimiter=',',skiprows=1)[:,1::2]
            concentration=np.cumsum(pdf,axis=0)*0.01
            self.assertGreaterEqual(concentration.min(),-1e-6)
            self.assertLessEqual(concentration.max(),1+1e-5)
            np.testing.assert_allclose(concentration[-1],1,atol=1e-4)

    def test_correlation_paths_with_two_processes(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            prefix='''grid g = nx:8, ny:4, Lx:1.0, Ly:1.0
g.create_field name=K, type=gaussian_sgs, lx=0.2, ly=0.2, stdev=1.0, seed=123
g.to_ranks input=K, output=u
'''
            command='g.extract_copula field=u, direction=x, range=0.01:0.1:3, samples=200, bins=2, save_corr='
            result=self.run_script(d,prefix+command+'corr.csv\n'+command+'nested/corr.csv\n',2)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            for name in ['corr.csv','nested/corr.csv']:
                values=np.loadtxt(d/name,delimiter=',')
                self.assertEqual(values.shape,(3,2))
                self.assertTrue(np.isfinite(values).all())
            # A real write failure must propagate to all ranks, not hang.
            (d/'blocked').write_text('regular file')
            result=self.run_script(d,prefix+command+'blocked/corr.csv\n',2)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('cannot write',result.stderr)

if __name__=='__main__': unittest.main()
