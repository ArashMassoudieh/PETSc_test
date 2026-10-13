// Read-only archive inputs; writes explicit PDF/CDF outputs under a new prefix.
// args: inverse_cdf D dt nx nu t_end output_prefix lc lambda_x lambda_y
#include "grid.h"
#include "petsc_init.h"
#include "TimeSeriesSet.h"
#include <iostream>
#include <cmath>
int main(int argc,char**argv) {
 PETScInit init(argc,argv);
 int ranks=0; MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
 if(ranks!=1) { std::cerr << "Run one MPI rank per realization; use the Python worker option for concurrency.\n"; return 2; }
 if(argc<11) { std::cerr << "Usage: ensemble_probe inverse_cdf D dt nx nu t_end output_prefix lc lambda_x lambda_y [PETSc options]\n"; return 2; }
 try {
 TimeSeries<double> cdf;
 if(!cdf.readfile(argv[1])) return 3;
 const double D=std::stod(argv[2]),dt=std::stod(argv[3]),end=std::stod(argv[6]);
 const int nx=std::stoi(argv[4]),nu=std::stoi(argv[5]);
 if(!std::isfinite(D) || !std::isfinite(dt) || !std::isfinite(end) || D<0 || dt<=0 || end<dt || nx<2 || nu<2) {
  std::cerr << "Invalid diffusion, time, or grid settings.\n"; return 2;
 }
 const std::string out=argv[7];
 Grid2D g(nx,nu,3.0,1.0);
 g.VelocityCorrelationModel=Grid2D::velocity_correlation_model::exponential;
 g.setDiffusionFactor(1.0);
 const double lc=std::stod(argv[8]),lx=std::stod(argv[9]),ly=std::stod(argv[10]);
 if(!std::isfinite(lc) || !std::isfinite(lx) || !std::isfinite(ly) || lc<=0 || lx<=0 || ly<=0) {
  std::cerr << "Mixing lengths must be positive and finite.\n"; return 2;
 }
 g.setMixingParams(lc,lx,ly);
 g.SetVal("diffusion",D);g.SetVal("porosity",1);g.SetVal("c_left",1);
 g.assignConstant("qx",Grid2D::ArrayKind::Fx,0);
 g.assignConstant("qy",Grid2D::ArrayKind::Fy,0);
 auto& q=g.flux("qx");
 for(int j=0;j<nu;++j) for(int i=0;i<=nx;++i)
  q[j*(nx+1)+i]=cdf.interpol(double(j)/(nu-1));
 g.computeMixingDiffusionCoefficient();
 g.assignConstant("C",Grid2D::ArrayKind::Cell,0);
 g.setBTCLocations({0.5,1.5,2.5});
 TimeSeriesSet<double> btc;
 g.SolveTransport(end,dt,"probe_",0,"","C",&btc);
 btc.write(out+"_cdf.csv");btc.derivative().write(out+"_pdf.csv");
 } catch(const std::exception& error) { std::cerr << error.what() << "\n"; return 1; }
}
