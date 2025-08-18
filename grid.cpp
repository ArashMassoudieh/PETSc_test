#include "grid.h"
#include <algorithm> // std::transform
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <random>
#include <limits>
#include "Matrix_arma.h"
#include "Vector_arma.h"
#include <fstream>
#include <iomanip>

// If you want PETSc interop, include your wrappers here:
#include "petscvector.h" // adjust include path to your project

Grid2D::Grid2D(int nx, int ny, double Lx, double Ly)
    : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly)
{
    assert(nx_ >= 2 && ny_ >= 2);
    dx_ = Lx_ / (nx_);
    dy_ = Ly_ / (ny_);
}

void Grid2D::forEachCell(const std::function<void(int,int)>& f) const {
    for (int j=0; j<ny_; ++j)
        for (int i=0; i<nx_; ++i)
            f(i,j);
}

void Grid2D::forEachInterior(const std::function<void(int,int)>& f) const {
    for (int j=1; j<ny_-1; ++j)
        for (int i=1; i<nx_-1; ++i)
            f(i,j);
}

void Grid2D::forEachFx(const std::function<void(int,int)>& f) const {
    for (int j=0; j<ny_; ++j)
        for (int i=1; i<nx_; ++i)
            f(i,j);
}

void Grid2D::forEachFy(const std::function<void(int,int)>& f) const {
    for (int j=1; j<ny_; ++j)
        for (int i=0; i<nx_; ++i)
            f(i,j);
}

void Grid2D::computeFacePropertyFromCell(const std::string& KcellName,
                                         std::vector<double>& KxFace,
                                         std::vector<double>& KyFace) const
{
    const auto it = fields_.find(KcellName);
    assert(it != fields_.end() && "Cell property field not found");
    const auto& Kc = it->second;

    KxFace.resize(FxSize());
    KyFace.resize(FySize());

    // vertical faces (between (i-1,j) and (i,j))
    forEachFx([&](int i, int j){
        const std::size_t I_L = cellIndex(i-1,j);
        const std::size_t I_R = cellIndex(i  ,j);
        const double k = harmonic(Kc[I_L], Kc[I_R]);
        KxFace[FxIndex(i,j)] = k;
    });

    // horizontal faces (between (i,j-1) and (i,j))
    forEachFy([&](int i, int j){
        const std::size_t I_D = cellIndex(i,j-1);
        const std::size_t I_U = cellIndex(i,j  );
        const double k = harmonic(Kc[I_D], Kc[I_U]);
        KyFace[FyIndex(i,j)] = k;
    });
}

// --------- PETSc interop (uses your PETScVector wrapper) ---------

void Grid2D::packFieldToPETSc(const std::string& fieldName, PETScVector& v) const {
    const auto it = fields_.find(fieldName);
    assert(it != fields_.end() && "Field not found");
    const auto& F = it->second;

    // v must already be created with global size = numCells()
    // We use VecSetValues with ADD/INSERT? The wrapper has setValue/assemble.
    // Simpler: set by indices [0..N-1].
    v.set(0.0);
    for (std::size_t I = 0; I < F.size(); ++I)
        v.setValue(static_cast<PetscInt>(I), static_cast<PetscScalar>(F[I]), INSERT_VALUES);
    v.assemble();
}

void Grid2D::unpackFieldFromPETSc(const PETScVector& v, const std::string& fieldName) {
    auto& F = field(fieldName); // creates if missing
    F.resize(numCells());

    // Pull values [0..N-1]
    const std::size_t N = F.size();
    std::vector<PetscInt> idx(N);
    for (std::size_t I=0; I<N; ++I) idx[I] = static_cast<PetscInt>(I);
    std::vector<PetscScalar> vals(N);
    // We rely on your PETScVector::raw() + VecGetValues in your project;
    // doing it directly here to keep dependency minimal:

    Vec rv = v.raw();
    PetscErrorCode ierr = VecGetValues(rv, static_cast<PetscInt>(N), idx.data(), vals.data());
    (void)ierr; // in your project, wrap with PetscCallAbort
    for (std::size_t I=0; I<N; ++I) F[I] = static_cast<double>(vals[I]);
}

// separable exponential covariance with unit variance (sill = 1)
static inline double cov_exp_sep(double dx, double dy, double lx, double ly) {
    // Guard tiny length scales
    const double eps = 1e-15;
    lx = (lx > eps ? lx : eps);
    ly = (ly > eps ? ly : eps);
    return std::exp(-std::abs(dx)/lx - std::abs(dy)/ly);
}

// metric consistent with the covariance for nearest-neighbor sorting
static inline double metric_sep(double dx, double dy, double lx, double ly) {
    const double eps = 1e-15;
    lx = (lx > eps ? lx : eps);
    ly = (ly > eps ? ly : eps);
    return std::abs(dx)/lx + std::abs(dy)/ly;
}

void Grid2D::makeGaussianFieldSGS(const std::string& name,
                                  double lx, double ly,
                                  int nneigh,
                                  unsigned seed,
                                  double nugget,
                                  int max_ring)
{
    assert(lx > 0.0 && ly > 0.0);
    assert(nneigh >= 0);
    const std::size_t N = numCells();

    // Allocate/clear field
    auto &F = fields_[name];
    F.assign(N, 0.0);

    // Flags: has this node been simulated yet?
    std::vector<unsigned char> known(N, 0);

    // Coordinates lambdas
    auto I = [&](int i,int j) -> std::size_t { return static_cast<std::size_t>(j)*nx_ + i; };
    auto ij = [&](std::size_t idx, int &i, int &j) {
        i = static_cast<int>(idx % nx_);
        j = static_cast<int>(idx / nx_);
    };
    auto xy = [&](int i,int j, double &x, double &y) {
        x = (i+0.5) * dx_; y = (j+0.5) * dy_;
    };

    // Random visiting order (path)
    std::vector<std::size_t> order(N);
    for (std::size_t k=0; k<N; ++k) order[k] = k;
    std::mt19937_64 path_rng(seed ^ 0xA57B'cafe'1234ULL);
    std::shuffle(order.begin(), order.end(), path_rng);

    // GSL RNG for normals
    gsl_rng_env_setup();
    gsl_rng *grng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(grng, seed);

    // Seed point: first in random order — unconditional standard normal
    {
        const std::size_t s = order[0];
        F[s] = gsl_ran_gaussian(grng, 1.0);  // mean 0, std 1
        known[s] = 1;
    }

    // Neighbor search: expand square "rings" around target until we collect >= nneigh
    auto gather_neighbors = [&](int it, int jt,
                                std::vector<std::size_t>& nidx,
                                std::vector<double>&      ndist)
    {
        nidx.clear(); ndist.clear();
        if (nneigh == 0) return;

        // optional cap on expansion
        const int ring_cap = (max_ring > 0 ? max_ring : std::max(nx_, ny_));
        int r = 1;
        double xt, yt; xy(it, jt, xt, yt);

        while ((int)nidx.size() < nneigh && r <= ring_cap) {
            const int i0 = std::max(0, it - r), i1 = std::min(nx_-1, it + r);
            const int j0 = std::max(0, jt - r), j1 = std::min(ny_-1, jt + r);

            // top edge (j=j0), left->right
            for (int i=i0; i<=i1; ++i) {
                const int j=j0;
                const std::size_t s = I(i,j);
                if (!known[s]) continue;
                double x,y; xy(i,j,x,y);
                nidx.push_back(s);
                ndist.push_back(metric_sep(x-xt,y-yt,lx,ly));
            }
            // bottom edge (j=j1), left->right (skip if same as top)
            if (j1 != j0) {
                for (int i=i0; i<=i1; ++i) {
                    const int j=j1;
                    const std::size_t s = I(i,j);
                    if (!known[s]) continue;
                    double x,y; xy(i,j,x,y);
                    nidx.push_back(s);
                    ndist.push_back(metric_sep(x-xt,y-yt,lx,ly));
                }
            }
            // left edge (i=i0), interior vertical (skip corners already added)
            for (int j=j0+1; j<=j1-1; ++j) {
                const int i=i0;
                const std::size_t s = I(i,j);
                if (!known[s]) continue;
                double x,y; xy(i,j,x,y);
                nidx.push_back(s);
                ndist.push_back(metric_sep(x-xt,y-yt,lx,ly));
            }
            // right edge (i=i1), interior vertical (skip corners)
            if (i1 != i0) {
                for (int j=j0+1; j<=j1-1; ++j) {
                    const int i=i1;
                    const std::size_t s = I(i,j);
                    if (!known[s]) continue;
                    double x,y; xy(i,j,x,y);
                    nidx.push_back(s);
                    ndist.push_back(metric_sep(x-xt,y-yt,lx,ly));
                }
            }

            // If we already have a lot, we can break
            if ((int)nidx.size() >= nneigh) break;
            ++r;
        }

        // Keep the closest nneigh (by metric consistent with covariance)
        if ((int)nidx.size() > nneigh) {
            std::vector<std::size_t> perm(nidx.size());
            for (std::size_t k=0;k<perm.size();++k) perm[k]=k;
            std::nth_element(perm.begin(), perm.begin()+nneigh, perm.end(),
                             [&](std::size_t a, std::size_t b){ return ndist[a] < ndist[b]; });
            std::vector<std::size_t> nidx2; nidx2.reserve(nneigh);
            std::vector<double>      ndist2; ndist2.reserve(nneigh);
            for (int k=0; k<nneigh; ++k) { nidx2.push_back(nidx[perm[k]]); ndist2.push_back(ndist[perm[k]]); }
            nidx.swap(nidx2); ndist.swap(ndist2);
        }
    };

    // Main SGS loop
    for (std::size_t t = 1; t < N; ++t) {
        const std::size_t s = order[t];
        if (known[s]) continue; // paranoia

        int it,jt; ij(s, it, jt);
        double xt, yt; xy(it, jt, xt, yt);

        // Collect nearest conditioned neighbors (up to nneigh)
        std::vector<std::size_t> nidx;
        std::vector<double> ndist;
        gather_neighbors(it, jt, nidx, ndist);
        const int m = static_cast<int>(nidx.size());

        double sample = 0.0;

        if (m == 0 || nneigh == 0) {
            // Unconditional if nothing to condition on yet
            sample = gsl_ran_gaussian(grng, 1.0);
        } else {
            // Build kriging system: K w = k
            CMatrix_arma K(m, m);
            CVector_arma k(m);
            CVector_arma z(m);

            // Fill K (cov among neighbors), with nugget on diag
            for (int a=0; a<m; ++a) {
                int ia, ja; ij(nidx[a], ia, ja);
                double xa, ya; xy(ia, ja, xa, ya);
                z(a) = F[nidx[a]];
                for (int b=0; b<m; ++b) {
                    int ib, jb; ij(nidx[b], ib, jb);
                    double xb, yb; xy(ib, jb, xb, yb);
                    const double cab = cov_exp_sep(xa-xb, ya-yb, lx, ly);
                    K(a,b) = (a==b) ? (cab + nugget) : cab;
                }
                // Fill k (cov to target)
                k(a) = cov_exp_sep(xa-xt, ya-yt, lx, ly);
            }

            // Solve K w = k (SPD-ish); fall back to regular solve if Cholesky fails
            CVector_arma w;
            bool ok = false;
            try {
                // Prefer Cholesky for speed/stability
                arma::mat cholU;
                if (arma::chol(cholU, K)) {
                    // solve via two triangular systems: K = U^T U
                    w = arma::solve(arma::trimatu(cholU), arma::solve(arma::trimatl(cholU.t()), k));
                    ok = true;
                }
            } catch (...) { /* ignore */ }

            if (!ok) {
                // Fallback (should rarely happen with nugget)
                w = arma::solve(K, k, arma::solve_opts::fast + arma::solve_opts::likely_sympd);
            }

            const double mu     = dotproduct(w, z);                 // kriging mean
            const double kKw    = dotproduct(k, w);
            double sigma2       = std::max(1e-14, 1.0 - kKw);      // kriging variance (unit sill)
            const double sigma  = std::sqrt(sigma2);
            const double epsN01 = gsl_ran_gaussian(grng, 1.0);
            sample = mu + sigma * epsN01;
        }

        F[s]     = sample;
        known[s] = 1;
    }

    gsl_rng_free(grng);
}


void Grid2D::writeNamedMatrix(const std::string& name, ArrayKind kind,
                              const std::string& filename,
                              char delimiter, int precision,
                              bool include_header, bool flipY) const
{
    const std::vector<double>* src = nullptr;
    int rows = 0, cols = 0;

    switch (kind) {
    case ArrayKind::Cell: {
        auto it = fields_.find(name);
        if (it == fields_.end()) throw std::runtime_error("cell field not found: " + name);
        src = &it->second;
        rows = ny_; cols = nx_;
        if (src->size() != static_cast<std::size_t>(rows * cols))
            throw std::runtime_error("size mismatch for field '" + name + "'");
    } break;

    case ArrayKind::Fx: { // vertical faces: rows=ny, cols=nx+1
        auto it = fluxes_.find(name);
        if (it == fluxes_.end()) throw std::runtime_error("flux not found: " + name);
        src = &it->second;
        rows = ny_; cols = nx_ + 1;
        if (src->size() != static_cast<std::size_t>(rows * cols))
            throw std::runtime_error("size mismatch for Fx '" + name + "' (expect ny*(nx+1))");
    } break;

    case ArrayKind::Fy: { // horizontal faces: rows=ny+1, cols=nx
        auto it = fluxes_.find(name);
        if (it == fluxes_.end()) throw std::runtime_error("flux not found: " + name);
        src = &it->second;
        rows = ny_ + 1; cols = nx_;
        if (src->size() != static_cast<std::size_t>(rows * cols))
            throw std::runtime_error("size mismatch for Fy '" + name + "' (expect (ny+1)*nx)");
    } break;
    }

    writeMatrixRaw(*src, rows, cols, name, filename, delimiter, precision, include_header, flipY);
}

void Grid2D::writeNamedMatrixAuto(const std::string& name,
                                  const std::string& filename,
                                  char delimiter, int precision,
                                  bool include_header, bool flipY) const
{
    auto itF = fields_.find(name);
    if (itF != fields_.end()) {
        if (itF->second.size() != static_cast<std::size_t>(nx_ * ny_))
            throw std::runtime_error("field '" + name + "' size mismatch");
        writeMatrixRaw(itF->second, ny_, nx_, name, filename, delimiter, precision, include_header, flipY);
        return;
    }

    auto itX = fluxes_.find(name);
    if (itX != fluxes_.end()) {
        const std::size_t sz = itX->second.size();
        if (sz == FxSize()) { writeMatrixRaw(itX->second, ny_, nx_ + 1, name, filename, delimiter, precision, include_header, flipY); return; }
        if (sz == FySize()) { writeMatrixRaw(itX->second, ny_ + 1, nx_, name, filename, delimiter, precision, include_header, flipY); return; }
        throw std::runtime_error("flux '" + name + "' unexpected size");
    }

    throw std::runtime_error("'" + name + "' not found in fields or fluxes");
}

void Grid2D::writeMatrixRaw(const std::vector<double>& data,
                            int rows, int cols,
                            const std::string& name,
                            const std::string& filename,
                            char delimiter,
                            int precision,
                            bool include_header,
                            bool flipY) const
{
    if (rows <= 0 || cols <= 0)
        throw std::runtime_error("writeMatrixRaw: nonpositive rows/cols");
    if (static_cast<std::size_t>(rows)*static_cast<std::size_t>(cols) != data.size())
        throw std::runtime_error("writeMatrixRaw: size mismatch for '" + name + "'");

    std::ofstream ofs(filename, std::ios::out | std::ios::trunc);
    if (!ofs) throw std::runtime_error("writeMatrixRaw: cannot open '" + filename + "'");

    if (include_header) {
        ofs << "# name=" << name
            << " kind=matrix"
            << " nx=" << nx_ << " ny=" << ny_
            << " Lx=" << std::setprecision(17) << Lx_
            << " Ly=" << std::setprecision(17) << Ly_
            << " dx=" << dx_ << " dy=" << dy_ << "\n";
    }

    ofs << std::setprecision(precision) << std::scientific;

    auto idxRC = [&](int r, int c) -> std::size_t {
        // stored row-major: r*cols + c
        return static_cast<std::size_t>(r)*cols + c;
    };

    if (!flipY) {
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                ofs << data[idxRC(r,c)];
                if (c + 1 < cols) ofs << delimiter;
            }
            ofs << '\n';
        }
    } else {
        for (int r = rows-1; r >= 0; --r) {
            for (int c = 0; c < cols; ++c) {
                ofs << data[idxRC(r,c)];
                if (c + 1 < cols) ofs << delimiter;
            }
            ofs << '\n';
        }
    }

    ofs.flush();
    if (!ofs) throw std::runtime_error("writeMatrixRaw: write failed for '" + filename + "'");
}

