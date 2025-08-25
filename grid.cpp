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
#include "petscmatrix.h"
#include "petscvector.h"

#ifdef GRID_USE_VTK
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkXMLImageDataWriter.h>
#endif

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
    return std::exp(-std::sqrt(std::pow(dx,2)/pow(lx,2) + std::pow(dy,2)/pow(ly,2)));
}

// metric consistent with the covariance for nearest-neighbor sorting
static inline double metric_sep(double dx, double dy, double lx, double ly) {
    const double eps = 1e-15;
    lx = (lx > eps ? lx : eps);
    ly = (ly > eps ? ly : eps);
    return std::pow(dx,2)/pow(lx,2) + std::pow(dy,2)/pow(ly,2);
}

#include "grid.h"
#include <armadillo>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <random>
#include <chrono>

// ---- helper already provided earlier (member) ----
// Grid2D::NeighborSet Grid2D::gatherNeighbors(...)
// (Use the version from my previous message; it enforces "same-row if exists".)

void Grid2D::makeGaussianFieldSGS(const std::string& name,
                                  double lx, double ly,
                                  int nneigh,
                                  unsigned long seed,
                                  double nugget,
                                  int max_ring,
                                  double tol_var0)
{
    // --- guards ---
    if (nx_ <= 0 || ny_ <= 0) throw std::runtime_error("SGS: empty grid");
    if (nneigh < 0) nneigh = 0;
    if (lx <= 0.0 || ly <= 0.0) throw std::runtime_error("SGS: lx,ly must be > 0");

    // --- ensure field storage ---
    auto& F = fields_[name];
    F.assign(static_cast<std::size_t>(nx_)*ny_, 0.0);

    // --- known bitmap ---
    std::vector<unsigned char> known(F.size(), 0);

    // --- RNGs ---
    gsl_rng* grng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(grng, seed ? seed : 12345UL);

    std::mt19937_64 path_rng(seed ? seed ^ 0x9e3779b97f4a7c15ULL : 0xC001D00D);
    std::uniform_int_distribution<int> Ui(0, nx_ - 1), Uj(0, ny_ - 1);

    // --- index helpers / coordinates ---
    auto I  = [&](int i,int j)->std::size_t { return static_cast<std::size_t>(j)*nx_ + i; };
    auto xy = [&](int i,int j,double& x,double& y){ cellCenter(i,j,x,y); };

    // covariance (anisotropic exponential, Euclidean metric)
    auto cov_exp_sep = [&](double dx, double dy)->double {
        const double rx = dx / lx;
        const double ry = dy / ly;
        const double r  = std::sqrt(rx*rx + ry*ry);
        return std::exp(-r);
    };

    // --- simulation path: random permutation of all cells ---
    std::vector<std::size_t> order(F.size());
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), path_rng);

    // --- seed the first site with N(0,1) ---
    {
        std::size_t s0 = order[0];
        F[s0] = gsl_ran_gaussian(grng, 1.0);
        known[s0] = 1;
    }

    // --- iterate remaining sites ---
    for (std::size_t t = 1; t < order.size(); ++t) {
        const std::size_t s   = order[t];
        const int it = static_cast<int>(s % nx_);
        const int jt = static_cast<int>(s / nx_);

        // gather neighbors (determined cells) near (it,jt)
        NeighborSet nbrs = gatherNeighbors(it, jt, nneigh, lx, ly, max_ring, known, /*ensure_same_row=*/true);
        const int m = static_cast<int>(nbrs.idx.size());

        double xt, yt; xy(it, jt, xt, yt);

        double mu = 0.0, sigma2 = 1.0;

        if (m == 0) {
            // No conditioning points yet (should only happen very early)
            mu = 0.0; sigma2 = 1.0;
        } else {
            // Build K (m×m), k (m), z (m)
            CMatrix_arma K(m, m);
            CVector_arma k(m);
            CVector_arma z(m);

            for (int a = 0; a < m; ++a) {
                const std::size_t sa = nbrs.idx[a];
                const int ia = static_cast<int>(sa % nx_), ja = static_cast<int>(sa / nx_);
                double xa, ya; xy(ia, ja, xa, ya);

                z(a) = F[sa];
                k(a) = cov_exp_sep(xa - xt, ya - yt);

                for (int b = 0; b < m; ++b) {
                    const std::size_t sb = nbrs.idx[b];
                    const int ib = static_cast<int>(sb % nx_), jb = static_cast<int>(sb / nx_);
                    double xb, yb; xy(ib, jb, xb, yb);
                    K(a,b) = cov_exp_sep(xa - xb, ya - yb);
                }
            }

            if (nugget > 0.0) {
                for (int a = 0; a < m; ++a) K(a,a) += nugget;
            }

            // Solve K w = k   (avoid forming K^{-1})
            arma::vec w;
            bool ok = arma::solve(w, K, k, arma::solve_opts::fast + arma::solve_opts::likely_sympd);
            if (!ok) {
                // small fallback jitter if needed
                arma::mat K2 = K;
                for (int a = 0; a < m; ++a) K2(a,a) += (nugget > 0.0 ? nugget : 1e-12);
                ok = arma::solve(w, K2, k, arma::solve_opts::fast + arma::solve_opts::likely_sympd);
                if (!ok) throw std::runtime_error("SGS: solve failed (K nearly singular)");
            }

            mu     = dotproduct(w, z);
            sigma2 = std::max(0.0, 1.0 - dotproduct(k, w));  // C(0)=1 for standardized field
        }

        // Draw (or clamp)
        double val;
        if (sigma2 < tol_var0) {
            val = mu;  // deterministic if variance is (near) zero
        } else {
            const double stddev = std::sqrt(sigma2);
            val = mu + stddev * gsl_ran_gaussian(grng, 1.0);
        }

        F[s] = val;
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

#ifdef GRID_USE_VTK

static void FillVTKCellArray(vtkDoubleArray* arr,
                             const std::vector<double>& data,
                             int cx, int cy, bool flipY)
{
    arr->SetNumberOfComponents(1);
    vtkIdType nT = static_cast<vtkIdType>(cx) * static_cast<vtkIdType>(cy);
    arr->SetNumberOfTuples(nT);

    vtkIdType id = 0;
    if (!flipY) {
        for (int j = 0; j < cy; ++j) {
            const std::size_t off = static_cast<std::size_t>(j) * cx;
            for (int i = 0; i < cx; ++i, ++id)
                arr->SetValue(id, data[off + i]);
        }
    } else {
        for (int j = cy - 1; j >= 0; --j) {
            const std::size_t off = static_cast<std::size_t>(j) * cx;
            for (int i = 0; i < cx; ++i, ++id)
                arr->SetValue(id, data[off + i]);
        }
    }
}

void Grid2D::writeNamedVTI(const std::string& name,
                           ArrayKind kind,
                           const std::string& filename,
                           const std::string& arrayName,
                           bool flipY) const
{
    // Select source and cell grid size
    const std::vector<double>* src = nullptr;
    int cellsX = 0, cellsY = 0;

    switch (kind) {
    case ArrayKind::Cell: {
        auto it = fields_.find(name);
        if (it == fields_.end()) throw std::runtime_error("Field not found: " + name);
        src = &it->second;
        cellsX = nx_; cellsY = ny_;
        if (src->size() != static_cast<std::size_t>(cellsX * cellsY))
            throw std::runtime_error("Field size mismatch for " + name);
    } break;

    case ArrayKind::Fx: {
        auto it = fluxes_.find(name);
        if (it == fluxes_.end()) throw std::runtime_error("Flux not found: " + name);
        src = &it->second;
        cellsX = nx_ + 1; cellsY = ny_;
        if (src->size() != static_cast<std::size_t>(cellsX * cellsY))
            throw std::runtime_error("Fx size mismatch for " + name);
    } break;

    case ArrayKind::Fy: {
        auto it = fluxes_.find(name);
        if (it == fluxes_.end()) throw std::runtime_error("Flux not found: " + name);
        src = &it->second;
        cellsX = nx_; cellsY = ny_ + 1;
        if (src->size() != static_cast<std::size_t>(cellsX * cellsY))
            throw std::runtime_error("Fy size mismatch for " + name);
    } break;
    }

    // Points = cells + 1 (inclusive extents in POINT indices)
    const int ptsX = cellsX + 1;
    const int ptsY = cellsY + 1;

    auto img = vtkSmartPointer<vtkImageData>::New();
    img->SetOrigin(0.0, 0.0, 0.0);
    img->SetSpacing(dx_, dy_, 1.0);
    // Extent is inclusive in POINTS: to get ptsX points, use 0..ptsX-1
    img->SetExtent(0, ptsX - 1, 0, ptsY - 1, 0, 0);

    auto arr = vtkSmartPointer<vtkDoubleArray>::New();
    arr->SetName(arrayName.empty() ? name.c_str() : arrayName.c_str());
    FillVTKCellArray(arr, *src, cellsX, cellsY, flipY);
    img->GetCellData()->SetScalars(arr);

    auto w = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    w->SetFileName(filename.c_str());
    w->SetInputData(img);
    // Force ASCII to avoid any binary/appended encoding issues:
    w->SetDataModeToAscii();
    // (Optional) w->EncodeAppendedDataOff(); // not used in ASCII mode
    if (w->Write() == 0) {
        throw std::runtime_error("VTK failed to write: " + filename);
    }
}

void Grid2D::writeNamedVTI_Auto(const std::string& name,
                                const std::string& filename,
                                const std::string& arrayName,
                                bool flipY) const
{
    auto itF = fields_.find(name);
    if (itF != fields_.end()) {
        if (itF->second.size() != static_cast<std::size_t>(nx_ * ny_))
            throw std::runtime_error("Field size mismatch for " + name);
        writeNamedVTI(name, ArrayKind::Cell, filename, arrayName, flipY);
        return;
    }
    auto itX = fluxes_.find(name);
    if (itX != fluxes_.end()) {
        const auto sz = itX->second.size();
        if (sz == FxSize()) { writeNamedVTI(name, ArrayKind::Fx, filename, arrayName, flipY); return; }
        if (sz == FySize()) { writeNamedVTI(name, ArrayKind::Fy, filename, arrayName, flipY); return; }
        throw std::runtime_error("Flux '"+name+"' has unexpected size");
    }
    throw std::runtime_error("Name '"+name+"' not found in fields or fluxes");
}
#endif // GRID_USE_VTK


Grid2D::NeighborSet Grid2D::gatherNeighbors(int it, int jt,
                                            int nneigh,
                                            double lx, double ly,
                                            int max_ring,
                                            const std::vector<unsigned char>& known,
                                            bool ensure_same_row) const
{
    NeighborSet out;
    if (nneigh <= 0) return out;

    // metric consistent with your covariance (monotone in sqrt'd version)
    auto metric_sep = [&](double dx, double dy) -> double {
        const double eps = 1e-15;
        lx = (lx > eps ? lx : eps);
        ly = (ly > eps ? ly : eps);
        return (dx*dx)/(lx*lx) + (dy*dy)/(ly*ly);
    };

    // helpers
    auto I  = [&](int i,int j) -> std::size_t { return static_cast<std::size_t>(j)*nx_ + i; };
    auto xy = [&](int i,int j, double& x,double& y){ cellCenter(i,j,x,y); };

    const int ring_cap = (max_ring > 0 ? max_ring : std::max(nx_, ny_));

    double xt, yt; xy(it, jt, xt, yt);

    std::vector<std::size_t> cand;
    cand.reserve(std::min(nneigh*4, nx_*ny_)); // heuristic

    auto push_if_known = [&](int ii, int jj) {
        const std::size_t s = I(ii, jj);
        if (known[s]) cand.push_back(s);
    };

    // Does the target row contain any known cells at all?
    auto row_has_known = [&]()->bool{
        for (int ii = 0; ii < nx_; ++ii) if (known[I(ii,jt)]) return true;
        return false;
    }();

    auto has_same_row_in = [&](const std::vector<std::size_t>& v)->bool{
        for (auto s : v) { int ii = s % nx_, jj = s / nx_; if (jj == jt) return true; }
        return false;
    };

    // Expand rings; do NOT break merely on count—also ensure same-row presence if required
    for (int r = 1; r <= ring_cap; ++r) {
        const int i0 = std::max(0, it - r), i1 = std::min(nx_-1, it + r);
        const int j0 = std::max(0, jt - r), j1 = std::min(ny_-1, jt + r);

        // top and bottom edges
        for (int i=i0; i<=i1; ++i) push_if_known(i, j0);
        if (j1 != j0) for (int i=i0; i<=i1; ++i) push_if_known(i, j1);

        // left and right edges (interior only)
        for (int j=j0+1; j<=j1-1; ++j) push_if_known(i0, j);
        if (i1 != i0) for (int j=j0+1; j<=j1-1; ++j) push_if_known(i1, j);

        // stop once we have enough AND (if needed) a same-row candidate
        if ((int)cand.size() >= nneigh) {
            if (!ensure_same_row || !row_has_known || has_same_row_in(cand))
                break;
        }
    }

    // De-duplicate (ring expansion won’t normally dup within a ring, but safe)
    std::sort(cand.begin(), cand.end());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

    if (cand.empty()) return out;

    // Compute distances once
    std::vector<double> dist(cand.size());
    for (std::size_t k=0; k<cand.size(); ++k) {
        int ii = cand[k] % nx_, jj = cand[k] / nx_;
        double x,y; xy(ii,jj,x,y);
        dist[k] = metric_sep(x - xt, y - yt);
    }

    // Keep true top-n by metric (partial select over an index indirection)
    if ((int)cand.size() > nneigh) {
        std::vector<std::size_t> idx(cand.size());
        std::iota(idx.begin(), idx.end(), 0);
        auto cmp = [&](std::size_t a, std::size_t b){ return dist[a] < dist[b]; };
        std::nth_element(idx.begin(), idx.begin()+nneigh, idx.end(), cmp);

        NeighborSet sel;
        sel.idx.reserve(nneigh);
        sel.dist.reserve(nneigh);
        for (int k=0; k<nneigh; ++k) {
            sel.idx.push_back(cand[idx[k]]);
            sel.dist.push_back(dist[idx[k]]);
        }
        return sel;
    } else {
        out.idx = std::move(cand);
        out.dist = std::move(dist);
        return out;
    }
}

static inline double harm(double a, double b) {
    // harmonic average, guarding zeros
    if (a <= 0 && b <= 0) return 0.0;
    if (a <= 0) return b;
    if (b <= 0) return a;
    return 2.0*a*b/(a+b);
}


void Grid2D::DarcySolve(double H_west, double H_east, const std::string &kx_field, const std::string &ky_field, const char* ksp_prefix)
{
    if (!hasField(kx_field) || !hasField(ky_field))
        throw std::runtime_error("DarcySolve: missing Kx and/or Ky fields");
    const auto& Kx = field(kx_field);  // size nx*ny
    const auto& Ky = field(ky_field);  // size nx*ny

    const int NX = nx(), NY = ny();
    const double DX = dx(), DY = dy();

    // FV transmissibility scalings on rectangles:
    //   T_x(face) = Kx_face * (DY/DX)
    //   T_y(face) = Ky_face * (DX/DY)
    const double sx = DY / DX;
    const double sy = DX / DY;

    auto I = [&](int i,int j)->PetscInt { return static_cast<PetscInt>(j*NX + i); };
    auto kx_at = [&](int i,int j)->double { return Kx[static_cast<std::size_t>(j)*NX + i]; };
    auto ky_at = [&](int i,int j)->double { return Ky[static_cast<std::size_t>(j)*NX + i]; };

    const PetscInt N = static_cast<PetscInt>(NX*NY);
    PETScMatrix A(N, N, /*nzPerRow=*/5);
    PETScVector b(N); b.set(0.0);

    PetscInt Istart=0, Iend=0;
    A.ownershipRange(Istart, Iend);

    for (PetscInt Irow = Istart; Irow < Iend; ++Irow) {
        const int j = static_cast<int>(Irow / NX);
        const int i = static_cast<int>(Irow - j*NX);

        double diag = 0.0;
        double rhs  = 0.0;

        // ---- X faces ----
        if (i == 0) {
            // West boundary: Dirichlet at x=0 -> half-cell face to boundary.
            // Transmissibility from cell center to boundary:
            //   T_W^D = Kx(i,j) * DY / (DX/2) = 2*Kx(i,j)*(DY/DX) = 2*Kx*sx
            const double TWD = 2.0 * kx_at(i,j) * sx;
            diag += TWD * 1.0;                 // contributes 1*TWD on the diagonal twice (see below)
            rhs  += TWD * H_west * 1.0;        // RHS gets 1*TWD*H_west twice
            // NOTE: using “twice” is equivalent to (2*T_W^D)*H_west and (2*T_W^D) on diag.
            // We’ll add the second contribution right after to keep the pattern clear:
            diag += TWD;
            rhs  += TWD * H_west;

            // East interior face
            const int iE = i+1;
            const double kxE = harm(kx_at(i,j), kx_at(iE,j));
            const double TE  = kxE * sx;
            diag += TE;
            const PetscInt colE = I(iE,j);
            A.setValue(Irow, colE, -TE, ADD_VALUES);

        } else if (i == NX-1) {
            // East boundary: Dirichlet at x=DX -> half-cell face
            const double TED = 2.0 * kx_at(i,j) * sx;
            diag += TED;
            rhs  += TED * H_east;
            diag += TED;
            rhs  += TED * H_east;

            // West interior face
            const int iW = i-1;
            const double kxW = harm(kx_at(i,j), kx_at(iW,j));
            const double TW  = kxW * sx;
            diag += TW;
            const PetscInt colW = I(iW,j);
            A.setValue(Irow, colW, -TW, ADD_VALUES);


        } else {
            // Interior in x: both west and east faces exist
            const int iW = i-1, iE = i+1;
            const double kxW = harm(kx_at(i,j), kx_at(iW,j));
            const double kxE = harm(kx_at(i,j), kx_at(iE,j));
            const double TW  = kxW * sx;
            const double TE  = kxE * sx;

            diag += TW + TE;

            const PetscInt cols[2] = { I(iW,j), I(iE,j) };
            const PetscScalar vals[2] = { -TW, -TE };
            A.setValues(1, &Irow, 2, cols, vals, ADD_VALUES);
        }

        // ---- Y faces (no-flux top/bottom => omit exterior faces) ----
        if (j > 0) {
            const int jS = j-1;
            const double kyS = harm(ky_at(i,j), ky_at(i,jS));
            const double TS  = kyS * sy;
            diag += TS;
            const PetscInt colS = I(i,jS);
            A.setValue(Irow, colS, -TS, ADD_VALUES);


        }
        if (j < NY-1) {
            const int jN = j+1;
            const double kyN = harm(ky_at(i,j), ky_at(i,jN));
            const double TN  = kyN * sy;
            diag += TN;
            const PetscInt colN = I(i,jN);
            A.setValue(Irow, colN, -TN, ADD_VALUES);

        }

        // finalize diagonal & RHS for this row
        A.setValue(Irow, Irow, diag, ADD_VALUES);
        if (rhs != 0.0)
            b.setValue(Irow, rhs, ADD_VALUES);
    }

    A.assemble();
    b.assemble();

    PetscCallAbort(PETSC_COMM_WORLD, MatSetOption(A.raw(), MAT_SPD, PETSC_TRUE));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetOption(A.raw(), MAT_SYMMETRIC, PETSC_TRUE));
    PetscCallAbort(PETSC_COMM_WORLD, MatSetOption(A.raw(), MAT_SYMMETRY_ETERNAL, PETSC_TRUE));

    PETScVector h = (ksp_prefix ? A.solveNew(b, ksp_prefix) : (A / b));

    // --- store head field ---
    auto& HEAD = field("head");
    HEAD.assign(static_cast<std::size_t>(NX)*NY, 0.0);
    {
        PetscInt Istart=0, Iend=0;
        A.ownershipRange(Istart, Iend);
        for (PetscInt Irow = Istart; Irow < Iend; ++Irow) {
            PetscScalar v;
            PetscCallAbort(PETSC_COMM_WORLD, VecGetValues(h.raw(), 1, &Irow, &v));
            HEAD[static_cast<std::size_t>(Irow)] = static_cast<double>(v);
        }
    }

    // --- compute fluxes (unchanged; consistent with BCs) ---
    auto& QX = flux("qx"); QX.assign(FxSize(), 0.0);
    auto& QY = flux("qy"); QY.assign(FySize(), 0.0);
    auto head_at = [&](int i,int j)->double { return HEAD[static_cast<std::size_t>(j)*NX + i]; };

    // vertical faces
    for (int j=0; j<NY; ++j) {
        // left boundary: one-sided gradient to Dirichlet
        {
            const int i = 0;
            const double kface = kx_at(i,j);
            const double grad  = (head_at(i,j) - H_west) / (0.5*DX);
            QX[FxIndex(0, j)] = -kface * grad;
        }
        // interior faces
        for (int i=1; i<=NX-1; ++i) {
            const double kface = harm(kx_at(i-1,j), kx_at(i,j));
            const double grad  = (head_at(i,j) - head_at(i-1,j)) / DX;
            QX[FxIndex(i, j)] = -kface * grad;
        }
        // right boundary
        {
            const int i = NX-1;
            const double kface = kx_at(i,j);
            const double grad  = (H_east - head_at(i,j)) / (0.5*DX);
            QX[FxIndex(NX, j)] = -kface * grad;
        }
    }
    // horizontal faces (no-flux at bottom/top)
    for (int i=0; i<NX; ++i) {
        QY[FyIndex(i, 0)] = 0.0;
        for (int j=1; j<=NY-1; ++j) {
            const double kface = harm(ky_at(i,j-1), ky_at(i,j));
            const double grad  = (head_at(i,j) - head_at(i,j-1)) / DY;
            QY[FyIndex(i, j)] = -kface * grad;
        }
        QY[FyIndex(i, NY)] = 0.0;
    }
}

void Grid2D::createExponentialField(const std::string& inputFieldName,
                            double a, double b,
                            const std::string& outputFieldName) {
    auto it = fields_.find(inputFieldName);
    if (it == fields_.end()) {
        throw std::runtime_error("Input field not found: " + inputFieldName);
    }

    const std::vector<double>& input = it->second;
    std::vector<double> output;
    output.reserve(input.size());

    for (double x : input) {
        output.push_back(std::exp(a * x + b));
    }

    // Overwrite if the output field already exists
    fields_[outputFieldName] = std::move(output);
}

double Grid2D::interpolate(const std::string& name, ArrayKind kind,
                           double x, double y, bool clamp) const
{
    const double Lx = nx_ * dx_;
    const double Ly = ny_ * dy_;

    auto clamp01 = [](double v, double lo, double hi) {
        return (v < lo ? lo : (v > hi ? hi : v));
    };

    // Clamp or reject outside-domain queries
    if (clamp) {
        x = clamp01(x, 0.0, Lx);
        y = clamp01(y, 0.0, Ly);
    } else {
        if (x < 0.0 || x > Lx || y < 0.0 || y > Ly)
            throw std::out_of_range("interpolate(): (x,y) outside domain");
    }

    // Accessors
    auto const& A_cell = (kind == ArrayKind::Cell) ? field(name)
                         : (kind == ArrayKind::Fx)  ? flux(name)
                                                   :                            flux(name); // Fy also in flux()
    // indexers
    auto idx_cell = [&](int i,int j) -> std::size_t {
        return static_cast<std::size_t>(j)*nx_ + static_cast<std::size_t>(i);
    };
    auto idx_fx = [&](int i,int j) -> std::size_t { // (nx+1)*ny
        return static_cast<std::size_t>(j)*(nx_+1) + static_cast<std::size_t>(i);
    };
    auto idx_fy = [&](int i,int j) -> std::size_t { // nx*(ny+1)
        return static_cast<std::size_t>(j)*nx_ + static_cast<std::size_t>(i);
    };

    // Bilinear helper
    auto bilerp = [](double v00, double v10, double v01, double v11,
                     double sx, double sy) -> double {
        // sx,sy in [0,1]; (0,0)=lower-left, (1,1)=upper-right
        const double vx0 = v00*(1.0 - sx) + v10*sx;
        const double vx1 = v01*(1.0 - sx) + v11*sx;
        return vx0*(1.0 - sy) + vx1*sy;
    };

    // Compute discrete coordinates and neighbors per kind
    if (kind == ArrayKind::Cell) {
        // logical coords in cell-center lattice
        double xi = x / dx_ - 0.5;          // i + frac; cell centers
        double yi = y / dy_ - 0.5;          // j + frac

        // For bilinear, base index must allow +1 neighbor
        int i0 = static_cast<int>(std::floor(xi));
        int j0 = static_cast<int>(std::floor(yi));

        // Clamp to interior so i1<=nx_-1, j1<=ny_-1
        i0 = std::max(0, std::min(nx_ - 2, i0));
        j0 = std::max(0, std::min(ny_ - 2, j0));

        const int i1 = i0 + 1;
        const int j1 = j0 + 1;

        const double sx = xi - i0;          // in [0,1]
        const double sy = yi - j0;

        const double v00 = field(name)[idx_cell(i0, j0)];
        const double v10 = field(name)[idx_cell(i1, j0)];
        const double v01 = field(name)[idx_cell(i0, j1)];
        const double v11 = field(name)[idx_cell(i1, j1)];

        return bilerp(v00, v10, v01, v11, sx, sy);
    }
    else if (kind == ArrayKind::Fx) {
        // vertical faces at x = i*dx (i=0..nx_), y = (j+0.5)dy (j=0..ny_-1)
        double ui = x / dx_;                // face index i with fractional
        double vj = y / dy_ - 0.5;          // row index j (cell-centered in y)

        // Base indices: i in [0..nx_-1] so i+1 exists; j in [0..ny_-2] so j+1 exists
        int i0 = static_cast<int>(std::floor(ui));
        int j0 = static_cast<int>(std::floor(vj));

        i0 = std::max(0, std::min(nx_ - 1, i0));    // note: faces go to nx_, but for bilinear we keep base <= nx_-1
        j0 = std::max(0, std::min(ny_ - 2, j0));

        const int i1 = std::min(nx_, i0 + 1);       // neighbor face to the right; clamp at nx_
        const int j1 = j0 + 1;

        const double sx = ui - i0;                  // [0,1]
        const double sy = vj - j0;

        const auto& FX = flux(name);                // (nx_+1)*ny_
        const double v00 = FX[idx_fx(i0, j0)];
        const double v10 = FX[idx_fx(i1, j0)];
        const double v01 = FX[idx_fx(i0, j1)];
        const double v11 = FX[idx_fx(i1, j1)];

        return bilerp(v00, v10, v01, v11, sx, sy);
    }
    else { // ArrayKind::Fy
        // horizontal faces at x = (i+0.5)dx (i=0..nx_-1), y = j*dy (j=0..ny_)
        double ui = x / dx_ - 0.5;          // column index i (cell-centered in x)
        double vj = y / dy_;                // face row j with fractional

        int i0 = static_cast<int>(std::floor(ui));
        int j0 = static_cast<int>(std::floor(vj));

        i0 = std::max(0, std::min(nx_ - 2, i0));
        j0 = std::max(0, std::min(ny_ - 1, j0));    // base <= ny_-1

        const int i1 = i0 + 1;
        const int j1 = std::min(ny_, j0 + 1);

        const double sx = ui - i0;
        const double sy = vj - j0;

        const auto& FY = flux(name);                // nx_*(ny_+1)
        const double v00 = FY[idx_fy(i0, j0)];
        const double v10 = FY[idx_fy(i1, j0)];
        const double v01 = FY[idx_fy(i0, j1)];
        const double v11 = FY[idx_fy(i1, j1)];

        return bilerp(v00, v10, v01, v11, sx, sy);
    }
}

double Grid2D::fieldMean(const std::string& name, ArrayKind kind) const
{
    const std::vector<double>* A = nullptr;
    std::size_t N = 0;

    if (kind == ArrayKind::Cell) {
        A = &field(name);
        N = static_cast<std::size_t>(nx_) * ny_;
    } else if (kind == ArrayKind::Fx) {
        A = &flux(name);
        N = static_cast<std::size_t>(nx_ + 1) * ny_;
    } else { // ArrayKind::Fy
        A = &flux(name);
        N = static_cast<std::size_t>(nx_) * (ny_ + 1);
    }

    double sum = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        sum += (*A)[i];

    return (N > 0 ? sum / static_cast<double>(N) : 0.0);
}

double Grid2D::fieldStdDev(const std::string& name, ArrayKind kind) const
{
    const std::vector<double>* A = nullptr;
    std::size_t N = 0;

    if (kind == ArrayKind::Cell) {
        A = &field(name);
        N = static_cast<std::size_t>(nx_) * ny_;
    } else if (kind == ArrayKind::Fx) {
        A = &flux(name);
        N = static_cast<std::size_t>(nx_ + 1) * ny_;
    } else { // ArrayKind::Fy
        A = &flux(name);
        N = static_cast<std::size_t>(nx_) * (ny_ + 1);
    }

    if (N == 0) return 0.0;

    // First pass: mean
    double mean = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        mean += (*A)[i];
    mean /= static_cast<double>(N);

    // Second pass: variance
    double var = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        double d = (*A)[i] - mean;
        var += d * d;
    }
    var /= static_cast<double>(N);

    return std::sqrt(var);
}

void Grid2D::normalizeField(const std::string& name, ArrayKind kind, double a, double b)
{
    if (b == 0.0)
        throw std::runtime_error("normalizeField: divisor b cannot be zero");

    std::vector<double>* A = nullptr;
    std::size_t N = 0;

    if (kind == ArrayKind::Cell) {
        A = &field(name);
        N = static_cast<std::size_t>(nx_) * ny_;
    } else if (kind == ArrayKind::Fx) {
        A = &flux(name);
        N = static_cast<std::size_t>(nx_ + 1) * ny_;
    } else { // ArrayKind::Fy
        A = &flux(name);
        N = static_cast<std::size_t>(nx_) * (ny_ + 1);
    }

    for (std::size_t i = 0; i < N; ++i)
        (*A)[i] = ((*A)[i] - a) / b;
}

double Grid2D::fieldAverageAtX(const std::string& name, double x) const
{
    if (!hasField(name))
        throw std::runtime_error("fieldAverageAtX: missing field " + name);

    const auto& F = field(name);
    if (F.size() != static_cast<std::size_t>(nx_ * ny_))
        throw std::runtime_error("fieldAverageAtX: only works with cell-centered fields");

    // map x to column index
    int i = static_cast<int>(std::floor(x / dx_));
    if (i < 0) i = 0;
    if (i > nx_-1) i = nx_-1;

    double sum = 0.0;
    for (int j=0; j<ny_; ++j)
        sum += F[static_cast<std::size_t>(j)*nx_ + i];

    return sum / static_cast<double>(ny_);
}

void Grid2D::transportStepUpwind(const std::string& name, double dt)
{
    if (dt <= 0.0) throw std::runtime_error("transportStepUpwind: dt must be > 0");

    // Ensure fluxes exist
    if (!hasFlux("qx") || !hasFlux("qy"))
        throw std::runtime_error("transportStepUpwind: fluxes qx/qy are missing; run Darcy first");

    const int NX = nx_, NY = ny_;
    const double DX = dx_, DY = dy_;
    const double cellArea = DX * DY; // unit thickness

    // Access/create concentration field
    auto& C = field(name);
    if (C.size() != static_cast<std::size_t>(NX) * NY)
        C.assign(static_cast<std::size_t>(NX) * NY, 0.0);

    // Read face fluxes (Darcy volumetric flux per unit thickness)
    const auto& QX = flux("qx"); // (NX+1)*NY
    const auto& QY = flux("qy"); // NX*(NY+1)

    // Helpers to index cell and faces
    auto idx_cell = [&](int i,int j)->std::size_t {
        return static_cast<std::size_t>(j)*NX + static_cast<std::size_t>(i);
    };
    auto idx_fx = [&](int i,int j)->std::size_t { // vertical faces: i=0..NX, j=0..NY-1
        return static_cast<std::size_t>(j)*(NX+1) + static_cast<std::size_t>(i);
    };
    auto idx_fy = [&](int i,int j)->std::size_t { // horizontal faces: i=0..NX-1, j=0..NY
        return static_cast<std::size_t>(j)*NX + static_cast<std::size_t>(i);
    };

    // Optional: quick CFL guard (upwind explicit)
    {
        double umax = 0.0, vmax = 0.0;
        // max |u_x|
        for (int j=0; j<NY; ++j) {
            for (int i=0; i<=NX; ++i) umax = std::max(umax, std::abs(QX[idx_fx(i,j)]));
        }
        // max |u_y|
        for (int j=0; j<=NY; ++j) {
            for (int i=0; i<NX; ++i)  vmax = std::max(vmax, std::abs(QY[idx_fy(i,j)]));
        }
        double dtCFL = std::numeric_limits<double>::infinity();
        if (umax > 0.0) dtCFL = std::min(dtCFL, DX / umax);
        if (vmax > 0.0) dtCFL = std::min(dtCFL, DY / vmax);
        if (dt > 0.99 * dtCFL && std::isfinite(dtCFL)) {
            // not throwing—just a gentle guardrail note
            fprintf(stderr, "Warning: dt exceeds CFL (dt=%.3e, dtCFL=%.3e)\n", dt, dtCFL);
        }
    }

    // One step: C^{n+1} = C^n - dt/Vol * (F_E - F_W + F_N - F_S)
    // with upwinded face concentrations.
    std::vector<double> Cnew(C.size(), 0.0);

    for (int j=0; j<NY; ++j) {
        for (int i=0; i<NX; ++i) {

            const std::size_t s = idx_cell(i,j);
            const double Cij = C[s];

            // West face (i-1/2): flux qx(i, j) defined at vertical face i
            const double qW = QX[idx_fx(i, j)]; // positive: left->right
            double CW; // upwind concentration at west face into cell (i,j)
            if (qW > 0.0) {
                // flow into cell from west boundary/neighbor: take upstream from left side
                if (i == 0) {
                    // Left boundary: Dirichlet C=1 applied on inflow
                    CW = 1.0;
                } else {
                    CW = C[idx_cell(i-1, j)];
                }
            } else {
                // flow leaving cell (to the west): upwind is the cell itself
                CW = Cij;
            }

            // East face (i+1/2): flux qx(i+1, j) at vertical face i+1
            const double qE = QX[idx_fx(i+1, j)]; // positive: left->right
            double CE;
            if (qE < 0.0) {
                // flow into cell from east side (right->left)
                if (i == NX-1) {
                    // Right boundary: "no-diffusion" outflow; for backflow we choose Cin = 0.
                    // (You can parameterize this if you prefer another inflow value.)
                    CE = 0.0;
                } else {
                    CE = C[idx_cell(i+1, j)];
                }
            } else {
                // flow leaving cell to the east: upwind is the cell
                CE = Cij;
            }

            // South face (i, j-1/2): qy(i, j)
            const double qS = QY[idx_fy(i, j)]; // positive: bottom->top
            double CS;
            if (qS > 0.0) {
                // into cell from south
                if (j == 0) {
                    // bottom boundary: Darcy no-flux ⇒ qS should be 0, but guard anyway
                    CS = Cij; // zero-gradient fallback (no inflow value specified)
                } else {
                    CS = C[idx_cell(i, j-1)];
                }
            } else {
                // leaving cell downward
                CS = Cij;
            }

            // North face (i, j+1/2): qy(i, j+1)
            const double qN = QY[idx_fy(i, j+1)]; // positive: bottom->top
            double CN;
            if (qN < 0.0) {
                // into cell from north
                if (j == NY-1) {
                    // top boundary: Darcy no-flux ⇒ qN should be 0, guard:
                    CN = Cij; // zero-gradient fallback
                } else {
                    CN = C[idx_cell(i, j+1)];
                }
            } else {
                // leaving cell upward
                CN = Cij;
            }

            // Face areas (unit thickness): vertical faces area = DY, horizontal faces area = DX
            const double FW = qW * CW * DY;
            const double FE = qE * CE * DY;
            const double FS = qS * CS * DX;
            const double FN = qN * CN * DX;

            const double net = (FE - FW) + (FN - FS); // net outflow
            Cnew[s] = Cij - (dt / cellArea) * net;
        }
    }

    // Overwrite the field with the new step
    C.swap(Cnew);
}

