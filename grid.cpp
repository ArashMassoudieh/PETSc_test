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
