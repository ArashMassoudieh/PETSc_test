#include "grid.h"
#include <algorithm> // std::transform
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <cmath>
#include <random>
#include <limits>
#include "Matrix_arma.h"
#include "Vector_arma.h"
#include <fstream>
#include <iomanip>
#include "petscmatrix.h"
#include "petscvector.h"
#include "TimeSeries.h"

#ifdef GRID_USE_VTK
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkXMLImageDataWriter.h>
#endif

// If you want PETSc interop, include your wrappers here:
#include "petscvector.h" // adjust include path to your project

// Add to grid.cpp constructor (Grid2D::Grid2D):

Grid2D::Grid2D(int nx, int ny, double Lx, double Ly)
    : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly),
    dx_(Lx/nx), dy_(Ly/ny),
    diffusion_coeff_(0.0),
    porosity_(1.0),
    c_left_(0.0),
    lc_(1.0),
    lambda_x_(0.1),
    lambda_y_(0.1),
    pdf_field_name_("c_u"),
    A(nullptr)
{
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
                                  char delimiter,
                                  int precision,
                                  bool include_header,
                                  bool flipY) const
{
    // --- Check for cell field ---
    auto itF = fields_.find(name);
    if (itF != fields_.end()) {
        const auto& f = itF->second;
        if (f.size() != static_cast<std::size_t>(nx_ * ny_))
            throw std::runtime_error("field '" + name + "' size mismatch (expected nx*ny)");
        writeMatrixRaw(f, ny_, nx_, name, filename, delimiter, precision,
                       include_header, flipY, Grid2D::ArrayKind::Cell);
        return;
    }

    // --- Check for flux ---
    auto itX = fluxes_.find(name);
    if (itX != fluxes_.end()) {
        const auto& f = itX->second;
        const std::size_t sz = f.size();

        if (sz == FxSize()) {
            // Fx = ny rows, nx+1 cols
            if (sz != static_cast<std::size_t>((nx_+1)*ny_))
                throw std::runtime_error("flux '" + name + "' Fx size mismatch");
            writeMatrixRaw(f, ny_, nx_+1, name, filename, delimiter, precision,
                           include_header, flipY, Grid2D::ArrayKind::Fx);
            return;
        }
        if (sz == FySize()) {
            // Fy = ny+1 rows, nx cols
            if (sz != static_cast<std::size_t>(nx_*(ny_+1)))
                throw std::runtime_error("flux '" + name + "' Fy size mismatch");
            writeMatrixRaw(f, ny_+1, nx_, name, filename, delimiter, precision,
                           include_header, flipY, Grid2D::ArrayKind::Fy);
            return;
        }

        throw std::runtime_error("flux '" + name + "' unexpected size = " +
                                 std::to_string(sz));
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
                            bool flipY,
                            ArrayKind kind) const
{
    if (rows <= 0 || cols <= 0)
        throw std::runtime_error("writeMatrixRaw: nonpositive rows/cols");
    if (static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols) != data.size())
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
        switch (kind) {
        case ArrayKind::Cell: return static_cast<std::size_t>(r)*nx_ + c;
        case ArrayKind::Fx:   return static_cast<std::size_t>(r)*(nx_+1) + c;
        case ArrayKind::Fy:   return static_cast<std::size_t>(r)*nx_ + c;
        default: throw std::runtime_error("Unknown ArrayKind in writeMatrixRaw");
        }
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

    auto col_has_known = [&]()->bool {
        for (int jj = 0; jj < ny_; ++jj) {
            if (known[I(it, jj)]) return true;
        }
        return false;
    }();

    auto has_same_row_in = [&](const std::vector<std::size_t>& v)->bool{
        for (auto s : v) { int ii = s % nx_, jj = s / nx_; if (jj == jt) return true; }
        return false;
    };

    auto has_same_col_in = [&](const std::vector<std::size_t>& v)->bool{
        for (auto s : v) { int ii = s % nx_, jj = s / nx_; if (ii == it) return true; }
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
            bool ok_row = !ensure_same_row || !row_has_known || has_same_row_in(cand);
            bool ok_col = !ensure_same_row || !col_has_known || has_same_col_in(cand);
            if (ok_row && ok_col)
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

    // --- vertical faces (QX) ---
    for (int j=0; j<NY; ++j) {
        // Left boundary face (Dirichlet head at west boundary)
        {
            const int i = 0;
            const double kface = kx_at(i,j);
            const double grad  = (head_at(i,j) - H_west) / (0.5*DX);
            QX[FxIndex(0, j)] = -kface * grad;
        }

        // Interior vertical faces
        for (int i=1; i<NX; ++i) {  // strictly < NX
            const double kface = harm(kx_at(i-1,j), kx_at(i,j));
            const double grad  = (head_at(i,j) - head_at(i-1,j)) / DX;
            QX[FxIndex(i, j)] = -kface * grad;
        }

        // Right boundary face (Dirichlet head at east boundary)
        {
            const int i = NX-1;
            const double kface = kx_at(i,j);
            const double grad  = (H_east - head_at(i,j)) / (0.5*DX);
            QX[FxIndex(NX, j)] = -kface * grad;
        }
    }

    // --- horizontal faces (QY) ---
    // --- horizontal faces (QY) ---
    for (int i=0; i<NX; ++i) {
        // Interior horizontal faces
        for (int j=1; j<NY; ++j) {
            const double kface = harm(ky_at(i,j-1), ky_at(i,j));
            const double grad  = (head_at(i,j) - head_at(i,j-1)) / DY;
            QY[FyIndex(i, j)] = -kface * grad;
        }
    }

    // Enforce no-flux at boundaries AFTER filling everything
    for (int i=0; i<NX; ++i) {
        QY[FyIndex(i, 0)]  = 0.0;   // bottom boundary
        QY[FyIndex(i, NY)] = 0.0;   // top boundary
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

TimeSeries<double> Grid2D::toSeries(const std::string& inputFieldName) const {
    auto it = fields_.find(inputFieldName);
    if (it == fields_.end()) {
        throw std::runtime_error("Input field not found: " + inputFieldName);
    }

    const std::vector<double>& input = it->second;
    TimeSeries<double> output;


    int counter = 0;
    for (double x : input) {
        output.append(counter, x);
    }

    return output;

}

void Grid2D::assignFromTimeSeries(
    const TimeSeries<double>& ts,
    const std::string& fieldName,
    ArrayKind kind
    )
{
    std::size_t expectedSize = 0;
    switch (kind) {
    case ArrayKind::Cell:
        expectedSize = static_cast<std::size_t>(nx_) * ny_;
        break;
    case ArrayKind::Fx:
        expectedSize = static_cast<std::size_t>(nx_ + 1) * ny_;
        break;
    case ArrayKind::Fy:
        expectedSize = static_cast<std::size_t>(nx_) * (ny_ + 1);
        break;
    default:
        throw std::runtime_error("assignFromTimeSeries: invalid ArrayKind");
    }

    if (ts.size() != expectedSize) {
        throw std::runtime_error(
            "assignFromTimeSeries: size mismatch (ts.size=" +
            std::to_string(ts.size()) + ", expected=" +
            std::to_string(expectedSize) + ")"
            );
    }

    std::vector<double>* target = nullptr;

    if (kind == ArrayKind::Cell) {
        auto& f = fields_[fieldName];
        f.resize(expectedSize);
        target = &f;
    } else {
        auto& f = flux(fieldName);
        f.resize(expectedSize);
        target = &f;
    }

    for (std::size_t i = 0; i < expectedSize; ++i) {
        (*target)[i] = ts[i].c;   // take the “c” values
    }
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

void Grid2D::normalizeField(const std::string& name, ArrayKind kind)
{
    double a = fieldMean(name, kind);
    double b = fieldStdDev(name,kind);

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

// Utility functions for transport matrix assembly
double Grid2D::getVelocityX(int i, int j, double porosity) const {
    // u at vertical face (i,j): between cells (i-1,j) and (i,j)
    const auto& qx = flux("qx");
    const std::size_t idx = static_cast<std::size_t>(j) * (nx_ + 1) + i;
    return qx[idx] / porosity;
}

double Grid2D::getVelocityY(int i, int j, double porosity) const {
    // v at horizontal face (i,j): between cells (i,j-1) and (i,j)
    const auto& qy = flux("qy");
    const std::size_t idx = static_cast<std::size_t>(j) * nx_ + i;
    return qy[idx] / porosity;
}

PetscInt Grid2D::cellToPetscIndex(int i, int j) const {
    return static_cast<PetscInt>(j * nx_ + i);
}


// -----------------------------
// X-direction contributions
// -----------------------------
void Grid2D::addTransportXTerms(int i, int j, double dt, double D, double porosity,
                                double& diag,
                                std::vector<PetscInt>& cols,
                                std::vector<PetscScalar>& vals) const
{
    const int NX = nx();
    const double DX = dx();
    const double inv_dx  = 1.0 / DX;
    const double inv_dx2 = 1.0 / (DX*DX);

    // face velocities
    double u_west = (i > 0)      ? getVelocityX(i, j, porosity)   : 0.0;
    double u_east = (i < NX - 1) ? getVelocityX(i+1, j, porosity) : 0.0;

    // ---- Interior ----
    if (i > 0 && i < NX-1) {
        double aW = std::max(u_west,0.0) * inv_dx + D*inv_dx2;
        double aE = -std::min(u_east,0.0) * inv_dx + D*inv_dx2;

        diag += aW + aE;

        cols.push_back(cellToPetscIndex(i-1,j));
        vals.push_back(-aW);

        cols.push_back(cellToPetscIndex(i+1,j));
        vals.push_back(-aE);
    }
    // --------------------
    // West boundary (i=0)
    // --------------------
    else if (i == 0) {
        // Diagonal: east outflow + contribution from ghost (-c0 term) + diffusion
        diag += std::max(u_east, 0.0) * inv_dx    // outflow east
                - std::min(u_west, 0.0) * inv_dx    // inflow from west ghost
                + 3.0 * D * inv_dx2;                // diffusion: (2cL - c0 - 2c0 + c1)

        // East neighbor (negative)
        if (NX > 1) {
            double coeff_east = - ( -std::min(u_east, 0.0) * inv_dx + D * inv_dx2 );
            cols.push_back(cellToPetscIndex(i+1, j));
            vals.push_back(coeff_east);
        }

        // NOTE: the constant source term from ghost (2*u_west/dx * c_left + 2*D/dx^2 * c_left)
        // is NOT added here — it must be added in assembleTransportRHS.
    }
    // ---- East boundary ----
    else if (i == NX-1) {
        double aW = std::max(u_west,0.0) * inv_dx + D*inv_dx2;
        diag += aW;
        cols.push_back(cellToPetscIndex(i-1,j));
        vals.push_back(-aW);
        // no east neighbor (outflow)
    }
}

// -----------------------------
// Y-direction contributions
// -----------------------------
void Grid2D::addTransportYTerms(int i, int j, double dt, double D, double porosity,
                                double& diag,
                                std::vector<PetscInt>& cols,
                                std::vector<PetscScalar>& vals) const
{
    const int NY = ny();
    const double DY = dy();
    const double inv_dy  = 1.0 / DY;
    const double inv_dy2 = 1.0 / (DY*DY);

    // Check if variable diffusion field exists
    bool has_D_y = hasFlux("D_y");

    // Get diffusion coefficients at faces
    double D_south = D;  // default
    double D_north = D;  // default

    if (has_D_y) {
        const auto& D_y_field = flux("D_y");
        // D_y is stored at vertical faces: index = j * nx + i
        if (j > 0) {
            D_south = D_y_field[j * nx() + i];
        }
        if (j < NY) {
            D_north = D_y_field[(j + 1) * nx() + i];
        }
    }

    // face velocities
    double v_south = (j > 0)     ? getVelocityY(i, j, porosity)   : 0.0;
    double v_north = (j < NY-1)  ? getVelocityY(i, j+1, porosity) : 0.0;

    // ---- Interior ----
    if (j > 0 && j < NY-1) {
        double aS = std::max(v_south, 0.0) * inv_dy + D_south * inv_dy2;
        double aN = -std::min(v_north, 0.0) * inv_dy + D_north * inv_dy2;
        diag += aS + aN;
        cols.push_back(cellToPetscIndex(i, j-1));
        vals.push_back(-aS);
        cols.push_back(cellToPetscIndex(i, j+1));
        vals.push_back(-aN);
    }
    // ---- Bottom boundary ----
    else if (j == 0) {
        double aN = -std::min(v_north, 0.0) * inv_dy + D_north * inv_dy2;
        diag += aN;
        if (NY > 1) {
            cols.push_back(cellToPetscIndex(i, j+1));
            vals.push_back(-aN);
        }
        // no south face (Neumann)
    }
    // ---- Top boundary ----
    else if (j == NY-1) {
        double aS = std::max(v_south, 0.0) * inv_dy + D_south * inv_dy2;
        diag += aS;
        cols.push_back(cellToPetscIndex(i, j-1));
        vals.push_back(-aS);
        // no north face (Neumann)
    }
}

void Grid2D::assembleTransportMatrix(PETScMatrix *A, double dt) const
{
    if (!hasFlux("qx") || !hasFlux("qy"))
        throw std::runtime_error("assembleTransportMatrix: missing velocity fields qx, qy");

    const int NX = nx(), NY = ny();
    const double inv_dt = 1.0 / dt;

    PetscInt Istart, Iend;
    A->ownershipRange(Istart, Iend);

    for (PetscInt Irow = Istart; Irow < Iend; ++Irow) {
        const int j = static_cast<int>(Irow / NX);
        const int i = static_cast<int>(Irow - j * NX);

        double diag = inv_dt;
        std::vector<PetscInt> cols;
        std::vector<PetscScalar> vals;

        addTransportXTerms(i, j, dt, diffusion_coeff_, porosity_, diag, cols, vals);
        addTransportYTerms(i, j, dt, diffusion_coeff_, porosity_, diag, cols, vals);

        if (!cols.empty()) {
            A->setValues(1, &Irow,
                         static_cast<PetscInt>(cols.size()),
                         cols.data(), vals.data(), ADD_VALUES);
        }
        A->setValue(Irow, Irow, diag, ADD_VALUES);

        /* -------- DEBUG PRINT --------
        if (NX <= 6 && NY <= 6) { // only for small grids
            double sum = diag;
            std::ostringstream oss;
            oss << "Row(" << i << "," << j << "): diag=" << diag;
            for (std::size_t k = 0; k < cols.size(); ++k) {
                sum += vals[k];
                oss << " , N" << k << "=" << vals[k] << "(col " << cols[k] << ")";
            }
            oss << " , RowSum=" << sum << " , ShouldBe=" << inv_dt;
            std::cout << oss.str() << std::endl;
        }
        /* -----------------------------*/
    }

    A->assemble();
}

void Grid2D::assembleTransportRHS(PETScVector& b,
                                  const std::string& c_field,
                                  double dt) const
{
    if (!hasField(c_field))
        throw std::runtime_error("assembleTransportRHS: missing concentration field " + c_field);

    const auto& C = field(c_field);
    const int NX = nx(), NY = ny();
    const double inv_dt = 1.0 / dt;
    const double DX = dx();
    const double inv_dx  = 1.0 / DX;
    const double inv_dx2 = 1.0 / (DX * DX);

    PetscInt Istart, Iend;
    b.ownershipRange(Istart, Iend);

    for (PetscInt Irow = Istart; Irow < Iend; ++Irow) {
        const int j = static_cast<int>(Irow / NX);
        const int i = static_cast<int>(Irow - j * NX);
        const std::size_t cell_idx = static_cast<std::size_t>(j) * NX + i;

        double rhs = C[cell_idx] * inv_dt;

        // --- West boundary: inject boundary concentration ---
        if (i == 0) {
            double u_west = getVelocityX(0, j, porosity_);
            // Advection inflow contribution
            rhs += (std::max(u_west,0.0) * inv_dx) * c_left_;
            // Diffusion ghost contribution
            if (diffusion_coeff_ > 0.0) {
                rhs += (2.0 * diffusion_coeff_ * inv_dx2) * c_left_;
            }
        }

        b.setValue(Irow, rhs, INSERT_VALUES);
    }

    b.assemble();
}



void Grid2D::transportStep(double dt, const char* ksp_prefix)
{
    if (!hasFlux("qx") || !hasFlux("qy"))
        throw std::runtime_error("transportStep: must solve flow first (missing qx, qy fluxes)");

    const int N = nx() * ny();

    // Initialize concentration field if it doesn't exist
    if (!hasField("C")) {
        auto& C = field("C");
        C.assign(static_cast<std::size_t>(N), 0.0);  // Initialize with zeros
    }

    // Add this in transportStep before matrix assembly:
    double u_boundary = getVelocityX(0, 0, porosity_);
    //std::cout << "Boundary velocity u[0,0] = " << u_boundary << std::endl;

    // Assemble RHS vector b
    PETScVector b(N);
    assembleTransportRHS(b, "C", dt);

    //b.saveToFile("RHS.txt");
    // Solve A * C^{n+1} = b
    PETScVector c_new = ksp_prefix ? A->solveNew(b, ksp_prefix) : (A->operator/(b));

    // Extract solution back to field "C" - get ALL values
    auto& C = field("C");
    PetscInt n;
    PetscCallAbort(PETSC_COMM_WORLD, VecGetSize(c_new.raw(), &n));

    std::vector<PetscInt> indices(n);
    std::vector<PetscScalar> values(n);
    for (PetscInt i = 0; i < n; ++i) indices[i] = i;

    PetscCallAbort(PETSC_COMM_WORLD, VecGetValues(c_new.raw(), n, indices.data(), values.data()));

    for (PetscInt i = 0; i < n; ++i) {
        C[static_cast<std::size_t>(i)] = static_cast<double>(values[i]);
    }
    //writeNamedMatrixAuto("C", "C.txt");
}

void Grid2D::SetVal(const std::string& prop, const double& value) {
    if (prop == "diffusion" || prop == "D") {
        if (value < 0.0) throw std::runtime_error("SetVal: diffusion coefficient must be non-negative");
        diffusion_coeff_ = value;
    }
    else if (prop == "porosity") {
        if (value <= 0.0 || value > 1.0)
            throw std::runtime_error("SetVal: porosity must be in range (0,1]");
        porosity_ = value;
    }
    else if (prop == "c_left" || prop == "left_bc") {
        c_left_ = value;
    }
    else {
        throw std::runtime_error("SetVal: unknown property '" + prop + "'. "
                                                                       "Valid properties: 'diffusion'/'D', 'porosity', 'c_left'/'left_bc'");
    }
}

void Grid2D::SolveTransport(const double& t_end,
                            const double& dt,
                            const char* ksp_prefix,
                            int output_interval,
                            const std::string& output_dir,
                            const std::string& filename,
                            TimeSeriesSet<double>* btc_data,   // <-- existing
                            int realization /* = -1 */)        // <-- NEW (optional)
{
    if (t_end <= 0.0) throw std::runtime_error("SolveTransport: t_end must be positive");
    if (dt <= 0.0) throw std::runtime_error("SolveTransport: dt must be positive");
    if (dt > t_end) throw std::runtime_error("SolveTransport: dt cannot be larger than t_end");
    if (!hasFlux("qx") || !hasFlux("qy"))
        throw std::runtime_error("SolveTransport: must call DarcySolve() first");

    const bool use_realization = (realization >= 0);

    // Helper: format realization tag (e.g., r0003)
    auto realizationTag = [&]() -> std::string {
        if (!use_realization) return "";
        std::ostringstream os;
        os << "r" << std::setw(4) << std::setfill('0') << realization;
        return os.str();
    };

    // Helper: build VTI filename for a given step (keeps old naming when realization is off)
    auto makeVtiName = [&](const std::string& base, int step) -> std::string {
        std::ostringstream os;

        if (!use_realization) {
            // OLD behavior: transport_C0000.vti
            os << base << std::setw(4) << std::setfill('0') << step << ".vti";
        } else {
            // NEW behavior: transport_C_r0003_0000.vti
            os << base << "_" << realizationTag() << "_"
               << std::setw(4) << std::setfill('0') << step << ".vti";
        }
        return os.str();
    };

    std::cout << "Transport parameters:\n"
              << "  c_left = " << c_left_ << "\n"
              << "  diffusion = " << diffusion_coeff_ << "\n"
              << "  porosity = " << porosity_ << "\n\n";

    // Initialize breakthrough curve recording if requested
    bool record_btc = (btc_data != nullptr && !BTCLocations_.empty());
    if (record_btc) {
        btc_data->resize(BTCLocations_.size());

        // Set names for each series based on x-location (+ optional realization tag)
        for (size_t i = 0; i < BTCLocations_.size(); i++) {
            std::ostringstream name;
            if (use_realization) {
                name << realizationTag() << "_";
            }
            name << "x=" << std::fixed << std::setprecision(2) << BTCLocations_[i];
            btc_data->setSeriesName(i, name.str());
        }

        std::cout << "Recording breakthrough curves at " << BTCLocations_.size()
                  << " locations\n";
    }

    double current_time = 0.0;
    int step_count = 0;
    const int total_steps = static_cast<int>(std::ceil(t_end / dt));

    std::cout << "Starting transport simulation:\n"
              << "  End time: " << t_end << "\n"
              << "  Time step: " << dt << "\n"
              << "  Total steps: " << total_steps << "\n"
              << "  Diffusion coeff: " << diffusion_coeff_ << "\n"
              << "  Porosity: " << porosity_ << "\n"
              << "  Left BC: " << c_left_ << "\n\n";

    // Assemble transport matrix A
    const int N = nx() * ny();
    A = new PETScMatrix(N, N, 5);
    assembleTransportMatrix(A, dt);

    std::string filename_;
    if (filename == "")
        filename_ = "transport_C";
    else
        filename_ = filename;

    // Write initial state
    writeNamedVTI_Auto("C", output_dir + "/" + makeVtiName(filename_, 0));

    // Record initial BTC values (t=0)
    if (record_btc) {
        std::vector<double> btc_values;
        btc_values.reserve(BTCLocations_.size());

        for (size_t i = 0; i < BTCLocations_.size(); i++) {
            double avg = getAverageAlongY("C", BTCLocations_[i]);
            btc_values.push_back(avg);
        }
        btc_data->append(current_time, btc_values);
    }

    while (current_time < t_end) {
        double this_dt = dt;

        try {
            transportStep(this_dt, ksp_prefix);
            step_count++;
            current_time += this_dt;

            // Record breakthrough curves at this timestep
            if (record_btc) {
                std::vector<double> btc_values;
                btc_values.reserve(BTCLocations_.size());

                for (size_t i = 0; i < BTCLocations_.size(); i++) {
                    double avg = getAverageAlongY("C", BTCLocations_[i]);
                    btc_values.push_back(avg);
                }
                btc_data->append(current_time, btc_values);
            }

            const int progress_interval = std::max(1, std::min(100, total_steps / 10));
            if (step_count % progress_interval == 0 || current_time >= t_end) {
                double progress = (current_time / t_end) * 100.0;
                std::cout << "Step " << step_count << "/" << total_steps
                          << ", Time: " << std::fixed << std::setprecision(6) << current_time
                          << " (" << std::setprecision(1) << progress << "%)\n";
                std::cout << "   Samples: ";
                printSampleC({{0,50}, {10,50}, {50,50}, {100,50}, {200,50}, {299,50}});
            }

            // Write snapshots
            if (output_interval > 0 && step_count % output_interval == 0) {
                writeNamedVTI_Auto("C", output_dir + "/" + makeVtiName(filename_, step_count));
            }
        } catch (const std::exception& e) {
            std::cerr << "Error in transport step " << step_count
                      << " at time " << current_time << ": " << e.what() << std::endl;
            throw;
        }
    }

    std::cout << "\nTransport simulation completed successfully!\n"
              << "Final time: " << current_time << "\n"
              << "Total steps taken: " << step_count << "\n";

    if (record_btc) {
        std::cout << "Recorded " << (*btc_data)[0].size()
                  << " time points for breakthrough curves\n";
    }

    if (hasField("C")) {
        const double mean_c = fieldMean("C", ArrayKind::Cell);
        const double std_c  = fieldStdDev("C", ArrayKind::Cell);
        std::cout << "Final concentration - Mean: " << std::scientific
                  << std::setprecision(4) << mean_c
                  << ", Std Dev: " << std_c << "\n";
    }
}


void Grid2D::printSampleC(const std::vector<std::pair<int,int>>& pts) const {
    if (!hasField("C")) return;
    const auto& C = field("C");
    for (auto [i,j] : pts) {
        if (i>=0 && i<nx_ && j>=0 && j<ny_) {
            std::size_t idx = static_cast<std::size_t>(j)*nx_ + i;
            std::cout << "C(" << i << "," << j << ")=" << C[idx] << "  ";
        }
    }
    std::cout << "\n";
}

std::pair<double,double> Grid2D::fieldMinMax(const std::string& name, ArrayKind kind) const
{
    const std::vector<double>* A = nullptr;
    std::size_t N = 0;

    if (kind == ArrayKind::Cell) {
        auto it = fields_.find(name);
        if (it == fields_.end()) throw std::runtime_error("fieldMinMax: field not found " + name);
        A = &it->second;
        N = static_cast<std::size_t>(nx_) * ny_;
    }
    else if (kind == ArrayKind::Fx) {
        auto it = fluxes_.find(name);
        if (it == fluxes_.end()) throw std::runtime_error("fieldMinMax: flux not found " + name);
        A = &it->second;
        N = static_cast<std::size_t>(nx_ + 1) * ny_;
    }
    else if (kind == ArrayKind::Fy) {
        auto it = fluxes_.find(name);
        if (it == fluxes_.end()) throw std::runtime_error("fieldMinMax: flux not found " + name);
        A = &it->second;
        N = static_cast<std::size_t>(nx_) * (ny_ + 1);
    }
    else {
        throw std::runtime_error("fieldMinMax: unknown ArrayKind");
    }

    if (N == 0) throw std::runtime_error("fieldMinMax: field is empty");

    auto [minIt, maxIt] = std::minmax_element(A->begin(), A->end());
    return { *minIt, *maxIt };
}

void Grid2D::assignConstant(const std::string& name, ArrayKind kind, double value)
{
    std::vector<double>* A = nullptr;
    std::size_t N = 0;

    if (kind == ArrayKind::Cell) {
        auto& f = field(name);  // creates if missing
        f.resize(static_cast<std::size_t>(nx_) * ny_);
        A = &f;
        N = A->size();
    }
    else if (kind == ArrayKind::Fx) {
        auto& f = flux(name);  // creates if missing
        f.resize(static_cast<std::size_t>(nx_ + 1) * ny_);
        A = &f;
        N = A->size();
    }
    else if (kind == ArrayKind::Fy) {
        auto& f = flux(name);
        f.resize(static_cast<std::size_t>(nx_) * (ny_ + 1));
        A = &f;
        N = A->size();
    }
    else {
        throw std::runtime_error("assignConstant: unknown ArrayKind");
    }

    std::fill(A->begin(), A->end(), value);
}

void Grid2D::computeMassBalanceError(const std::string& fieldName) {
    const int NX = nx();
    const int NY = ny();

    // allocate/overwrite field
    auto& MB = field(fieldName);
    MB.assign(static_cast<std::size_t>(NX)*NY, 0.0);

    const auto& QX = flux("qx"); // size: (NX+1)*NY
    const auto& QY = flux("qy"); // size: NX*(NY+1)

    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            // fluxes at the four faces of cell (i,j)
            double q_west  = QX[j*(NX+1) + i];     // inflow (+ into cell)
            double q_east  = QX[j*(NX+1) + (i+1)]; // outflow (+ out of cell)
            double q_south = QY[j*NX + i];         // inflow (+ into cell)
            double q_north = QY[(j+1)*NX + i];     // outflow (+ out of cell)

            // net divergence = inflows - outflows
            double balance = (q_west - q_east) + (q_south - q_north);

            MB[static_cast<std::size_t>(j)*NX + i] = balance;
        }
    }
}

void Grid2D::computeRowSumErrorField(PETScMatrix* A,
                                     double dt,
                                     const std::string& fieldName)
{
    const int NX = nx(), NY = ny();
    const PetscInt N = NX * NY;

    auto& ERR = field(fieldName);
    ERR.assign(static_cast<std::size_t>(N), 0.0);

    PetscInt Istart, Iend;
    A->ownershipRange(Istart, Iend);

    const double inv_dt = 1.0 / dt;

    for (PetscInt Irow = Istart; Irow < Iend; ++Irow) {
        PetscInt ncols;
        const PetscInt *cols;
        const PetscScalar *vals;

        // Grab nonzeros in this row
        PetscErrorCode ierr = MatGetRow(A->raw(), Irow, &ncols, &cols, &vals);
        if (ierr) throw std::runtime_error("MatGetRow failed");

        double sum = 0.0;
        for (PetscInt k = 0; k < ncols; ++k) {
            sum += static_cast<double>(vals[k]);
        }

        // Compare to expected 1/dt
        ERR[static_cast<std::size_t>(Irow)] = sum - inv_dt;

        MatRestoreRow(A->raw(), Irow, &ncols, &cols, &vals);
    }
}

TimeSeries<double> Grid2D::sampleSecondDerivative(
    const std::string& fieldName,
    ArrayKind kind,
    DerivDir dir,
    int nPoints,
    double delta,
    unsigned long seed) const
{
    if (!hasField(fieldName) && !hasFlux(fieldName)) {
        throw std::runtime_error("Field/flux not found: " + fieldName);
    }

    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    std::uniform_real_distribution<double> Ux(0.0, Lx_);
    std::uniform_real_distribution<double> Uy(0.0, Ly_);

    TimeSeries<double> ts;

    for (int k = 0; k < nPoints; ++k) {
        double x = Ux(rng);
        double y = Uy(rng);

        // central location
        double g0 = interpolate(fieldName, kind, x, y, /*clamp=*/true);

        double gplus, gminus;
        if (dir == DerivDir::X) {
            gplus  = interpolate(fieldName, kind, x+delta, y, true);
            gminus = interpolate(fieldName, kind, x-delta, y, true);
        } else {
            gplus  = interpolate(fieldName, kind, x, y+delta, true);
            gminus = interpolate(fieldName, kind, x, y-delta, true);
        }

        double secondDeriv = (gminus - 2.0*g0 + gplus) / (delta*delta);

        ts.append(g0, secondDeriv);
    }

    return ts;
}

TimeSeries<double> Grid2D::exportFieldToTimeSeries(
    const std::string& fieldName,
    ArrayKind kind
    ) const
{
    const std::vector<double>* src = nullptr;

    if (kind == ArrayKind::Cell) {
        src = &field(fieldName);
    }
    else if (kind == ArrayKind::Fx || kind == ArrayKind::Fy) {
        src = &flux(fieldName);
    }
    else {
        throw std::runtime_error("exportFieldToTimeSeries: unknown ArrayKind");
    }

    TimeSeries<double> ts;

    for (std::size_t i = 0; i < src->size(); ++i) {
        ts.append(static_cast<double>(i), (*src)[i]);
    }

    return ts;
}


TimeSeries<double> Grid2D::sampleGaussianPerturbation(
    const std::string& fieldName,
    ArrayKind kind,
    int nSamples,
    double delta,
    unsigned long seed,
    PerturbDir dir  // new argument
    ) const {
    if (!hasField(fieldName) && !hasFlux(fieldName)) {
        throw std::runtime_error("sampleGaussianPerturbation: field/flux not found: " + fieldName);
    }

    // RNGs
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    std::uniform_real_distribution<double> Ux(0.0, Lx_);
    std::uniform_real_distribution<double> Uy(0.0, Ly_);
    std::uniform_real_distribution<double> Utheta(0.0, 2.0 * M_PI);
    std::normal_distribution<double> N01(0.0, 1.0);

    TimeSeries<double> ts;

    for (int k = 0; k < nSamples; ++k) {
        // 1. Random location
        double x0 = Ux(rng);
        double y0 = Uy(rng);

        // 2. Original field value
        double v0 = interpolate(fieldName, kind, x0, y0, /*clamp=*/true);

        // 3. Random Gaussian perturbation
        double eps = N01(rng);

        // 4. Direction handling
        double dx = 0.0, dy = 0.0;
        switch (dir) {
        case PerturbDir::Radial: {
            double theta = Utheta(rng);
            dx = delta * eps * std::cos(theta);
            dy = delta * eps * std::sin(theta);
            break;
        }
        case PerturbDir::XOnly:
            dx = delta * eps;
            dy = 0.0;
            break;
        case PerturbDir::YOnly:
            dx = 0.0;
            dy = delta * eps;
            break;
        }

        // 5. Perturbed coordinates
        double x1 = std::max(0.0, std::min(Lx_, x0 + dx));
        double y1 = std::max(0.0, std::min(Ly_, y0 + dy));

        // 6. Perturbed field value
        double v1 = interpolate(fieldName, kind, x1, y1, /*clamp=*/true);

        // 7. Add to TimeSeries (t = v0, c = v1)
        ts.append(v0, v1);
    }

    return ts;
}

// Add to grid.cpp
std::pair<double, double> Grid2D::getVelocityAt(double x, double y,
                                                const std::string& qx_name,
                                                const std::string& qy_name) const
{
    // Check bounds
    if (x < 0.0 || x > Lx_ || y < 0.0 || y > Ly_) {
        throw std::runtime_error("getVelocityAt: point (" + std::to_string(x) + ", "
                                 + std::to_string(y) + ") is outside domain");
    }

    // Check if flux fields exist
    if (!hasFlux(qx_name)) {
        throw std::runtime_error("getVelocityAt: flux field '" + qx_name + "' not found");
    }
    if (!hasFlux(qy_name)) {
        throw std::runtime_error("getVelocityAt: flux field '" + qy_name + "' not found");
    }

    // Get flux arrays - changed from const double* to const std::vector<double>&
    const std::vector<double>& qx = flux(qx_name);
    const std::vector<double>& qy = flux(qy_name);

    // Find cell indices (which cell contains this point)
    // Cell (i,j) has corners at (i*dx, j*dy) and ((i+1)*dx, (j+1)*dy)
    int i = static_cast<int>(x / dx_);
    int j = static_cast<int>(y / dy_);

    // Clamp to valid cell indices
    i = std::min(i, nx_ - 1);
    j = std::min(j, ny_ - 1);

    // Local coordinates within the cell [0,1] x [0,1]
    double xi = (x - i * dx_) / dx_;   // normalized x position in cell
    double eta = (y - j * dy_) / dy_;  // normalized y position in cell

    // Clamp to [0,1] to handle boundary cases
    xi = std::max(0.0, std::min(1.0, xi));
    eta = std::max(0.0, std::min(1.0, eta));

    // Interpolate qx (defined at x-faces)
    // qx is defined at faces: left face of cell (i,j) and right face
    // For cell (i,j):
    //   - qx at left face (x = i*dx) has index: j*(nx_+1) + i
    //   - qx at right face (x = (i+1)*dx) has index: j*(nx_+1) + (i+1)

    double vx = 0.0;
    if (i < nx_) {
        // Bilinear interpolation for qx
        // We need 4 qx values at the corners of a cell-centered grid
        int idx_left = j * (nx_ + 1) + i;
        int idx_right = j * (nx_ + 1) + (i + 1);

        if (j < ny_ - 1) {
            int idx_left_top = (j + 1) * (nx_ + 1) + i;
            int idx_right_top = (j + 1) * (nx_ + 1) + (i + 1);

            // Bilinear interpolation in eta direction first, then xi
            double qx_bottom = (1.0 - xi) * qx[idx_left] + xi * qx[idx_right];
            double qx_top = (1.0 - xi) * qx[idx_left_top] + xi * qx[idx_right_top];
            vx = (1.0 - eta) * qx_bottom + eta * qx_top;
        } else {
            // On top boundary, only interpolate in xi direction
            vx = (1.0 - xi) * qx[idx_left] + xi * qx[idx_right];
        }
    }

    // Interpolate qy (defined at y-faces)
    // qy is defined at faces: bottom face of cell (i,j) and top face
    // For cell (i,j):
    //   - qy at bottom face (y = j*dy) has index: j*nx_ + i
    //   - qy at top face (y = (j+1)*dy) has index: (j+1)*nx_ + i

    double vy = 0.0;
    if (j < ny_) {
        // Bilinear interpolation for qy
        int idx_bottom = j * nx_ + i;
        int idx_top = (j + 1) * nx_ + i;

        if (i < nx_ - 1) {
            int idx_bottom_right = j * nx_ + (i + 1);
            int idx_top_right = (j + 1) * nx_ + (i + 1);

            // Bilinear interpolation in xi direction first, then eta
            double qy_left = (1.0 - eta) * qy[idx_bottom] + eta * qy[idx_top];
            double qy_right = (1.0 - eta) * qy[idx_bottom_right] + eta * qy[idx_top_right];
            vy = (1.0 - xi) * qy_left + xi * qy_right;
        } else {
            // On right boundary, only interpolate in eta direction
            vy = (1.0 - eta) * qy[idx_bottom] + eta * qy[idx_top];
        }
    }

    return std::make_pair(vx, vy);
}

// Standard normal PDF
double Grid2D::phi(double z)
{
    return gsl_ran_gaussian_pdf(z, 1.0);
}

// Standard normal CDF
double Grid2D::Phi(double z)
{
    return gsl_cdf_gaussian_P(z, 1.0);
}

// Inverse standard normal CDF (probit function)
double Grid2D::Phi_inv(double u)
{
    // Clamp u to valid range (0, 1) to avoid numerical issues
    if (u <= 0.0) u = 1e-10;
    if (u >= 1.0) u = 1.0 - 1e-10;

    return gsl_cdf_gaussian_Pinv(u, 1.0);
}

// Mixing kernel: κ(v) = v/lc + D/λ_x² + D/λ_y²
double Grid2D::kappa(double v, double lc, double lambda_x, double lambda_y) const
{
    double D = diffusion_coeff_;  // Use the existing diffusion coefficient

    double term1 = v / lc;
    double term2 = D / (lambda_x * lambda_x);
    double term3 = D / (lambda_y * lambda_y);

    return term1 + term2 + term3;
}

// Set mixing parameters
void Grid2D::setMixingParams(double lc, double lambda_x, double lambda_y)
{
    if (lc <= 0.0 || lambda_x <= 0.0 || lambda_y <= 0.0) {
        throw std::runtime_error("setMixingParams: all parameters must be positive");
    }

    lc_ = lc;
    lambda_x_ = lambda_x;
    lambda_y_ = lambda_y;

    std::cout << "Mixing parameters set:\n"
              << "  lc = " << lc_ << "\n"
              << "  lambda_x = " << lambda_x_ << "\n"
              << "  lambda_y = " << lambda_y_ << "\n";
}

// Add to grid.cpp:

void Grid2D::convertToUniformScore(const std::string& c_field_name,
                                   const std::string& u_field_name)
{
    const auto& c_field = field(c_field_name);
    auto& u_field = field(u_field_name);

    if (c_field.empty()) {
        throw std::runtime_error("convertToUniformScore: field " + c_field_name + " is empty");
    }

    const int N = numCells();

    // Create sorted indices for ranking
    std::vector<std::pair<double, int>> sorted_vals;
    sorted_vals.reserve(N);

    for (int i = 0; i < N; ++i) {
        sorted_vals.push_back({c_field[i], i});
    }

    // Sort by concentration value
    std::sort(sorted_vals.begin(), sorted_vals.end());

    // Assign uniform scores based on rank
    for (int rank = 0; rank < N; ++rank) {
        int idx = sorted_vals[rank].second;
        // Use (rank + 0.5) / N to avoid exactly 0 or 1
        u_field[idx] = (rank + 0.5) / N;
    }

    std::cout << "Converted " << c_field_name << " to uniform scores in "
              << u_field_name << std::endl;
}

void Grid2D::convertFromUniformScore(const std::string& u_field_name,
                                     const std::string& c_field_name)
{
    const auto& u_field = field(u_field_name);
    auto& c_field = field(c_field_name);

    if (u_field.empty()) {
        throw std::runtime_error("convertFromUniformScore: field " + u_field_name + " is empty");
    }

    const int N = numCells();

    // For now, just convert using normal scores
    // In practice, you'd need to store the original marginal distribution
    for (int i = 0; i < N; ++i) {
        double u = u_field[i];
        // Clamp to valid range
        u = std::max(1e-10, std::min(1.0 - 1e-10, u));

        // Convert to normal score
        c_field[i] = Phi_inv(u);
    }

    std::cout << "Converted uniform scores " << u_field_name
              << " back to " << c_field_name << std::endl;
}

double Grid2D::velocityAtU(double x, double u) const
{
    // This requires interpolation from the flux field
    // For now, we'll assume a simple approach where we find the y-location
    // corresponding to u and then get velocity at (x, y)

    // This is a placeholder - you'll need to implement proper interpolation
    // based on your specific needs

    throw std::runtime_error("velocityAtU: not yet implemented - need to map u to y");
}


void Grid2D::SolveMixingPDF(const double& t_end,
                            const double& dt,
                            const char* ksp_prefix,
                            int output_interval,
                            const std::string& output_dir)
{
    if (t_end <= 0.0) throw std::runtime_error("SolveMixingPDF: t_end must be positive");
    if (dt <= 0.0) throw std::runtime_error("SolveMixingPDF: dt must be positive");
    if (dt > t_end) throw std::runtime_error("SolveMixingPDF: dt cannot be larger than t_end");
    if (!hasFlux("qx") || !hasFlux("qy"))
        throw std::runtime_error("SolveMixingPDF: must call DarcySolve() first");



    std::cout << "Mixing PDF parameters:\n"
              << "  c_left = " << c_left_ << "\n"
              << "  diffusion = " << diffusion_coeff_ << "\n"
              << "  lc = " << lc_ << "\n"
              << "  lambda_x = " << lambda_x_ << "\n"
              << "  lambda_y = " << lambda_y_ << "\n\n";

    // Helper function for output paths
    auto makeOutputPath = [&output_dir](const std::string& filename) {
        if (output_dir.empty()) return filename;
        if (output_dir.back() == '/' || output_dir.back() == '\\') {
            return output_dir + filename;
        }
        return output_dir + "/" + filename;
    };

    double current_time = 0.0;
    int step_count = 0;
    const int total_steps = static_cast<int>(std::ceil(t_end / dt));

    std::cout << "Starting mixing PDF simulation:\n"
              << "  End time: " << t_end << "\n"
              << "  Time step: " << dt << "\n"
              << "  Total steps: " << total_steps << "\n"
              << "  Grid: " << nx_ << " x " << ny_ << " (x is spatial, y is u-space)\n\n";

    // Initialize PDF field if not exists
    pdf_field_name_ = "c_u";
    if (!hasField(pdf_field_name_)) {
        // Initialize with uniform distribution at left boundary
        auto& c_u = field(pdf_field_name_);
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int idx = j * nx_ + i;
                if (i == 0) {
                    // Left boundary: uniform distribution
                    c_u[idx] = 1.0;  // ∫c_u du = 1
                } else {
                    c_u[idx] = 0.0;
                }
            }
        }
    }

    // Assemble mixing PDF matrix
    const int N = nx() * ny();
    A = new PETScMatrix(N, N, 9);  // 9-point stencil for 2D with cross-derivatives
    assembleMixingPDFMatrix(A, dt);

    // Write initial state
    writeNamedVTI_Auto(pdf_field_name_, makeOutputPath("mixing_pdf_step000.vti"));

    while (current_time < t_end) {
        double this_dt = dt;
        try {
            mixingPDFStep(this_dt, ksp_prefix);
            step_count++;
            current_time += this_dt;

            const int progress_interval = std::max(1, std::min(100, total_steps / 10));
            if (step_count % progress_interval == 0 || current_time >= t_end) {
                double progress = (current_time / t_end) * 100.0;
                std::cout << "Step " << step_count << "/" << total_steps
                          << ", Time: " << std::fixed << std::setprecision(6) << current_time
                          << " (" << std::setprecision(1) << progress << "%)\n";
            }

            // Write snapshots
            if (output_interval > 0 && step_count % output_interval == 0) {
                std::ostringstream fname;
                fname << "mixing_pdf_step"
                      << std::setw(4) << std::setfill('0') << step_count << ".vti";
                writeNamedVTI_Auto(pdf_field_name_, makeOutputPath(fname.str()));
            }
        } catch (const std::exception& e) {
            std::cerr << "Error in mixing PDF step " << step_count
                      << " at time " << current_time << ": " << e.what() << std::endl;
            throw;
        }
    }

    std::cout << "\nMixing PDF simulation completed successfully!\n"
              << "Final time: " << current_time << "\n"
              << "Total steps taken: " << step_count << "\n";
}

void Grid2D::assembleMixingPDFMatrix(PETScMatrix* A, double dt)
{
    const auto& qx = flux("qx");
    const auto& qy = flux("qy");
    const double dx = dx_;
    const double du = dy_;
    const double D = diffusion_coeff_;
    const double inv_dx = 1.0 / dx;
    const double inv_dx2 = 1.0 / (dx * dx);

    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = j * nx_ + i;
            double u = static_cast<double>(j) / (ny_ - 1);

            // Boundary conditions at u=0 and u=1
            if (j == 0 || j == ny_ - 1) {
                // Top/bottom in u-space: zero flux dc/du = 0
                A->setValue(idx, idx, 1.0);
                if (j == 0 && j + 1 < ny_) {
                    A->setValue(idx, idx + nx_, -1.0);
                } else if (j == ny_ - 1 && j - 1 >= 0) {
                    A->setValue(idx, idx - nx_, -1.0);
                }
                continue;
            }

            // Right boundary: zero gradient dc/dx = 0
            if (i == nx_ - 1) {
                A->setValue(idx, idx, 1.0);
                A->setValue(idx, idx - 1, -1.0);
                continue;
            }

            // Interior AND left boundary (i=0) - treat the same way!
            int idx_qx_left = j * (nx_ + 1) + i;
            int idx_qx_right = j * (nx_ + 1) + i + 1;
            double vx = 0.5 * (qx[idx_qx_left] + qx[idx_qx_right]);

            // Get velocity at left face (for i=0, this is the ghost boundary)
            double vx_west = (i > 0) ? qx[idx_qx_left] : qx[idx_qx_left];
            double vx_east = qx[idx_qx_right];

            // Calculate φ[Φ^{-1}(u)]
            double z = Phi_inv(u);
            double phi_z = phi(z);
            double phi_z_sq = phi_z * phi_z;

            // Calculate κ(v)
            double kappa_v = kappa(vx, lc_, lambda_x_, lambda_y_);

            // Diffusion coefficient in u-space
            double D_u = phi_z_sq * kappa_v;

            double beta_x = dt * D / (dx * dx);
            double beta_u = dt * D_u / (du * du);

            // Diagonal term
            double diag = 1.0;

            // X-direction advection-diffusion
            if (i == 0) {
                // Left boundary - ghost cell approach
                diag += std::max(vx_east, 0.0) * dt * inv_dx    // outflow east
                        - std::min(vx_west, 0.0) * dt * inv_dx   // inflow from ghost
                        + 2.0 * beta_x;                          // diffusion

                // East neighbor
                double coeff_east = - ( -std::min(vx_east, 0.0) * dt * inv_dx + beta_x );
                A->setValue(idx, idx + 1, -coeff_east);

            } else {
                // Interior: upwind advection
                double alpha_xm, alpha_xp;
                if (vx_west > 0) {
                    alpha_xm = vx_west * dt * inv_dx + beta_x;
                    alpha_xp = -beta_x;
                } else {
                    alpha_xm = beta_x;
                    alpha_xp = -vx_west * dt * inv_dx - beta_x;
                }

                diag += alpha_xm + alpha_xp;
                A->setValue(idx, idx - 1, -alpha_xm);
                A->setValue(idx, idx + 1, -alpha_xp);
            }

            // U-direction diffusion (central differences)
            diag += 2.0 * beta_u;
            A->setValue(idx, idx - nx_, -beta_u);
            A->setValue(idx, idx + nx_, -beta_u);

            if (i == 0 && j >= 1 && j <= 3) {
                double coeff_east = - ( -std::min(vx_east, 0.0) * dt * inv_dx + beta_x );
                std::cout << "Matrix i=0, j=" << j
                          << ", vx_west=" << vx_west
                          << ", vx_east=" << vx_east
                          << ", diag (FINAL)=" << diag
                          << ", coeff_east=" << coeff_east
                          << ", beta_u=" << beta_u
                          << std::endl;
            }

            A->setValue(idx, idx, diag);
        }
    }
    A->assemble();
}

void Grid2D::mixingPDFStep(double dt, const char* ksp_prefix)
{
    auto& c_u = field(pdf_field_name_);
    const int N = numCells();
    const auto& qx = flux("qx");
    const double dx = dx_;
    const double inv_dx = 1.0 / dx;
    const double inv_dx2 = 1.0 / (dx * dx);
    const double D = diffusion_coeff_;


    // Debug: Check current state
    double sum_before = 0.0;
    for (int i = 0; i < N; ++i) {
        sum_before += c_u[i];
    }
    std::cout << "  Before solve: sum(c_u) = " << sum_before << std::endl;

    // Create RHS vector
    PETScVector b(N);

    static int step_count = 0;
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = j * nx_ + i;
            double rhs = c_u[idx];  // Current concentration

            if (i == nx_ - 1) {
                rhs = 0.0;  // RHS for constraint equation
                b.setValue(idx, rhs);
                continue;
            }

            // Add boundary contribution at i=0
            if (i == 0) {
                int idx_qx_left = j * (nx_ + 1);
                double vx_west = qx[idx_qx_left];

                // **ADD DEBUG HERE**
                if (step_count < 1) {
                    double adv_contrib = std::max(vx_west, 0.0) * dt * inv_dx;
                    double diff_contrib = 2.0 * D * dt * inv_dx2;
                    std::cout << "  j=" << j << ", vx_west=" << vx_west
                              << ", c_u[" << idx << "]=" << c_u[idx]  // <-- ADD THIS
                              << ", adv_contrib=" << adv_contrib
                              << ", diff_contrib=" << diff_contrib
                              << ", total_contrib=" << (adv_contrib + diff_contrib) * c_left_
                              << std::endl;
                }

                // Advection inflow from ghost
                rhs += (std::max(vx_west, 0.0) * dt * inv_dx) * c_left_;
                // Diffusion from ghost
                rhs += (2.0 * D * dt * inv_dx2) * c_left_;
            }
            b.setValue(idx, rhs);
        }
    }
    b.assemble();

    // **DEBUG: Print the RHS for the first few steps**

    if (step_count < 3) {
        std::cout << "\n  === RHS vector (step " << step_count << ") ===" << std::endl;
        std::vector<PetscInt> idx(N);
        std::vector<PetscScalar> vals(N);
        for (int i = 0; i < N; ++i) idx[i] = i;
        PetscCallAbort(PETSC_COMM_WORLD, VecGetValues(b.raw(), N, idx.data(), vals.data()));

        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int ii = j * nx_ + i;
                std::cout << std::setw(8) << std::setprecision(3) << vals[ii] << " ";
            }
            std::cout << std::endl;
        }

        // Also print the solution after solving
        PETScVector x(N);
        A->solve(b, x, ksp_prefix);

        std::cout << "\n  === Solution vector (step " << step_count << ") ===" << std::endl;
        std::vector<PetscScalar> vals_x(N);
        PetscCallAbort(PETSC_COMM_WORLD, VecGetValues(x.raw(), N, idx.data(), vals_x.data()));

        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                int ii = j * nx_ + i;
                std::cout << std::setw(8) << std::setprecision(3) << vals_x[ii] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        // Extract solution and update field
        for (int i = 0; i < N; ++i) {
            c_u[i] = vals_x[i];
        }

        step_count++;
    } else {
        // Normal solve without debug output
        PETScVector x(N);
        A->solve(b, x, ksp_prefix);

        // Extract solution and update field
        std::vector<PetscInt> idx(N);
        std::vector<PetscScalar> vals_x(N);
        for (int i = 0; i < N; ++i) idx[i] = i;
        PetscCallAbort(PETSC_COMM_WORLD, VecGetValues(x.raw(), N, idx.data(), vals_x.data()));

        for (int i = 0; i < N; ++i) {
            c_u[i] = vals_x[i];
        }
    }

    // Debug: Check new state and compute change
    double sum_after = 0.0;
    for (int i = 0; i < N; ++i) {
        sum_after += c_u[i];
    }
    double change = sum_after - sum_before;
    std::cout << "  After solve: sum(c_u) = " << sum_after << std::endl;
    std::cout << "  Change: " << change << std::endl;
}


TimeSeries<double> Grid2D::extractFieldCDF(const std::string& field_name,
                                           ArrayKind kind,
                                           int num_bins,
                                           double threshold) const
{
    // Export field to TimeSeries
    TimeSeries<double> field_data = exportFieldToTimeSeries(field_name, kind);

    // Collect values above threshold
    std::vector<double> sorted_values;
    sorted_values.reserve(field_data.size());
    for (size_t i = 0; i < field_data.size(); ++i) {
        double val = field_data[i].c;
        if (val >= threshold) {  // Only include values >= threshold
            sorted_values.push_back(val);
        }
    }

    if (sorted_values.empty()) return TimeSeries<double>();

    std::sort(sorted_values.begin(), sorted_values.end());

    // Create inverse CDF from the filtered data
    TimeSeries<double> inverse_cdf;
    size_t n = sorted_values.size();
    size_t step = std::max(1UL, n / num_bins);

    for (size_t i = 0; i < n; i += step) {
        double u = static_cast<double>(i) / static_cast<double>(n - 1);
        inverse_cdf.addPoint(u, sorted_values[i]);
    }

    // Always include the last point at u=1.0
    inverse_cdf.addPoint(1.0, sorted_values[n - 1]);

    return inverse_cdf;
}

double Grid2D::getFieldValueAtCDF(const std::string& field_name, ArrayKind kind, double u) const
{
    TimeSeries<double> inverse_cdf = extractFieldCDF(field_name, kind);
    return inverse_cdf.interpol(u);
}

void Grid2D::computeMixingDiffusionCoefficient()
{
    // Create or get the D_y flux array
    // D_y is located at vertical faces: (nx) × (ny+1)
    if (!hasFlux("D_y")) {
        std::vector<double> D_y_data(nx_ * (ny_ + 1), 0.0);
        fluxes_["D_y"] = D_y_data;
    }

    auto& D_y = flux("D_y");

    // For each vertical face
    for (int j = 0; j <= ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = j * nx_ + i;

            // Compute u at this face (same as for q_y)
            double u;
            if (j == 0) {
                u = 0.0;
            } else if (j == ny_) {
                u = 1.0;
            } else {
                // Average of cells above and below
                u = static_cast<double>(j - 0.5) / (ny_ - 1);
            }

            // Get velocity at this location (average from neighboring cells)
            double vx;
            if (i == 0) {
                // At left boundary, use velocity from first cell
                int idx_qx = j * (nx_ + 1);
                vx = flux("qx")[idx_qx];
            } else if (i == nx_ - 1) {
                // At right boundary, use velocity from last cell
                int idx_qx = j * (nx_ + 1) + nx_;
                vx = flux("qx")[idx_qx];
            } else {
                // Interior: average from left and right faces
                int idx_qx_left = j * (nx_ + 1) + i;
                int idx_qx_right = j * (nx_ + 1) + i + 1;
                vx = 0.5 * (flux("qx")[idx_qx_left] + flux("qx")[idx_qx_right]);
            }

            // Compute phi(Phi^{-1}(u))
            double z = Phi_inv(u);
            double phi_z = phi(z);
            double phi_z_sq = phi_z * phi_z;

            // Compute kappa(v)
            double kappa_v = kappa(vx, lc_, lambda_x_, lambda_y_);

            // Compute D_y = phi^2(z) * kappa(v)
            D_y[idx] = phi_z_sq * kappa_v;
        }
    }
}

double Grid2D::getFieldValueAt(const std::string& field_name, double x, double y) const
{
    // Check if field exists
    auto it = fields_.find(field_name);
    if (it == fields_.end()) {
        throw std::runtime_error("Field '" + field_name + "' not found");
    }

    // Check bounds
    if (x < 0.0 || x > Lx_ || y < 0.0 || y > Ly_) {
        throw std::out_of_range("Coordinates out of domain range in getFieldValueAt");
    }

    // Find grid indices that bracket (x,y)
    double xi = x / dx_;
    double yj = y / dy_;

    int i0 = static_cast<int>(std::floor(xi));
    int j0 = static_cast<int>(std::floor(yj));
    int i1 = std::min(i0 + 1, nx_ - 1);
    int j1 = std::min(j0 + 1, ny_ - 1);

    // Interpolation weights
    double alpha_x = xi - i0;
    double alpha_y = yj - j0;

    // Bilinear interpolation
    const auto& field = it->second;
    double val00 = field[idx(i0, j0)];
    double val10 = field[idx(i1, j0)];
    double val01 = field[idx(i0, j1)];
    double val11 = field[idx(i1, j1)];

    double val_y0 = (1.0 - alpha_x) * val00 + alpha_x * val10;
    double val_y1 = (1.0 - alpha_x) * val01 + alpha_x * val11;

    return (1.0 - alpha_y) * val_y0 + alpha_y * val_y1;
}

double Grid2D::getAverageAlongY(const std::string& field_name, double x) const
{
    // Check if field exists
    auto it = fields_.find(field_name);
    if (it == fields_.end()) {
        throw std::runtime_error("Field '" + field_name + "' not found");
    }

    // Check x bounds
    if (x < 0.0 || x > Lx_) {
        throw std::out_of_range("x coordinate out of domain range in getAverageAlongY");
    }

    // Sample along y at ny_ points
    double sum = 0.0;

    for (int j = 0; j < ny_; j++) {
        double y = j * dy_;
        sum += getFieldValueAt(field_name, x, y);
    }

    return sum / static_cast<double>(ny_);
}

void Grid2D::setBTCLocations(const std::vector<double>& locations)
{
    BTCLocations_ = locations;
}

const std::vector<double>& Grid2D::getBTCLocations() const
{
    return BTCLocations_;
}

std::vector<double>& Grid2D::getBTCLocations()
{
    return BTCLocations_;
}
