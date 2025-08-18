#pragma once
#include <string>
#include <vector>
#include <functional>
#include <unordered_map>
#include <cassert>
#include <cstddef>

/**
 * @file grid.h
 * @brief Lightweight 2D structured, cell-centered grid with named fields and face flux storage.
 *
 * Design goals:
 *  - General-purpose container for PDEs (Darcy head/pressure/concentration, etc.)
 *  - Cell-centered scalar fields, registered by name (e.g., "h", "p", "c").
 *  - Face-centered flux arrays:
 *      Fx at vertical faces between (i-1,j) and (i,j), size (nx-1)*ny
 *      Fy at horizontal faces between (i,j-1) and (i,j), size nx*(ny-1)
 *  - Utilities:
 *      * index helpers (cell <-> linear; face indices)
 *      * geometry access (dx, dy, coordinates)
 *      * iteration over cells/interior/boundary
 *      * harmonic means for face properties (e.g., conductivity)
 *      * pack/unpack a cell field to a PETScVector (optional hook)
 *
 * This grid does NOT depend on PETSc; to interoperate, call pack/unpack helpers
 * you’ll find at the bottom (declared but if you don’t include petsc headers,
 * they will be compiled only in grid.cpp where PETSc headers are present).
 */

class PETScVector; // fwd-decl to match your wrapper name via include order

class Grid2D {
public:
    // --------- Construction & geometry ---------

    /**
   * @brief Create a cell-centered grid on [0,Lx] × [0,Ly] with nx × ny cells.
   * The cell (i,j) (0-based) centers are at (i*dx, j*dy), where:
   *  dx = Lx / (nx - 1), dy = Ly / (ny - 1).
   */
    Grid2D(int nx, int ny, double Lx, double Ly);

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    double Lx() const { return Lx_; }
    double Ly() const { return Ly_; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }

    /// total number of cells
    std::size_t numCells() const { return static_cast<std::size_t>(nx_)*ny_; }

    // --------- Indexing helpers ---------

    /// Flattened index for cell (i,j) in row-major order: I = j*nx + i
    std::size_t cellIndex(int i, int j) const {
        assert(i>=0 && i<nx_ && j>=0 && j<ny_);
        return static_cast<std::size_t>(j)*nx_ + i;
    }

    // --- Index helper ---
    inline int idx(int i, int j) const { return j*nx_ + i; }

    /// (i,j) from flattened index I
    void cellIJ(std::size_t I, int& i, int& j) const {
        i = static_cast<int>(I % nx_);
        j = static_cast<int>(I / nx_);
    }

    /// X,Y coordinates of cell center (i,j)
    void cellCenter(int i, int j, double& x, double& y) const {
        x = i * dx_; y = j * dy_;
    }

    /// true if (i,j) lies on the domain boundary
    bool isBoundary(int i, int j) const {
        return (i==0 || j==0 || i==nx_-1 || j==ny_-1);
    }

    // --------- Face indexing (for flux arrays) ---------
    //
    // Fx is stored on vertical faces between (i-1,j) and (i,j), i in [1..nx-1], j in [0..ny-1]
    // index: FxIndex(i,j) with layout (nx-1) * ny
    //
    // Fy is stored on horizontal faces between (i,j-1) and (i,j), i in [0..nx-1], j in [1..ny-1]
    // index: FyIndex(i,j) with layout nx * (ny-1)

    std::size_t FxSize() const { return static_cast<std::size_t>(nx_-1)*ny_; }
    std::size_t FySize() const { return static_cast<std::size_t>(nx_)*(ny_-1); }

    std::size_t FxIndex(int i, int j) const {
        // valid i: 1..nx-1, j: 0..ny-1
        assert(i>=1 && i<=nx_-1 && j>=0 && j<ny_);
        return static_cast<std::size_t>(j)*(nx_-1) + (i-1);
    }

    std::size_t FyIndex(int i, int j) const {
        // valid i: 0..nx-1, j: 1..ny-1
        assert(i>=0 && i<nx_ && j>=1 && j<=ny_-1);
        return static_cast<std::size_t>( (j-1)*nx_ + i );
    }

    // --------- Field registry (cell-centered scalars) ---------

    /**
   * @brief Create (or get) a named cell-centered scalar field.
   * If the name does not exist, it is created and zero-initialized.
   * @return reference to underlying std::vector<double> of size numCells().
   */
    std::vector<double>& field(const std::string& name) {
        auto& v = fields_[name];
        if (v.empty()) v.resize(numCells(), 0.0);
        return v;
    }

    /// const access
    const std::vector<double>& field(const std::string& name) const {
        auto it = fields_.find(name);
        assert(it != fields_.end());
        return it->second;
    }

    /// true if a field with @p name exists
    bool hasField(const std::string& name) const {
        return fields_.find(name) != fields_.end();
    }

    /// erase a field if it exists
    void dropField(const std::string& name) {
        fields_.erase(name);
    }

    // --- Flux fields (defined at interfaces) ---
    std::vector<double>& flux(const std::string& name) {
        return fluxes_[name];
    }
    const std::vector<double>& flux(const std::string& name) const {
        auto it = fluxes_.find(name);
        assert(it != fluxes_.end());
        return it->second;
    }

    // --- Initialization helpers ---
    void addField(const std::string& name, double initVal = 0.0) {
        fields_[name] = std::vector<double>(numCells(), initVal);
    }
    void addFlux(const std::string& name, std::size_t numInterfaces, double initVal = 0.0) {
        fluxes_[name] = std::vector<double>(numInterfaces, initVal);
    }

    // --------- Iteration helpers ---------

    /// visit all cells (i,j)
    void forEachCell(const std::function<void(int i,int j)>& f) const;

    /// visit interior cells only (i ∈ [1..nx-2], j ∈ [1..ny-2])
    void forEachInterior(const std::function<void(int i,int j)>& f) const;

    /// visit vertical faces: (i,j) with i∈[1..nx-1], j∈[0..ny-1]
    void forEachFx(const std::function<void(int i,int j)>& f) const;

    /// visit horizontal faces: (i,j) with i∈[0..nx-1], j∈[1..ny-1]
    void forEachFy(const std::function<void(int i,int j)>& f) const;

    // --------- Face properties (e.g., conductivity) ---------

    /**
   * @brief Harmonic mean of cell-centered properties across a face.
   * Guards small values to keep SPD when used as a weight.
   */
    static inline double harmonic(double a, double b, double floor = 1e-16) {
        if (a < floor) a = floor;
        if (b < floor) b = floor;
        return 2.0*a*b/(a+b);
    }

    // Convenience: compute face-centered properties from a cell field
    // (Kx on vertical faces and Ky on horizontal faces)
    void computeFacePropertyFromCell(const std::string& KcellName,
                                     std::vector<double>& KxFace, // size FxSize()
                                     std::vector<double>& KyFace  // size FySize()
                                     ) const;

    // --------- PETSc interop (optional) ---------
    // Defined in grid.cpp to avoid forcing PETSc headers here.

    /// Pack a cell field into a PETScVector (same row-major ordering).
    void packFieldToPETSc(const std::string& fieldName, PETScVector& v) const;

    /// Unpack a PETScVector back into a named cell field.
    void unpackFieldFromPETSc(const PETScVector& v, const std::string& fieldName);

    /**
 * @brief Sequential Gaussian Simulation (SGS) with exponential covariance.
 *
 * - Marginal: standard normal (mean 0, var 1).
 * - Covariance: C(dx,dy) = exp(-|dx|/lx - |dy|/ly)   (separable exponential).
 * - Path: random permutation of all cells (first seed point is random).
 * - Each new node is conditioned on up to @p nneigh *already simulated* nearest
 *   neighbors using simple kriging (known mean = 0).
 *
 * @param name      Field name to create/overwrite (size = nx*ny).
 * @param lx, ly    Correlation lengths (physical units; >0).
 * @param nneigh    Max # of nearest neighbors to use (e.g. 8, 12, 16).
 * @param seed      RNG seed (used for path & Gaussian draws).
 * @param nugget    Small diagonal added to the kriging system for numerical stability.
 * @param max_ring  (optional) cap on neighbor search ring radius in *cells*;
 *                  if <=0, it expands until enough neighbors or domain boundary.
 *
 * Notes:
 * - Uses GSL for N(0,1) draws; uses Matrix_Arma / Vector_Arma for solves.
 * - For very large grids, keep nneigh modest; complexity is roughly O(N * nneigh^3)
 *   for the kriging solves, with cheap neighbor search per node.
 */
    void makeGaussianFieldSGS(const std::string& name,
                              double lx, double ly,
                              int nneigh,
                              unsigned seed = 0,
                              double nugget = 1e-10,
                              int max_ring = 0);

    // How the named array is laid out on the grid
    enum class ArrayKind { Cell, Fx, Fy };

    /**
   * @brief Write a named array (cell field or flux) as a dense matrix file.
   *
   * For Cell:  rows = ny,     cols = nx.
   * For Fx:    rows = ny,     cols = nx-1   (vertical faces).
   * For Fy:    rows = ny-1,   cols = nx     (horizontal faces).
   *
   * @param name            Field/flux name.
   * @param kind            Layout kind (Cell, Fx, Fy).
   * @param filename        Output file path.
   * @param delimiter       Column separator (default ',').
   * @param precision       FP precision (default 8, scientific).
   * @param include_header  Prepend a one-line header with grid info.
   * @param flipY           If true, write j = ny-1 .. 0 (top→bottom).
   *
   * Throws std::runtime_error on missing name or size mismatch.
   */
    void writeNamedMatrix(const std::string& name, ArrayKind kind,
                          const std::string& filename,
                          char delimiter = ',',
                          int precision = 8,
                          bool include_header = true,
                          bool flipY = false) const;

    /**
   * @brief Same as writeNamedMatrix, but auto-detects kind:
   *        - if name in fields_  → Cell
   *        - else if name in fluxes_ and size==FxSize() → Fx
   *        - else if name in fluxes_ and size==FySize() → Fy
   *        - otherwise throws.
   */
    void writeNamedMatrixAuto(const std::string& name,
                              const std::string& filename,
                              char delimiter = ',',
                              int precision = 8,
                              bool include_header = true,
                              bool flipY = false) const;

    // Convenience (kept for compatibility): write a cell field explicitly.
    void writeFieldMatrix(const std::string& fieldName,
                          const std::string& filename,
                          char delimiter = ',',
                          int precision = 8,
                          bool include_header = true,
                          bool flipY = false) const {
        writeNamedMatrix(fieldName, ArrayKind::Cell, filename, delimiter, precision, include_header, flipY);
    }


private:
    int nx_, ny_;
    double Lx_, Ly_, dx_, dy_;

    // Registry of cell-centered scalar fields by name
    std::unordered_map<std::string, std::vector<double>> fields_;

    // fluxes (can be per-cell, per-face, or per-interface depending on convention)
    std::unordered_map<std::string, std::vector<double>> fluxes_;

    // generic writer
    void writeMatrixRaw(const std::vector<double>& data,
                        int rows, int cols,
                        const std::string& name,
                        const std::string& filename,
                        char delimiter,
                        int precision,
                        bool include_header,
                        bool flipY) const;
};


// separable exponential covariance with unit variance (sill = 1)
static inline double cov_exp_sep(double dx, double dy, double lx, double ly);

// metric consistent with the covariance for nearest-neighbor sorting
static inline double metric_sep(double dx, double dy, double lx, double ly);
