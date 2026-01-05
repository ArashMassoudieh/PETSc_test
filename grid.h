#pragma once
#include <string>
#include <vector>
#include <functional>
#include <unordered_map>
#include <cassert>
#include <cstddef>
#include <petscsys.h>
#include "TimeSeriesSet.h"

class PETScMatrix;


enum class PerturbDir {
    Radial,
    XOnly,
    YOnly
};


class PETScVector; // fwd-decl to match your wrapper name via include order


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
    enum class DerivDir { X, Y };
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

    // Vertical faces (qx): size = (nx_+1) * ny_
    // Valid: i = 0..nx_, j = 0..ny_-1
    inline std::size_t FxIndex(int i, int j) const {
        assert(i >= 0 && i <= nx_ && j >= 0 && j < ny_);
        return static_cast<std::size_t>(j) * (nx_ + 1) + static_cast<std::size_t>(i);
    }

    // Horizontal faces (qy): size = nx_ * (ny_+1)
    // Valid: i = 0..nx_-1, j = 0..ny_
    inline std::size_t FyIndex(int i, int j) const {
        assert(i >= 0 && i < nx_ && j >= 0 && j <= ny_);
        return static_cast<std::size_t>(j) * nx_ + static_cast<std::size_t>(i);
    }

    // And ensure these match:
    inline std::size_t FxSize() const { return static_cast<std::size_t>(nx_ + 1) * ny_; }
    inline std::size_t FySize() const { return static_cast<std::size_t>(nx_) * (ny_ + 1); }


    void computeMassBalanceError(const std::string& fieldName);
    void computeRowSumErrorField(PETScMatrix* A,
                                         double dt,
                                         const std::string& fieldName);
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

    /// true if a flux with @p name exists
    bool hasFlux(const std::string& name) const {
        return fluxes_.find(name) != fluxes_.end();
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

    // grid.h  (inside class Grid2D)

    // Add to grid.h (in the public section of Grid2D class)
    std::pair<double, double> getVelocityAt(double x, double y,
                                            const std::string& qx_name = "qx",
                                            const std::string& qy_name = "qy") const;
public:
    struct NeighborSet {
        std::vector<std::size_t> idx;   // global linear cell indices (j*nx_ + i)
        std::vector<double>      dist;  // anisotropic metric distance to target (same order as idx)
    };


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
                              unsigned long seed = 12345UL,
                              double nugget = 0.0,
                              int max_ring = 0,
                              double tol_var0 = 1e-12);

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

    /// Bilinear interpolation of a stored array at physical (x,y).
    /// - Cell: centers at ((i+0.5)dx, (j+0.5)dy), size nx*ny
    /// - Fx  : vertical faces at (i*dx,      (j+0.5)dy), size (nx+1)*ny
    /// - Fy  : horizontal faces at ((i+0.5)dx, j*dy     ), size nx*(ny+1)
    /// If clamp=true (default), x,y are clamped to the domain box before sampling.
    /// If clamp=false and (x,y) is outside, throws std::out_of_range.
    double interpolate(const std::string& name, ArrayKind kind,
                       double x, double y, bool clamp = true) const;

    // Compute mean and standard deviation of a field
    double fieldMean(const std::string& name, ArrayKind kind) const;
    double fieldStdDev(const std::string& name, ArrayKind kind) const;

    // Normalize values of a field: v -> (v - a)/b
    void normalizeField(const std::string& name, ArrayKind kind, double a, double b);

    // Average of field values at given x (cell field only)
    double fieldAverageAtX(const std::string& name, double x) const;

    void transportStepUpwind(const std::string& name, double dt);

    void DarcySolve(double H_west, double H_east, const std::string &kx_field, const std::string &ky_field, const char* ksp_prefix = nullptr);


    void createExponentialField(const std::string& inputFieldName,
                                        double a, double b,
                                        const std::string& outputFieldName);

    // Converts a field to time-series;
    TimeSeries<double> toSeries(const std::string& inputFieldName) const;

    /**
     * Sample perturbed field values using Gaussian perturbations.
     * Each sample takes a random point in the field, evaluates it,
     * then perturbs the point according to a Gaussian displacement
     * and evaluates the field again.
     *
     * @param fieldName Name of the field or flux
     * @param kind Whether it is a field or flux
     * @param nSamples Number of perturbation samples
     * @param delta Scaling factor for Gaussian perturbation
     * @param seed Random seed (0 = random_device)
     * @param dir Perturbation direction (Radial, XOnly, YOnly)
     *
     * @return A TimeSeries where (t = original value, c = perturbed value)
     */
    TimeSeries<double> sampleGaussianPerturbation(
        const std::string& fieldName,
        ArrayKind kind,
        int nSamples,
        double delta,
        unsigned long seed,
        PerturbDir dir = PerturbDir::Radial  // default
        ) const;

    // Main transport matrix assembly function
    void assembleTransportMatrix(PETScMatrix *A, double dt) const;

    void assembleTransportRHS(PETScVector& b, const std::string& c_field, double dt) const;

    void transportStep(double dt, const char* ksp_prefix);

    void SetVal(const std::string& prop, const double& value);

    // Getter functions for accessing transport properties
    double getDiffusion() const { return diffusion_coeff_; }
    double getPorosity() const { return porosity_; }
    double getLeftBC() const { return c_left_; }

    void SolveTransport(const double& t_end, const double& dt, const char* ksp_prefix = nullptr, int output_interval = 1, const std::string& output_dir = "", const std::string& filename = "", TimeSeriesSet<double>* btc_data = nullptr);

    void printSampleC(const std::vector<std::pair<int,int>>& pts) const;

    std::pair<double,double> fieldMinMax(const std::string& name, ArrayKind kind) const;

    void assignConstant(const std::string& name, ArrayKind kind, double value);

    /**
     * @brief Sample second derivatives of a field at random points.
     *
     * Picks nPoints random (x,y) in the domain, evaluates the field g at
     * (x-Δ, x, x+Δ) or (y-Δ, y, y+Δ) depending on dir, and approximates the
     * second derivative using central differences.
     *
     * Result is returned as a TimeSeries<double>, where:
     *   - t = g(x,y) (field value at the sample point)
     *   - c = second derivative wrt x or y
     *
     * @param fieldName  Name of field or flux to sample.
     * @param kind       ArrayKind (Cell, Fx, Fy).
     * @param dir        DerivDir::X or DerivDir::Y.
     * @param nPoints    Number of random samples.
     * @param delta      Step size Δ for finite difference.
     * @param seed       RNG seed (optional, 0 = random device).
     * @return           TimeSeries<double> of (value, second derivative).
     */
    TimeSeries<double> sampleSecondDerivative(
        const std::string& fieldName,
        ArrayKind kind,
        DerivDir dir,
        int nPoints,
        double delta,
        unsigned long seed = 0
        ) const;

    /**
 * @brief Export all values of a field into a TimeSeries<double>.
 *
 * Each entry has:
 *   - t = counter index (0,1,2,...)
 *   - c = field value
 *
 * @param fieldName  Name of field or flux
 * @param kind       ArrayKind (Cell, Fx, Fy)
 * @return           TimeSeries<double> of values
 */
    TimeSeries<double> exportFieldToTimeSeries(
        const std::string& fieldName,
        ArrayKind kind
        ) const;


    /**
 * @brief Assign values from a TimeSeries<double> into a new field/flux.
 *
 * The TimeSeries length must match the size of the target ArrayKind:
 *   - Cell: nx * ny
 *   - Fx:   (nx+1) * ny
 *   - Fy:   nx * (ny+1)
 *
 * @param ts         Input TimeSeries<double> (values in .c).
 * @param fieldName  Name of the field/flux to create/overwrite.
 * @param kind       ArrayKind (Cell, Fx, or Fy).
 */
    void assignFromTimeSeries(
        const TimeSeries<double>& ts,
        const std::string& fieldName,
        ArrayKind kind
        );
#ifdef GRID_USE_VTK
    /**
 * @brief Write a named array (cell field or face flux) to a .vti file using VTK.
 *
 * ArrayKind::Cell → cells = (nx × ny)
 * ArrayKind::Fx   → cells = ((nx+1) × ny)      // vertical faces
 * ArrayKind::Fy   → cells = (nx × (ny+1))      // horizontal faces
 *
 * The data are written as CellData (one scalar per cell).
 *
 * @param name       Key in fields_ (Cell) or fluxes_ (Fx/Fy).
 * @param kind       Layout of the named array.
 * @param filename   Path to .vti file.
 * @param arrayName  Optional name for the VTK array (defaults to 'name').
 * @param flipY      If true, reverse Y row order when writing (visual preference).
 */
    void writeNamedVTI(const std::string& name,
                       ArrayKind kind,
                       const std::string& filename,
                       const std::string& arrayName = "",
                       bool flipY = false) const;

    /** Convenience: auto-detect Cell / Fx / Fy from maps + sizes. */
    void writeNamedVTI_Auto(const std::string& name,
                            const std::string& filename,
                            const std::string& arrayName = "",
                            bool flipY = false) const;
#endif // GRID_USE_VTK

    // Add to grid.h - in the Grid2D class public section:

    // Normal distribution functions (using GSL)
    static double phi(double z);           // Standard normal PDF
    static double Phi(double z);           // Standard normal CDF
    static double Phi_inv(double u);       // Inverse standard normal CDF (probit)

    // Mixing kernel
    double kappa(double v, double lc, double lambda_x, double lambda_y) const;

    // Parameters for mixing equation
    void setMixingParams(double lc, double lambda_x, double lambda_y);
    double lc() const { return lc_; }
    double lambda_x() const { return lambda_x_; }
    double lambda_y() const { return lambda_y_; }

    // Add to grid.h public section:

    // Convert concentration field to uniform scores (u-space)
    void convertToUniformScore(const std::string& c_field_name,
                               const std::string& u_field_name);

    // Convert uniform scores back to concentration
    void convertFromUniformScore(const std::string& u_field_name,
                                 const std::string& c_field_name);

    // Get velocity at a given u-value (for advection in u-space)
    double velocityAtU(double x, double u) const;

    // Solve mixing PDF equation in (x,u) space
    void SolveMixingPDF(const double& t_end,
                        const double& dt,
                        const char* ksp_prefix,
                        int output_interval = 0,
                        const std::string& output_dir = "");

    // Assemble the mixing PDF matrix
    void assembleMixingPDFMatrix(PETScMatrix* A, double dt);

    // Single time step for mixing PDF
    void mixingPDFStep(double dt, const char* ksp_prefix);

     /**
     * @brief Extract empirical CDF from a field by integrating the PDF
     * @param field_name Name of the field
     * @param kind Field kind (Cell, Fx, Fy)
     * @param num_bins Number of bins for the distribution
     * @return TimeSeries where time=value and value=cumulative probability
     */
    TimeSeries<double> extractFieldCDF(const std::string& field_name, ArrayKind kind, int num_bins = 100) const;

    /**
    * @brief Get field value at a given cumulative probability
    * @param field_name Name of the field
    * @param kind Field kind (Cell, Fx, Fy)
    * @param u Cumulative probability [0,1]
    * @return Field value at probability u
    */
    double getFieldValueAtCDF(const std::string& field_name, ArrayKind kind, double u) const;


    void computeMixingDiffusionCoefficient();

    double getFieldValueAt(const std::string& property_name, double x, double y) const;
    double getAverageAlongY(const std::string& property_name, double x) const;

    void setBTCLocations(const std::vector<double>& locations);
    const std::vector<double>& getBTCLocations() const;
    std::vector<double>& getBTCLocations();

private:
    int nx_, ny_;
    double Lx_, Ly_, dx_, dy_;

    // Transport properties
    double diffusion_coeff_;     // Diffusion coefficient D
    double porosity_;           // Porosity
    double c_left_;            // Left boundary concentration
    double lc_;          // Correlation length scale
    double lambda_x_;    // Transverse dispersion length in x
    double lambda_y_;    // Transverse dispersion length in y
    std::vector<double> BTCLocations_; // Breakthrough curve monitoring locations

    // Registry of cell-centered scalar fields by name
    std::unordered_map<std::string, std::vector<double>> fields_;

    std::string pdf_field_name_;

    // fluxes (can be per-cell, per-face, or per-interface depending on convention)
    std::unordered_map<std::string, std::vector<double>> fluxes_;

    // generic writer
    void writeMatrixRaw(const std::vector<double>& data,
                                int rows, int cols,
                                const std::string& name,
                                const std::string& filename,
                                char delimiter = ',',
                                int precision = 3,
                                bool include_header = true,
                                bool flipY = false,
                                ArrayKind kind = Grid2D::ArrayKind::Cell) const;

    /**
   * @brief Collect the n nearest *determined* neighbors of (it,jt) using ring expansion,
   *        then select the true top-n by the anisotropic metric.
   *
   * The function keeps expanding rings until:
   *  - we have at least `nneigh` candidates, and
   *  - if the target row has any known cells, at least one **same-row** neighbor
   *    is included in the candidate set (ensures correctness for lx >> Lx).
   *
   * @param it,jt       target cell indices
   * @param nneigh      number of neighbors to return (<=0 returns empty)
   * @param lx,ly       correlation lengths (for metric)
   * @param max_ring    if >0, cap ring radius in cells; if 0, expand up to max(nx_,ny_)
   * @param known       bitmap of which cells are already simulated (size nx_*ny_)
   * @param ensure_same_row  if true, enforces presence of a same-row neighbor when one exists
   * @return NeighborSet {idx, dist}, both sized <= nneigh
   */
    NeighborSet gatherNeighbors(int it, int jt,
                                int nneigh,
                                double lx, double ly,
                                int max_ring,
                                const std::vector<unsigned char>& known,
                                bool ensure_same_row = true) const;

    // Transport equation utility functions
    double getVelocityX(int i, int j, double porosity) const;
    double getVelocityY(int i, int j, double porosity) const;
    PetscInt cellToPetscIndex(int i, int j) const;

    void addTransportXTerms(int i, int j, double dt, double D, double porosity,
                            double& diag, std::vector<PetscInt>& cols,
                            std::vector<PetscScalar>& vals) const;

    void addTransportYTerms(int i, int j, double dt, double D, double porosity,
                            double& diag, std::vector<PetscInt>& cols,
                            std::vector<PetscScalar>& vals) const;
    PETScMatrix *A = nullptr; //Transport Matrix
};

// harmonic mean
static inline double harm(double a, double b);

// separable exponential covariance with unit variance (sill = 1)
static inline double cov_exp_sep(double dx, double dy, double lx, double ly);

// metric consistent with the covariance for nearest-neighbor sorting
static inline double metric_sep(double dx, double dy, double lx, double ly);
