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

    // --------- Flux storage (face-centered) ---------

    /// allocate/zero the face-centered flux arrays
    void allocateFluxes() {
        Fx_.assign(FxSize(), 0.0);
        Fy_.assign(FySize(), 0.0);
    }

    std::vector<double>& Fx() { return Fx_; }
    std::vector<double>& Fy() { return Fy_; }
    const std::vector<double>& Fx() const { return Fx_; }
    const std::vector<double>& Fy() const { return Fy_; }

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

private:
    int nx_, ny_;
    double Lx_, Ly_, dx_, dy_;

    // Registry of cell-centered scalar fields by name
    std::unordered_map<std::string, std::vector<double>> fields_;

    // Face-centered flux arrays (optional; allocate when needed)
    std::vector<double> Fx_, Fy_;
};
