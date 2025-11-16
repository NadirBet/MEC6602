# O-Grid Mesh System Documentation

## Table of Contents
1. [Overview](#overview)
2. [Topology and Coordinate System](#topology-and-coordinate-system)
3. [Memory Layout and Indexing](#memory-layout-and-indexing)
4. [Ghost Cell Philosophy](#ghost-cell-philosophy)
5. [Geometric Computations](#geometric-computations)
6. [Face Geometry and Normals](#face-geometry-and-normals)
7. [Spectral Radius](#spectral-radius)
8. [Implementation Flow](#implementation-flow)
9. [Key Equations Reference](#key-equations-reference)

---

## 1. Overview

### Purpose
This mesh system provides a **structured O-grid** for 2D Euler equation solvers around airfoils. It implements:
- Finite volume cell-centered storage
- Ghost layers for boundary conditions
- Shared-face architecture for flux computation
- Full geometric metric calculations

### Key Features
- **O-grid topology**: Wraps circumferentially around airfoil
- **Periodic i-direction**: Trailing edge connectivity
- **Radial j-direction**: Body surface to far-field
- **Ghost layers**: 2 layers (`NGHOST = 2`) on all sides
- **Normalization**: Automatic non-dimensionalization

---

## 2. Topology and Coordinate System

### 2.1 O-Grid Structure

```
                    Far-field boundary (j_max)
                   ╱────────────────────────╲
                 ╱          Radial            ╲
               ╱           direction            ╲
             ╱              (j)                   ╲
           ╱                                        ╲
          │              ┌─────────┐                │
          │            ╱             ╲              │
          │           │   Airfoil     │             │
          │           │   (Body)      │  ←──────────┼─ i wraps (periodic)
          │            ╲             ╱              │
          │              └─────────┘                │
           ╲                                        ╱
             ╲                                    ╱
               ╲                                ╱
                 ╲                            ╱
                   ╲────────────────────────╱
                    Body surface (j = 0)
```

### 2.2 Coordinate Directions

| Direction | Index | Physical Meaning | Boundary Condition |
|-----------|-------|------------------|-------------------|
| **i** | Circumferential | Wraps around airfoil | **Periodic** (connects at trailing edge) |
| **j** | Radial | Body → Far-field | **Wall at j=0**, **Freestream at j=max** |

### 2.3 Domain Dimensions

```cpp
// Physical domain (actual flow field)
int ni;  // Circumferential cells
int nj;  // Radial cells

// Total domain (with ghosts)
int niTotal = ni + 2*NGHOST;  // ni + 4
int njTotal = nj + 2*NGHOST;  // nj + 4

// Nodes (vertices)
int niNodes = niTotal + 1;    // One more than cells
int njNodes = njTotal + 1;
```

### 2.4 Physical vs Ghost Regions

```
j-direction (radial):
┌─────────────────────────────────────┐
│  Top ghost (far-field)              │  j ∈ [nj+2, nj+3]
├─────────────────────────────────────┤
│                                     │
│  Physical domain                    │  j ∈ [2, nj+1]  (NGHOST to NGHOST+nj-1)
│  (airfoil flow field)               │
│                                     │
├─────────────────────────────────────┤
│  Bottom ghost (body)                │  j ∈ [0, 1]
└─────────────────────────────────────┘

i-direction (circumferential):
  Left ghost  │    Physical    │  Right ghost
  (periodic)  │     domain     │  (periodic)
    i ∈ [0,1] │  i ∈ [2,ni+1] │  i ∈ [ni+2,ni+3]
```

---

## 3. Memory Layout and Indexing

### 3.1 Index Mapping Functions

All arrays use **row-major (C-style) ordering**:

#### Node Index
```cpp
int nodeIndex(int i, int j) const {
    return j * niNodes + i;
}
```
- **Nodes**: Corner points of cells
- **Range**: `i ∈ [0, niTotal]`, `j ∈ [0, njTotal]`
- **Total count**: `niNodes × njNodes = (niTotal+1) × (njTotal+1)`

#### Cell Index
```cpp
int cellIndex(int i, int j) const {
    return j * niTotal + i;
}
```
- **Cells**: Control volumes (quadrilaterals)
- **Range**: `i ∈ [0, niTotal-1]`, `j ∈ [0, njTotal-1]`
- **Total count**: `niTotal × njTotal`

### 3.2 Shared-Face Indexing

Each face is stored **once** and shared between adjacent cells.

#### I-Face Index (Circumferential faces)
```cpp
int iFaceIndex(int i, int j) const {
    return j * niTotal + i;
}
```
- **Connects**: Node `(i,j)` → Node `(i+1,j)`
- **Orientation**: Horizontal (constant j)
- **Range**: `i ∈ [0, niTotal-1]`, `j ∈ [0, njTotal]`
- **Total count**: `niTotal × (njTotal+1)`
- **Separates**: Cells `(i, j-1)` and `(i, j)` in radial direction

```
Cell layout for i-faces:
         j+1  o────────o
              │  (i,j)  │
         j    o────f────o  ← i-face at (i,j)
              │ (i,j-1) │
         j-1  o────────o
              i       i+1
```

#### J-Face Index (Radial faces)
```cpp
int jFaceIndex(int i, int j) const {
    return j * niNodes + i;
}
```
- **Connects**: Node `(i,j)` → Node `(i,j+1)`
- **Orientation**: Vertical (constant i)
- **Range**: `i ∈ [0, niNodes-1]`, `j ∈ [0, njTotal-1]`
- **Total count**: `niNodes × njTotal`
- **Separates**: Cells `(i-1, j)` and `(i, j)` in circumferential direction

```
Cell layout for j-faces:
              i+1
               │
         j+1   o
               │
               f  ← j-face at (i,j)
               │
         j     o
               │
         (i-1,j) | (i,j)
```

### 3.3 Interior Bounds Helper

```cpp
void getInteriorBounds(int& imin, int& imax, int& jmin, int& jmax) const {
    imin = NGHOST;              // 2
    imax = NGHOST + ni;         // ni + 2
    jmin = NGHOST;              // 2
    jmax = NGHOST + nj;         // nj + 2
}
```
Returns loop bounds for **physical cells only**.

---

## 4. Ghost Cell Philosophy

### 4.1 Why Ghost Cells Exist

Ghost cells are **numerical devices** for:
1. **Boundary condition implementation**: Provide "other side" for flux calculation
2. **Uniform stencil**: Same flux formula at boundaries and interior
3. **High-order reconstruction**: Support MUSCL, gradient calculations
4. **Avoid special cases**: Code treats all cells identically

### 4.2 Ghost Cells Are NOT Physical Space

**Critical Understanding**:
- Ghost cells at `j=0,1` (below body) represent **fictitious space "inside" the airfoil**
- They are **mathematical constructs** only
- Their geometry supports numerical operations, not physical interpretation
- The **wall face at j=2** has true airfoil geometry from the mesh file

### 4.3 What Ghost Cells Provide

| Property | Computed? | Purpose |
|----------|-----------|---------|
| Ghost **nodes** (x,y) | ✅ Yes | Foundation for all geometry |
| Ghost **cell centers** | ✅ Yes | Distance calculations, reconstruction |
| Ghost **cell areas** | ✅ Yes | Finite volume scheme consistency |
| Ghost **face geometry** | ✅ Yes | Flux calculation at boundaries |
| Ghost **cell states** (ρ,u,v,p) | Set by solver | Enforce boundary conditions |
| Ghost **spectral radius** | ❌ No | Not time-marched |

### 4.4 Ghost Node Construction

Ghost nodes are filled by **extrapolation** or **copying** from physical nodes:

```
Step 1: Read physical nodes from PLOT3D
Step 2: Fill i-periodic ghosts (copy)
Step 3: Fill j-extrapolated ghosts (symmetric)
Step 4: Compute geometry for ALL cells using ghost nodes
```

---

## 5. Geometric Computations

### 5.1 Mesh Initialization Sequence

```cpp
// From readPlot3D():
1. allocate(ni, nj);              // Allocate arrays
2. Read nodes into physical region
3. normalizePhysicalNodes(Lref);  // Non-dimensionalize
4. fillGhostNodes();              // Extrapolate/copy
5. computeCellCenters();          // Average of 4 corners
6. computeCellAreas();            // Shoelace formula
7. computeFaceGeometry();         // Normals and lengths
```

### 5.2 Node Filling: I-Direction (Periodic)

**Executed for ALL j-levels (including undefined ghost rows at this stage)**:

```cpp
for (int j = 0; j < njNodes; ++j) {
    // Left ghosts (i=0,1)
    for (int g = 0; g < NGHOST; ++g) {
        int ig = imin - 1 - g;     // i = 1, 0
        int is = imax - g;         // i = ni+2, ni+1
        xNodes[nodeIndex(ig, j)] = xNodes[nodeIndex(is, j)];
        yNodes[nodeIndex(ig, j)] = yNodes[nodeIndex(is, j)];
    }
    
    // Right ghosts (i=ni+3, ni+4)
    for (int g = 0; g < NGHOST; ++g) {
        int ig = imax + 1 + g;     // i = ni+3, ni+4
        int is = imin + g;         // i = 2, 3
        xNodes[nodeIndex(ig, j)] = xNodes[nodeIndex(is, j)];
        yNodes[nodeIndex(ig, j)] = yNodes[nodeIndex(is, j)];
    }
}
```

**Result**: Circumferential wraparound connectivity established.

**Mathematical interpretation**:
```
Node(0, j) ≡ Node(ni+2, j)  (trailing edge connection)
Node(1, j) ≡ Node(ni+1, j)
```

### 5.3 Node Filling: J-Direction (Symmetric Extrapolation)

**Bottom ghosts** (near body, j=0,1):

```cpp
for (int i = 0; i < niNodes; ++i) {
    for (int g = 0; g < NGHOST; ++g) {
        int jg = jmin - 1 - g;       // j = 1, 0
        int j1 = jmin;               // j = 2 (body surface)
        int j2 = jmin + 1 + g;       // j = 3, 4
        
        xNodes[nodeIndex(i, jg)] = 2.0*xNodes[nodeIndex(i, j1)] 
                                    - xNodes[nodeIndex(i, j2)];
        yNodes[nodeIndex(i, jg)] = 2.0*yNodes[nodeIndex(i, j1)] 
                                    - yNodes[nodeIndex(i, j2)];
    }
}
```

**Equation**:
$$
\vec{r}_{\text{ghost}} = 2\vec{r}_{\text{wall}} - \vec{r}_{\text{interior}}
$$

**Geometric interpretation**:
```
j=4   •             Interior point
      │
j=3   •             Interior point
      │
j=2   •  ───────    Body surface (WALL)
      │
j=1   •             Ghost (mirrored)
      │
j=0   •             Ghost (mirrored further)
```

This creates a **symmetric reflection** about the wall at j=2.

**Top ghosts** (far-field, j=nj+3, nj+4):

```cpp
for (int i = 0; i < niNodes; ++i) {
    for (int g = 0; g < NGHOST; ++g) {
        int jg = jmax + 1 + g;       // j = nj+3, nj+4
        int j1 = jmax - g;           // j = nj+2, nj+1
        int j2 = jmax - 1 - g;       // j = nj+1, nj
        
        xNodes[nodeIndex(i, jg)] = 2.0*xNodes[nodeIndex(i, j1)] 
                                    - xNodes[nodeIndex(i, j2)];
        yNodes[nodeIndex(i, jg)] = 2.0*yNodes[nodeIndex(i, j1)] 
                                    - yNodes[nodeIndex(i, j2)];
    }
}
```

Same extrapolation outward to far-field.

### 5.4 Cell Centers

**For ALL cells** (physical + ghost):

```cpp
void computeCellCenters() {
    for (int j = 0; j < njTotal; ++j) {
        for (int i = 0; i < niTotal; ++i) {
            int n00 = nodeIndex(i,     j);
            int n10 = nodeIndex(i + 1, j);
            int n01 = nodeIndex(i,     j + 1);
            int n11 = nodeIndex(i + 1, j + 1);
            
            int c = cellIndex(i, j);
            xCells[c] = 0.25 * (xNodes[n00] + xNodes[n10] + 
                                xNodes[n01] + xNodes[n11]);
            yCells[c] = 0.25 * (yNodes[n00] + yNodes[n10] + 
                                yNodes[n01] + yNodes[n11]);
        }
    }
}
```

**Equation**:
$$
\vec{r}_{\text{cell}} = \frac{1}{4}\left(\vec{r}_{00} + \vec{r}_{10} + \vec{r}_{01} + \vec{r}_{11}\right)
$$

```
Quadrilateral cell:
  n01 •────────• n11
      │        │
      │   •xc  │  ← Cell center
      │        │
  n00 •────────• n10
```

### 5.5 Cell Areas (Shoelace Formula)

**For ALL cells** (physical + ghost):

```cpp
void computeCellAreas() {
    for (int j = 0; j < njTotal; ++j) {
        for (int i = 0; i < niTotal; ++i) {
            // Get 4 corner coordinates
            double x00 = xNodes[n00], y00 = yNodes[n00];
            double x10 = xNodes[n10], y10 = yNodes[n10];
            double x01 = xNodes[n01], y01 = yNodes[n01];
            double x11 = xNodes[n11], y11 = yNodes[n11];
            
            double A = 0.5 * std::abs(
                x00*(y10 - y01) +
                x10*(y11 - y00) +
                x11*(y01 - y10) +
                x01*(y00 - y11)
            );
            cellArea[cellIndex(i,j)] = A;
        }
    }
}
```

**Shoelace formula** (Green's theorem for polygon area):

$$
A = \frac{1}{2}\left|\sum_{k=0}^{n-1}(x_k y_{k+1} - x_{k+1} y_k)\right|
$$

For a quadrilateral traversed as 00→10→11→01→00:

$$
A = \frac{1}{2}\left|x_{00}(y_{10}-y_{01}) + x_{10}(y_{11}-y_{00}) + x_{11}(y_{01}-y_{10}) + x_{01}(y_{00}-y_{11})\right|
$$

**Works for arbitrary convex/concave quads** (essential for curved O-grid cells).

### 5.6 Normalization

```cpp
void normalizePhysicalNodes(double Lref) {
    if (Lref <= 0.0) return;
    for (double &x : xNodes) x /= Lref;
    for (double &y : yNodes) y /= Lref;
}
```

Where:
$$
L_{\text{ref}} = \max(L_x, L_y)
$$

$$
L_x = x_{\max} - x_{\min}, \quad L_y = y_{\max} - y_{\min}
$$

**Purpose**: 
- Non-dimensional coordinates ∈ O(1)
- Improves numerical stability in Euler solver
- Avoids roundoff errors from large coordinate values

---

## 6. Face Geometry and Normals

### 6.1 Face Normal Convention

**Critical design choice**: Normals computed using right-hand rule relative to local coordinate direction.

### 6.2 I-Faces (Circumferential)

```cpp
void computeFaceGeometry() {
    for (int j = 0; j <= njTotal; ++j) {
        for (int i = 0; i < niTotal; ++i) {
            int f = iFaceIndex(i, j);
            int n0 = nodeIndex(i,     j);
            int n1 = nodeIndex(i + 1, j);
            
            double dx = xNodes[n1] - xNodes[n0];
            double dy = yNodes[n1] - yNodes[n0];
            double L  = std::sqrt(dx*dx + dy*dy);
            iFaceLen[f] = L;
            
            if (L > 0.0) {
                iFaceNormal[f] = { -dy/L,  dx/L };  // Unit normal
            }
        }
    }
}
```

**Face vector and normal**:

```
       n1
       •
      ╱│
  dy ╱ │ 
    ╱  │ n̂ = (-dy/L, dx/L)
   •───→
  n0  dx

Tangent vector: t̂ = (dx/L, dy/L)
Normal vector:  n̂ = (-dy/L, dx/L) = 90° CCW rotation of tangent
```

**Equations**:
$$
\Delta x = x_1 - x_0, \quad \Delta y = y_1 - y_0
$$

$$
L = \sqrt{\Delta x^2 + \Delta y^2}
$$

$$
\hat{n} = \frac{1}{L}\begin{pmatrix} -\Delta y \\ \Delta x \end{pmatrix}
$$

**Normal orientation**: Points toward **increasing j** (radially outward).

### 6.3 J-Faces (Radial)

```cpp
for (int j = 0; j < njTotal; ++j) {
    for (int i = 0; i < niNodes; ++i) {
        int f = jFaceIndex(i, j);
        int n0 = nodeIndex(i, j);
        int n1 = nodeIndex(i, j + 1);
        
        double dx = xNodes[n1] - xNodes[n0];
        double dy = yNodes[n1] - yNodes[n0];
        double L  = std::sqrt(dx*dx + dy*dy);
        jFaceLen[f] = L;
        
        if (L > 0.0) {
            jFaceNormal[f] = {  dy/L, -dx/L };  // Unit normal
        }
    }
}
```

**Face vector and normal**:

```
  n1
  •
  │╲
  │ ╲ n̂ = (dy/L, -dx/L)
  │  ╲
  •   
  n0

Tangent vector: t̂ = (dx/L, dy/L)
Normal vector:  n̂ = (dy/L, -dx/L) = 90° CW rotation of tangent
```

**Equations**:
$$
\hat{n} = \frac{1}{L}\begin{pmatrix} \Delta y \\ -\Delta x \end{pmatrix}
$$

**Normal orientation**: Points toward **increasing i** (circumferentially).

### 6.4 Normal Direction Verification

For a cell at (i,j):

```
Cell (i,j) bounded by 4 faces:

  j+1  o────f_top────o
       │             │
  f_L  │  Cell (i,j) │ f_R
       │             │
  j    o───f_bot─────o
       i           i+1

where:
  f_bot = i-face at (i,   j)   → normal points UP    (toward j+1)
  f_top = i-face at (i, j+1)   → normal points UP    (toward j+2)
  f_L   = j-face at (i,   j)   → normal points RIGHT (toward i+1)
  f_R   = j-face at (i+1, j)   → normal points RIGHT (toward i+2)
```

**All normals point "outward" from the cell in the positive coordinate direction.**

---

## 7. Spectral Radius

### 7.1 Purpose

The spectral radius `λ` represents the **maximum characteristic wave speed** at a cell, used for:
1. **CFL condition**: `Δt ≤ CFL × V/λ`
2. **Time-step calculation** for explicit Euler solver
3. **Artificial dissipation scaling**

### 7.2 Computation

```cpp
void computeSpectralRadius(const std::vector<Primitive>& W, double gamma) {
    lambdaI.assign(niTotal * njTotal, 0.0);
    lambdaJ.assign(niTotal * njTotal, 0.0);
    lambda.assign(niTotal * njTotal, 0.0);
    
    int imin, imax, jmin, jmax;
    getInteriorBounds(imin, imax, jmin, jmax);
    
    for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
            int c = cellIndex(i, j);
            
            // I-direction contribution
            int fiL = iFaceIndex(i, j);
            int fiR = iFaceIndex(i, j + 1);
            
            double nxL = iFaceNormal[fiL][0], nyL = iFaceNormal[fiL][1];
            double nxR = iFaceNormal[fiR][0], nyR = iFaceNormal[fiR][1];
            double LiL = iFaceLen[fiL],        LiR = iFaceLen[fiR];
            
            double unL = std::abs(W[c].u * nxL + W[c].v * nyL);
            double unR = std::abs(W[c].u * nxR + W[c].v * nyR);
            double a   = W[c].a;  // Speed of sound
            
            lambdaI[c] = 0.5 * ((unL + a)*LiL + (unR + a)*LiR);
            
            // J-direction contribution
            int fjB = jFaceIndex(i,     j);
            int fjT = jFaceIndex(i + 1, j);
            
            double nxB = jFaceNormal[fjB][0], nyB = jFaceNormal[fjB][1];
            double nxT = jFaceNormal[fjT][0], nyT = jFaceNormal[fjT][1];
            double LjB = jFaceLen[fjB],        LjT = jFaceLen[fjT];
            
            double unB = std::abs(W[c].u * nxB + W[c].v * nyB);
            double unT = std::abs(W[c].u * nxT + W[c].v * nyT);
            
            lambdaJ[c] = 0.5 * ((unB + a)*LjB + (unT + a)*LjT);
            
            lambda[c] = lambdaI[c] + lambdaJ[c];
        }
    }
}
```

### 7.3 Theory

For 2D Euler equations, the eigenvalues of the flux Jacobian are:

$$
\lambda = u \cdot \hat{n} \pm a, \quad u \cdot \hat{n}
$$

where:
- $u \cdot \hat{n} = u n_x + v n_y$ = normal velocity component
- $a = \sqrt{\gamma p/\rho}$ = speed of sound

**Maximum wave speed** through a face:

$$
\lambda_{\text{max}} = |u \cdot \hat{n}| + a
$$

### 7.4 Cell Spectral Radius

For each cell, sum contributions from all faces:

**I-direction** (2 i-faces per cell):
$$
\lambda_I = \frac{1}{2}\left[(|u \cdot \hat{n}_L| + a)L_L + (|u \cdot \hat{n}_R| + a)L_R\right]
$$

**J-direction** (2 j-faces per cell):
$$
\lambda_J = \frac{1}{2}\left[(|u \cdot \hat{n}_B| + a)L_B + (|u \cdot \hat{n}_T| + a)L_T\right]
$$

**Total**:
$$
\lambda_{\text{total}} = \lambda_I + \lambda_J
$$

**Note**: Only computed for **physical cells** (ghosts not time-marched).

### 7.5 Time-Step Calculation

In the Euler solver:

$$
\Delta t = \text{CFL} \times \min_{\text{all cells}} \left(\frac{V_{\text{cell}}}{\lambda_{\text{cell}}}\right)
$$

where:
- $V_{\text{cell}}$ = cell area (2D)
- $\text{CFL} \in (0, 1]$ = Courant number

---

## 8. Implementation Flow

### 8.1 Complete Initialization Sequence

```
1. allocate(ni, nj)
   ├─ Allocate node arrays:  (ni+5) × (nj+5)
   ├─ Allocate cell arrays:  (ni+4) × (nj+4)
   ├─ Allocate i-face arrays: (ni+4) × (nj+5)
   └─ Allocate j-face arrays: (ni+5) × (nj+4)

2. readPlot3D(filename)
   ├─ Read physical nodes → region [2:ni+2, 2:nj+2]
   ├─ Find bounding box (Lx, Ly)
   ├─ Normalize: Lref = max(Lx, Ly)
   └─ Scale all nodes by 1/Lref

3. fillGhostNodes()
   ├─ I-periodic: copy nodes from opposite side
   │  ├─ Left ghost [0:1, :] ← Right physical [ni+1:ni+2, :]
   │  └─ Right ghost [ni+3:ni+4, :] ← Left physical [2:3, :]
   │
   └─ J-extrapolate: symmetric reflection
      ├─ Bottom ghost [:, 0:1] ← Extrapolate from [:, 2:4]
      └─ Top ghost [:, nj+3:nj+4] ← Extrapolate from [:, nj:nj+2]

4. computeCellCenters()
   └─ For all cells: average 4 corner nodes

5. computeCellAreas()
   └─ For all cells: shoelace formula on 4 corners

6. computeFaceGeometry()
   ├─ For all i-faces: compute length and normal
   └─ For all j-faces: compute length and normal
```

### 8.2 Usage in Euler Solver

```cpp
// In solver initialization:
Mesh mesh;
mesh.readPlot3D("ogrid.xyz");

// Before time loop:
mesh.computeSpectralRadius(W, gamma);

// Time-stepping:
double dt = computeTimeStep(mesh.lambda, mesh.cellArea, CFL);

// Flux computation at face (i,j):
int face = mesh.iFaceIndex(i, j);
double nx = mesh.iFaceNormal[face][0];
double ny = mesh.iFaceNormal[face][1];
double L  = mesh.iFaceLen[face];

Flux F = computeEulerFlux(WL, WR, nx, ny);
Residual[c] += F * L;  // Multiply by face length
```

---

## 9. Key Equations Reference

### 9.1 Index Mappings

| Entity | Formula | Range |
|--------|---------|-------|
| Node | `j × niNodes + i` | `i ∈ [0, niTotal]`, `j ∈ [0, njTotal]` |
| Cell | `j × niTotal + i` | `i ∈ [0, niTotal-1]`, `j ∈ [0, njTotal-1]` |
| I-face | `j × niTotal + i` | `i ∈ [0, niTotal-1]`, `j ∈ [0, njTotal]` |
| J-face | `j × niNodes + i` | `i ∈ [0, niNodes-1]`, `j ∈ [0, njTotal-1]` |

### 9.2 Ghost Node Extrapolation

**Symmetric reflection**:
$$
\vec{r}_{\text{ghost}} = 2\vec{r}_{\text{boundary}} - \vec{r}_{\text{interior}}
$$

**Periodic copy** (i-direction):
$$
\vec{r}(0, j) = \vec{r}(n_i + 2, j), \quad \vec{r}(n_i + 3, j) = \vec{r}(2, j)
$$

### 9.3 Cell Center

$$
\vec{r}_c = \frac{1}{4}\sum_{k=0}^{3} \vec{r}_k
$$

### 9.4 Cell Area (Shoelace)

$$
A = \frac{1}{2}\left|x_0(y_1 - y_3) + x_1(y_2 - y_0) + x_2(y_3 - y_1) + x_3(y_0 - y_2)\right|
$$

### 9.5 Face Normal and Length

For edge from $\vec{r}_0$ to $\vec{r}_1$:

$$
\vec{\Delta} = \vec{r}_1 - \vec{r}_0 = (\Delta x, \Delta y)
$$

$$
L = \|\vec{\Delta}\| = \sqrt{\Delta x^2 + \Delta y^2}
$$

**I-face normal** (CCW rotation):
$$
\hat{n} = \frac{1}{L}(-\Delta y, \Delta x)
$$

**J-face normal** (CW rotation):
$$
\hat{n} = \frac{1}{L}(\Delta y, -\Delta x)
$$

### 9.6 Spectral Radius

$$
\lambda = \sum_{\text{faces}} \frac{1}{2}(|\vec{u} \cdot \hat{n}| + a) L
$$

where $a = \sqrt{\gamma p / \rho}$ is the speed of sound.

### 9.7 CFL Time-Step

$$
\Delta t = \text{CFL} \times \min_{\text{cells}} \left(\frac{A_{\text{cell}}}{\lambda_{\text{cell}}}\right)
$$

---

## 10. Important Implementation Notes

### 10.1 Ghost Cell States (Set by Euler Solver, Not Mesh)

The mesh provides **geometry only**. The Euler solver must set ghost cell **flow states**:

**Wall BC** (j=2 face, between ghost j=1 and physical j=2):
```cpp
// Mirror velocity, extrapolate scalars
W_ghost.rho = W_physical.rho;
W_ghost.p   = W_physical.p;
vec2 u_phys = {W_physical.u, W_physical.v};
vec2 u_ghost = u_phys - 2.0*(dot(u_phys, n_wall))*n_wall;  // Reflect
W_ghost.u = u_ghost[0];
W_ghost.v = u_ghost[1];
```

**Far-field BC** (j=nj+2 face):
```cpp
// Set to freestream
W_ghost.rho = rho_inf;
W_ghost.u   = u_inf;
W_ghost.v   = v_inf;
W_ghost.p   = p_inf;
```

**Periodic BC** (automatic from geometry):
- Ghost cells at i=0,1 get states from physical i=ni+1, ni+2
- Ghost cells at i=ni+3, ni+4 get states from physical i=2, 3

### 10.2 VTK Output

**Physical cells only**:
```cpp
mesh.writeVTKPhysical("output.vtk");  // Only i ∈ [2,ni+1], j ∈ [2,nj+1]
```

**Face visualization** (for debugging):
```cpp
mesh.writeFaceVTKPhysical("faces.vtk");  // Shows normals as vectors
```

### 10.3 Consistency Checks

After mesh initialization, verify:
1. ✅ All cell areas > 0
2. ✅ All face lengths > 0  
3. ✅ Normals are unit vectors: `|n| = 1 ± 1e-12`
4. ✅ Periodic connectivity: `xNodes[0,j] ≈ xNodes[ni+2,j]`
5. ✅ Face orientation: normals point away from lower-index cells

---

## 11. Summary

This O-grid mesh system provides:

✅ **Topologically consistent** structured grid for airfoils  
✅ **Ghost layers** with proper geometric extrapolation  
✅ **Shared-face storage** for efficient flux computation  
✅ **Unit normals** with consistent orientation  
✅ **Spectral radius** for CFL time-stepping  
✅ **Non-dimensional coordinates** for numerical stability  

**The mesh is fully geometric** - flow physics (boundary conditions, flux calculation) are handled by the Euler solver using this geometric foundation.

---

## 12. Quick Reference: Common Operations

```cpp
// Get physical cell bounds
int imin, imax, jmin, jmax;
mesh.getInteriorBounds(imin, imax, jmin, jmax);

// Loop over physical cells
for (int j = jmin; j < jmax; ++j) {
    for (int i = imin; i < imax; ++i) {
        int c = mesh.cellIndex(i, j);
        double xc = mesh.xCells[c];
        double yc = mesh.yCells[c];
        double A  = mesh.cellArea[c];
    }
}

// Get face geometry for flux
int f = mesh.iFaceIndex(i, j);
double nx = mesh.iFaceNormal[f][0];
double ny = mesh.iFaceNormal[f][1];
double L  = mesh.iFaceLen[f];

// Access spectral radius for time-step
double dt_cell = CFL * mesh.cellArea[c] / mesh.lambda[c];
```

---

**End of Documentation**
