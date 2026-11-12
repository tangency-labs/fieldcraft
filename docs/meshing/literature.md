# Volume meshing literature survey for FieldCraft (tet → hex)

Compiled 2026-10-13. Scope: watertight triangulated B-rep tessellation in, then (1) tetrahedral mesh by
Delaunay methods with Voronoi/ODT improvement, then (2) all-hex and hex-dominant meshes.

## How this was verified

- **DOIs.** Every DOI below resolved through the DOI handle API (`doi.org/api/handles/<doi>` → responseCode 1).
  Each title, first author, year and venue was also checked against Crossref metadata. One wrong DOI was
  caught and fixed: `10.1145/276884.276889` is Amenta's *Voronoi filtering*, and Shewchuk's SoCG'98 paper is
  `.276894`.
- **Other URLs** (author PDFs, arXiv, project pages, repos) returned HTTP 200 on 2026-10-13. Two exceptions:
  `cs.nyu.edu/~yixinhu/tetwild.pdf` (403) and `zenodo.org/record/6642020` (403 to a scripted client, though the
  Zenodo REST API returned the record). Neither is used as a primary link; I cite the DOI instead.
- **Section numbers** were read from the PDFs for: Shewchuk 1997, 1998, 2002; Klingner–Shewchuk; Labelle–Shewchuk;
  Diazzi 2021 and 2023; fTetWild; Attene 2020; Livesu 2021; Garanzha 2021; Pietroni survey; HexMe; DualCube. For any
  other paper I give no section number rather than guess one.
- **Licences** come from GitHub/GitLab licence metadata, LICENSE files or README licence statements. I read no
  source code. "none" means the repo has no licence, which makes it all-rights-reserved: treat it as unreadable
  for a clean-room project, same as GPL.

### Licence legend for associated code

| Status | Meaning |
|---|---|
| **READ OK** | Permissive (MIT, BSD, Apache-2.0, ISC, MPL-2.0 file-level, public domain). May be read for reference. Still re-implement from the paper. |
| **FORBIDDEN** | GPL, LGPL or AGPL. Do not open. |
| **NO LICENCE** | All rights reserved. Do not open. |

---

## Part A: Tetrahedral meshing

### A0. Robust predicates and exactness (foundation for everything)

**A0.1 Shewchuk 1997, adaptive-precision predicates.** **Verdict: build on.**
- **Citation.** J. R. Shewchuk. *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.* Discrete & Computational Geometry 18(3):305–363, 1997.
- **Links.** DOI [10.1007/PL00009321](https://doi.org/10.1007/PL00009321). PDF: https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf
- **Method.** Floating-point *expansions* (non-overlapping sums of doubles) with error-free transforms:
  - Two-Sum and Fast-Two-Sum (§2.3).
  - Expansion-Sum / Fast-Expansion-Sum (§2.4).
  - Two-Product with Dekker splitting (§2.5).
  - Scale-Expansion (§2.6).
  - Compression (§2.7).
  - Adaptive evaluation (§3): compute the determinant in plain doubles with a certified error bound. Escalate through
    increasingly exact stages only when the sign is uncertain.
- **Predicates.** `orient2d`, `orient3d`, `incircle` and `insphere`, each with stage-wise error bounds (§4.2).
  §5 lists caveats: overflow, underflow, and the assumption of round-to-even IEEE arithmetic.
- **Problem solved.** Exact sign of the core Delaunay predicates at near-double speed.
- **Robustness.** Exact for inputs in double precision. It does not handle constructed points such as
  intersections or Steiner points; that is what A0.3 is for. It does not resolve true degeneracies (cospherical
  points); that is what A0.2 is for.
- **Complexity.** Low to moderate: about 1–2k lines. The error-bound constants must be derived exactly as in the paper.
- **Code.**
  - `predicates.c` at https://www.cs.cmu.edu/~quake/robust.html is **public domain** (the page says so explicitly):
    READ OK. The same page notes that *Triangle* is copyrighted and not public domain; don't read it.
  - Rust port `georust/robust` (Apache-2.0, https://github.com/georust/robust): READ OK. Useful as a differential
    test oracle.
- **Why.** Every downstream algorithm (Bowyer–Watson, CDT, refinement, octree/surface classification) assumes
  exact `orient3d` and `insphere`.

**A0.2 Symbolic perturbation (Simulation of Simplicity).** **Verdict: build on.**
- **Citations.**
  - H. Edelsbrunner, E. P. Mücke. *Simulation of Simplicity: a technique to cope with degenerate cases in geometric
    algorithms.* ACM TOG 9(1):66–104, 1990. DOI [10.1145/77635.77639](https://doi.org/10.1145/77635.77639).
  - O. Devillers, M. Teillaud. *Perturbations for Delaunay and weighted Delaunay 3D triangulations.* CGTA 44(3):160–168,
    2011. DOI [10.1016/j.comgeo.2010.09.010](https://doi.org/10.1016/j.comgeo.2010.09.010).
- **Method.** When `insphere` is exactly zero, break the tie with a consistent symbolic perturbation. Devillers–Teillaud
  give a perturbation for `insphere` that needs only `orient3d` evaluations on sorted point indices. That is much
  simpler than general SoS.
- **Problem solved.** Cospherical and coplanar inputs, which are ubiquitous in CAD tessellations (grids, cylinders,
  regular polygons). The Delaunay triangulation stays unique and the Bowyer–Watson cavity stays a star-shaped ball.
- **Complexity.** Low, given A0.1.
- **Code.** No code was read.

**A0.3 Attene 2020, indirect predicates.** **Verdict: build on (from the paper only).**
- **Citation.** M. Attene. *Indirect Predicates for Geometric Constructions.* Computer-Aided Design 126:102856, 2020.
- **Links.** DOI [10.1016/j.cad.2020.102856](https://doi.org/10.1016/j.cad.2020.102856). arXiv: https://arxiv.org/abs/2105.09772
- **Method.** Represent constructed points implicitly, never rounding them:
  - LPI: line–plane intersection.
  - TPI: three-plane intersection.
  - LNC: point on a segment.
- **Evaluation.** Predicates on implicit points are expanded symbolically into homogeneous polynomials. They are then
  evaluated with a three-level filter (§5):
  - static/semi-static floating-point filter (§5.1);
  - interval arithmetic (§5.2);
  - exact expansions (§5.3).
- The paper also covers caching (§5.4) and generating predicate code automatically from a description (§6).
- **Problem solved.** Exact predicates on Steiner and intersection points without rational or bignum arithmetic.
  This is the enabling trick of A2.3 and A2.4.
- **Robustness.** Exact. The cost grows with the "depth" of the construction, so the method is designed for one
  level of implicit points defined by explicit inputs.
- **Complexity.** Moderate: about 3–5k lines, mostly generated. Write a small Rust code generator, as the paper does.
- **Code.** `MarcoAttene/Indirect_Predicates` is **LGPL-2.1**: FORBIDDEN.
- **Reference alternative.** Lévy's PCK (next entry) is BSD.

**A0.4 Lévy 2016, Predicate Construction Kit (PCK).** **Verdict: reference.**
- **Citation.** B. Lévy. *Robustness and efficiency of geometric programs: The Predicate Construction Kit (PCK).*
  Computer-Aided Design 72:3–12, 2016. DOI [10.1016/j.cad.2015.10.004](https://doi.org/10.1016/j.cad.2015.10.004).
- **Method.** A code generator from a predicate formula to:
  - a semi-static filter;
  - an exact expansion fallback;
  - SoS for degenerate cases.
- It is used in geogram's Delaunay and restricted-Voronoi code.
- **Code.** `BrunoLevy/geogram` is **BSD-3-Clause** (LICENSE file): READ OK, with one caveat. geogram vendors
  `third_party/tetgen` (AGPL) and `third_party/triangle` (non-free), and has wrappers `delaunay/delaunay_tetgen.*` and
  `delaunay/delaunay_triangle.*`. **Never open those.** Geogram's own `delaunay_3d.cpp`, `parallel_delaunay_3d.cpp`,
  `voronoi/RVD*.cpp`, `CVT.cpp` and its numerics/PCK files are BSD.
- **Verdict reason.** It is the best permissively licensed reference for a production Delaunay + RVD + predicate stack.

---

**A0.5 Boundary surface preparation (implemented in `fieldcraft-mesh::surface`).** **Verdict: BUILD ON.**
- **Guigue and Devillers 2003.** P. Guigue, O. Devillers. *Fast and Robust Triangle–Triangle Overlap Test Using
  Orientation Predicates.* Journal of Graphics Tools 8(1):25–32, 2003. Triangle overlap decided from `orient3d` /
  `orient2d` signs alone, so it is exact on the A0.1 predicates. We use the equivalent edge–triangle reduction (two
  closed triangles meet exactly when an edge of one meets the other) and add the cases of triangles that share a
  vertex or an edge, which may meet only there.
- **Ericson 2005.** C. Ericson. *Real-Time Collision Detection.* Morgan Kaufmann, 2005, ch. 6. Bounding volume
  hierarchies: top-down median-split construction and the self-collision traversal that pairs the triangles to test.
- **Jacobson, Kavan, Sorkine-Hornung 2013.** *Robust Inside-Outside Segmentation using Generalized Winding Numbers.*
  ACM TOG 32(4):33, 2013. DOI [10.1145/2461912.2461916](https://doi.org/10.1145/2461912.2461916). The winding number
  of a vertex of one shell with respect to another tells nesting, so voids must face inwards.
- **van Oosterom and Strackee 1983.** *The Solid Angle of a Plane Triangle.* IEEE Trans. Biomed. Eng. 30(2):125–126,
  1983. DOI [10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207). The triangle solid angle the
  winding number sums.

---

**A0.6 Edge winding of a finished mesh (implemented in `fieldcraft-mesh::tet::winding`).** **Verdict: BUILD ON.**
- **Hormann and Agathos 2001.** K. Hormann, A. Agathos. *The Point in Polygon Problem for Arbitrary Polygons.*
  Computational Geometry 20(3):131–144, 2001. DOI [10.1016/S0925-7721(01)00012-8](https://doi.org/10.1016/S0925-7721(01)00012-8).
  The winding number of a closed polygon about a point is the signed count of its edges that cross a ray from it, each
  edge counted on a half-open interval. Round a tetrahedral mesh's edge, the link of positive tetrahedra is such a
  polygon (each tetrahedron one step of less than half a turn); the gate counts the steps that cross the half-plane
  through one link vertex with exact `orient3d` signs (A0.1): one for an interior edge, none for a boundary edge.

---

### A1. Delaunay construction

**A1.1 Bowyer 1981 / Watson 1981, incremental cavity insertion.** **Verdict: build on.**
- **Citations.**
  - A. Bowyer. *Computing Dirichlet tessellations.* The Computer Journal 24(2):162–166, 1981.
    DOI [10.1093/comjnl/24.2.162](https://doi.org/10.1093/comjnl/24.2.162).
  - D. F. Watson. *Computing the n-dimensional Delaunay tessellation with application to Voronoi polytopes.*
    The Computer Journal 24(2):167–172, 1981. DOI [10.1093/comjnl/24.2.167](https://doi.org/10.1093/comjnl/24.2.167).
- **Method.** Insert a point, find all tetrahedra whose circumsphere contains it (the cavity), delete them, and star the
  cavity boundary to the new point.
- **Robustness.** With exact `insphere`/`orient3d` (A0.1) and SoS (A0.2) the cavity is always a star-shaped ball and
  the method never fails. With floating-point predicates it corrupts the mesh on degenerate CAD input.
- **Complexity.** Low for the core. Practical speed needs:
  - BRIO insertion order: N. Amenta, S. Choi, G. Rote, *Incremental constructions con BRIO*, SoCG 2003, 211–219,
    DOI [10.1145/777792.777824](https://doi.org/10.1145/777792.777824). Use it together with Hilbert/Morton sorting.
  - Visibility or straight walking for point location: O. Devillers, S. Pion, M. Teillaud, *Walking in a
    triangulation*, IJFCS 13(2):181–199, 2002, DOI [10.1142/S0129054102001047](https://doi.org/10.1142/S0129054102001047).
  - A compact tet-adjacency array.
  - Termination of the visibility walk in a Delaunay triangulation: H. Edelsbrunner, *An acyclicity theorem for cell
    complexes in d dimensions*, Combinatorica 10(3):251–260, 1990, DOI
    [10.1007/BF02122779](https://doi.org/10.1007/BF02122779).
  - Hilbert index: J. Skilling, *Programming the Hilbert curve*, AIP Conf. Proc. 707:381–387, 2004, DOI
    [10.1063/1.1751381](https://doi.org/10.1063/1.1751381).
- **Code.** For a permissive reference, see geogram `delaunay_3d.cpp` (BSD).
- **Implemented** in `fieldcraft-mesh::delaunay` (ghost vertex for the hull, Devillers–Teillaud perturbation for
  cospherical points, exact verifier).

**A1.2 Marot, Pellerin, Remacle 2019, "One machine, one minute, three billion tetrahedra".** **Verdict: reference (later).**
- **Citation.** IJNME 117(9):967–990, 2019. DOI [10.1002/nme.5987](https://doi.org/10.1002/nme.5987).
- **Method.** Parallel Delaunay kernel: domain partitioning on a space-filling curve, lock-free cavity work,
  cache-friendly data layout.
- **Use.** Only once single-threaded correctness is nailed down.
- **Code.** HXT lives inside Gmsh (GPL): FORBIDDEN.

---

### A2. Boundary recovery: constrained or conforming Delaunay

**A2.1 Shewchuk 2002, CDTs and provably good boundary recovery.** **Verdict: build on (theory).**
- **Citation.** J. R. Shewchuk. *Constrained Delaunay Tetrahedralizations and Provably Good Boundary Recovery.*
  Proc. 11th International Meshing Roundtable, 2002.
- **Links.** No DOI exists. PDF: https://people.eecs.berkeley.edu/~jrs/papers/cdtbasic.pdf
- **Method.** A PLC does not always have a CDT. If every segment is *strongly Delaunay*, a CDT exists. The paper
  protects edges by splitting them (§3 Edge Protection, §5 Provably Good Edge Protection) until this holds. Then it
  builds the CDT (§6), for example by incremental facet insertion (§7).
- **Problem solved.** The existence theory behind every practical CDT code.
- **Robustness.** Proven termination and existence. The provable protection inserts more Steiner points than needed.
- **Complexity.** Medium.
- **Companion papers.**
  - Shewchuk, *General-dimensional constrained Delaunay and constrained regular triangulations I*, DCG 39:580–637,
    2008. DOI [10.1007/s00454-008-9060-3](https://doi.org/10.1007/s00454-008-9060-3).
    PDF: https://people.eecs.berkeley.edu/~jrs/papers/cdtj1.pdf
  - Shewchuk, *Updating and constructing CDTs by flips* (SoCG 2003). PDF: https://people.eecs.berkeley.edu/~jrs/papers/cdtflip.pdf

**A2.2 Si and Gärtner 2005; Si 2010; Si 2015; Shewchuk and Si 2014.** **Verdict: reference (papers only; TetGen is AGPL).**
- **Citations.**
  - H. Si, K. Gärtner. *Meshing Piecewise Linear Complexes by Constrained Delaunay Tetrahedralizations.* IMR 14,
    147–163, 2005. DOI [10.1007/3-540-29090-7_9](https://doi.org/10.1007/3-540-29090-7_9).
  - H. Si. *Constrained Delaunay tetrahedral mesh generation and refinement.* Finite Elements in Analysis and Design
    46(1–2):33–46, 2010. DOI [10.1016/j.finel.2009.06.017](https://doi.org/10.1016/j.finel.2009.06.017).
  - H. Si. *TetGen, a Delaunay-based quality tetrahedral mesh generator.* ACM TOMS 41(2):11, 2015.
    DOI [10.1145/2629697](https://doi.org/10.1145/2629697).
  - J. R. Shewchuk, H. Si. *Incrementally constructing and updating constrained Delaunay tetrahedralizations with
    finite-precision coordinates.* Engineering with Computers 30:253–269, 2014.
    DOI [10.1007/s00366-013-0331-0](https://doi.org/10.1007/s00366-013-0331-0).
    PDF (IMR'12 version): https://people.eecs.berkeley.edu/~jrs/papers/cdtupdate.pdf
- **Method.**
  1. Delaunay tetrahedralization of the vertices.
  2. Segment recovery, with Steiner points placed using a protecting-ball rule that needs far fewer points than
     Shewchuk 2002.
  3. Face recovery by cavity retriangulation.
  4. Optional quality refinement (radius–edge bound) under the CDT.
- **Limitations.**
  - Steiner points are rounded to doubles, which needs tolerances and can fail.
  - Diazzi et al. 2023 §3.5 ("Possible failures in tetgen") also show a **theoretical** gap: cavity expansion can fail
    even in exact arithmetic (their Fig. 4).
- **Code.** TetGen is **AGPL**: FORBIDDEN. Read the papers only.

**A2.3 Diazzi, Panozzo, Vaxman, Attene 2023, robust CDT.** **Verdict: BUILD ON. This is the recommended boundary-recovery algorithm.**
- **Citation.** L. Diazzi, D. Panozzo, A. Vaxman, M. Attene. *Constrained Delaunay Tetrahedrization: A Robust and
  Practical Approach.* ACM TOG 42(6):181 (SIGGRAPH Asia), 2023.
- **Links.** DOI [10.1145/3618352](https://doi.org/10.1145/3618352). arXiv: https://arxiv.org/abs/2309.09805
- **Method.** Si–Gärtner's pipeline with every Steiner point represented **implicitly** with indirect predicates (A0.3),
  so no geometric tolerance is ever needed.
  - **Background, §3:**
    - §3.1 problem statement;
    - §3.2 PLC characteristics;
    - §3.3 segment recovery;
    - §3.4 face recovery;
    - §3.5 TetGen failure modes.
  - **Algorithm, §4:**
    - §4.1 why exact rational numbers are not enough;
    - §4.2 rational Steiner points;
    - §4.3 implicit Steiner points (§4.3.1 indirect `orient3d`/`inSphere`);
    - §4.4 a modified gift-wrapping algorithm that replaces Si's cavity expansion where that can fail;
    - §4.5 interior/exterior labelling;
    - §4.6 the full algorithm.
  - **Rounding, §5.1:** output can be rounded to floating-point coordinates afterwards.
- **Results.** **100% success on the 4,408 valid Thingi10K models.** The output conforms exactly to the input. It is
  parameter-free.
- **Problem solved.** Exact boundary conformity, which CAE needs for face groups, contact pairs, shared interfaces
  between parts and boundary-condition tagging, together with a termination guarantee.
- **Limitations.**
  - The input must be a valid PLC (no self-intersections). Fine for a B-rep tessellator that guarantees
    watertight, non-self-intersecting output. Otherwise run A2.4 or A2.5 first.
  - Rounding to doubles "almost always" stays valid but is not guaranteed. Keep implicit points through the quality
    phase, or round with a validity check.
  - Steiner counts can be large on nasty inputs. The paper notes about 9% of models needing over 50k.
- **Complexity.** High: about 8–12k lines of Rust including the predicate generator. It is the single most valuable
  investment in the tet pipeline.
- **Code.** `MarcoAttene/CDT` is **GPL, or LGPL by build option** (README): FORBIDDEN.
- **Implemented** in `fieldcraft-mesh::cdt` from the paper only, with these choices:
  - Steiner points are LNC points (§4.3) on the original segment, `a + t (b − a)` with `t` a double; the in-sphere test
    on them is exact (interval filter, then exact dyadic integers) and perturbed with the same lexicographic lift as
    the Delaunay step, so Bowyer–Watson insertion of Steiner points stays consistent with it.
  - Facet recovery uses the modified gift-wrapping (§4.4) for every cavity, not only after Si's cavity expansion
    fails: it is the always-correct path, and Si's expansion is a speed-up. Condition (ii) is an exact
    proper-intersection test of closed simplices (separating plane through three of the vertices, then recursion
    into the contact plane), which handles the coplanar configurations of CAD tessellations; recovered facets inside
    a cavity stay constraints.
  - Reference points (§3.3) are searched among the tetrahedra whose boxes meet the missing subsegment's box (they
    contain the tetrahedra the subsegment crosses); on the corpus this picks the same point as searching the whole
    diametral ball.
  - Rounding (§5.1): the flips of §5.1 run on the rounded coordinates (`cdt::repair`), decided by exact `orient3d` on
    the doubles: every tetrahedron that rounding (or nearly cospherical input) leaves flat, inverted or below
    FieldCraft's validity criterion is targeted by edge removal and multi-face removal that never touch a segment, and
    a flat tetrahedron lying on one input triangle is deleted. What survives (a flat tetrahedron spanning two nearly
    coplanar input triangles: 1–65 per part on 7 corpus parts) is counted by the gate and removed by the quality stage,
    whose boundary is the B-rep faces rather than the input triangles.

**A2.4 Diazzi and Attene 2021, convex polyhedral meshing.** **Verdict: reference.**
- **Citation.** L. Diazzi, M. Attene. *Convex polyhedral meshing for robust solid modeling.* ACM TOG 40(6):259, 2021.
- **Links.** DOI [10.1145/3478513.3480564](https://doi.org/10.1145/3478513.3480564). arXiv: https://arxiv.org/abs/2109.14434
- **Method.**
  1. Delaunay of the vertices.
  2. Map tets to the constraints they intersect (§3.3).
  3. Binary space partition of the tets by constraint planes, with implicit points (§3.4).
  4. Virtual constraints (§3.5).
  5. Identify constrained facets (§3.6).
  6. Partition cells into inside/outside (§3.7).
- **Output.** A convex polyhedral mesh that conforms exactly even to **self-intersecting** input. It can be
  tetrahedralized afterwards.
- **Use for us.** A repair path for dirty B-reps: overlapping faces, booleans, imported STEP/STL with intersections.
- **Code.** `MarcoAttene/VolumeMesher` is **GPL-3**: FORBIDDEN.

**A2.5 Hu et al. 2018, TetWild; Hu et al. 2020, fTetWild (envelope approach).** **Verdict: reference; use as a fallback design for dirty input.**
- **Citations.**
  - Y. Hu, Q. Zhou, X. Gao, A. Jacobson, D. Zorin, D. Panozzo. *Tetrahedral Meshing in the Wild.* ACM TOG 37(4):60,
    2018. DOI [10.1145/3197517.3201353](https://doi.org/10.1145/3197517.3201353).
  - Y. Hu, T. Schneider, B. Wang, D. Zorin, D. Panozzo. *Fast Tetrahedral Meshing in the Wild.* ACM TOG 39(4):117,
    2020. DOI [10.1145/3386569.3392385](https://doi.org/10.1145/3386569.3392385). arXiv: https://arxiv.org/abs/1908.03581
- **Method.** The input surface may move within an ε-envelope:
  - §3.2 envelope test;
  - §3.3 preprocessing;
  - §3.4 incremental triangle insertion into a floating-point tet mesh. Each insertion is rejected if it would
    create an inverted tet;
  - §3.5 AMIPS-driven split/collapse/swap/smooth improvement that never leaves the envelope.
- **Results.** fTetWild guarantees a **valid floating-point** output (TetWild produced rational coordinates that may
  invert after rounding). It meshes 98.7% of Thingi10K in under 2 minutes.
- **Limitations.**
  - fTetWild does **not** guarantee every input triangle is in the output.
  - The output surface only approximates the input. Sharp features and face IDs must be re-projected.
  - That is acceptable for "meshes in the wild" but awkward for CAE, where face groups must match B-rep faces exactly.
- **Complexity.** High.
- **Code.**
  - `wildmeshing/fTetWild` is **MPL-2.0**: READ OK (its own `src/`; it pulls in geogram and libigl).
  - Original `Yixin-Hu/TetWild` is **GPL-3.0** per GitHub metadata: FORBIDDEN.
  - `wildmeshing/wildmeshing-toolkit` is MIT: READ OK.

**A2.6 Cherchi, Livesu, Scateni, Attene 2020, mesh arrangements.** **Verdict: reference (needed only for self-intersecting input).**
- **Citation.** *Fast and robust mesh arrangements using floating-point arithmetic.* ACM TOG 39(6):250, 2020.
  DOI [10.1145/3414685.3417818](https://doi.org/10.1145/3414685.3417818).
- **Method.** Resolves all self-intersections of a triangle soup with indirect predicates.
- **Code.** `mlivesu/InteractiveAndRobustMeshBooleans` is MIT per GitHub, a related later work. **Caution:** check
  whether it bundles Attene's LGPL predicate headers before opening any file.

---

### A3. Quality: Delaunay refinement, ODT/CVT, sliver removal

**A3.1 Shewchuk 1998, tetrahedral Delaunay refinement; thesis 1997.** **Verdict: build on.**
- **Citation.** J. R. Shewchuk. *Tetrahedral Mesh Generation by Delaunay Refinement.* SoCG 1998, 86–95.
- **Links.**
  - DOI [10.1145/276884.276894](https://doi.org/10.1145/276884.276894).
  - PDF: https://people.eecs.berkeley.edu/~jrs/papers/delref3d.pdf
  - Thesis: https://www.cs.cmu.edu/~quake-papers/delaunay-refinement.pdf
- **Method.** Insert the circumcentres of tets whose radius–edge ratio exceeds B. Encroached subsegments are split
  first, then encroached subfacets.
- **Guarantee.** Termination with B = 2 when input angles are at least 90°.
- **Limitations.** It bounds radius–edge ratio, not dihedral angles. **Slivers survive**, so pair it with A3.3–A3.6.
  Small input angles, common in CAD, need the CDT-based variants: Si 2010 (A2.2), and Shewchuk & Si, SoCG 2014,
  DOI [10.1145/2582112.2582138](https://doi.org/10.1145/2582112.2582138).
- **Complexity.** Medium, on top of A2.3.
- **Code.** None read.

- **Implemented** in `fieldcraft-mesh::tet::refine`: circumcentres of tetrahedra whose radius–edge ratio exceeds B = 2
  (when their circumradius is at least half the local size) or whose circumradius exceeds 0.75 of the local size; a
  circumcentre outside the domain or inside a boundary face's diametral sphere splits that face instead (at its
  circumcentre when inside it, else at its longest edge's midpoint, a feature edge split when that edge is one), every
  boundary point projected onto its B-rep face or curve. Small input angles are handled by bounds rather than Si's
  theory: faces below half the local size are not split, no vertex closer than 0.05 of the size to another, a
  tetrahedron whose point could not go in is not retried, caps on insertions and tetrahedra. Before refinement, edge
  collapse coarsens the boundary recovery's dense Steiner points (X. Li, M. S. Shephard, M. W. Beall, *3D anisotropic
  mesh adaptation by mesh modification*, CMAME 194(48–49):4915–4950, 2005). A collapse keeps the link condition
  (T. K. Dey, H. Edelsbrunner, S. Guha, D. V. Nekhayev, *Topology preserving edge contraction*, Publ. Inst. Math.
  (Beograd) (N.S.) 66(80):23–45, 1999), the boundary coned to one extra vertex as in the paper, so it never lays an edge
  or face onto one the mesh has.

**A3.2 Cheng, Dey, Shewchuk, *Delaunay Mesh Generation* (CRC, 2012).** **Verdict: reference (the textbook).**
- **Links.** DOI [10.1201/b12987](https://doi.org/10.1201/b12987).
- **Use.** The consolidated treatment of 3D refinement, small-angle handling and sliver theory. Buy it.
  Its chapters on smooth surfaces and piecewise-smooth complexes insert ridge (feature-curve) vertices on the exact
  curve. `fieldcraft-mesh::tet` follows that for its feature curves (A3.2a).

**A3.2a Boivin and Ollivier-Gooch 2002, refinement on curved boundaries.** **Verdict: build on (curve splits).**
- **Citation.** C. Boivin, C. Ollivier-Gooch, *Guaranteed-quality triangular mesh generation for domains with curved
  boundaries*, IJNME 55(10):1185–1213, 2002. DOI [10.1002/nme.542](https://doi.org/10.1002/nme.542).
- **Method.** Ruppert-style Delaunay refinement whose boundary is curved: an encroached boundary segment is split at a
  point of the true curve, not of its chord, and the encroachment and quality guarantees are kept for curves
  discretized finely enough.
- **Use in FieldCraft.** `fieldcraft-mesh::tet` (refinement, ODT smoothing and improvement) places every feature-curve
  vertex on the exact curve from its insertion on, wherever both faces' exact surfaces are known (#674). Placed on the
  tessellation's chord, a curve refined below the facet size left its vertices up to the chord's sag off the curve, with
  the faces beside them meshed around that, and the final snap could not move them out without flattening elements.
  Face vertices stay on the tessellation until the final snap: projected onto the exact surfaces during improvement,
  vertices on both sides of a blend's tangent edge come out too close to flat to separate.

**A3.3 Alliez, Cohen-Steiner, Yvinec, Desbrun 2005, Variational Tetrahedral Meshing.** **Verdict: build on (ODT smoothing).**
- **Citation.** ACM TOG 24(3):617–625, 2005. DOI [10.1145/1073204.1073238](https://doi.org/10.1145/1073204.1073238).
- **Method.** Minimize the ODT energy: the L1 error between the paraboloid and its piecewise-linear interpolant.
  Alternate between:
  - moving each interior vertex to the volume-weighted average of its incident tets' circumcentres (closed form);
  - re-Delaunaying the connectivity.
- **Findings.** ODT outperforms CVT for tet quality. Results improve with a sizing field.
- **Limitations.**
  - Boundary vertices need special handling, and the paper's boundary treatment is the weak point.
  - Slivers remain at the end. The paper uses a final perturbation/exudation pass.
- **Complexity.** Low to moderate, given a Delaunay kernel.

**A3.4 Chen and Xu 2004; Chen and Holst 2011, Optimal Delaunay Triangulations.** **Verdict: build on (theory and update formulas).**
- **Citations.**
  - L. Chen, J. Xu. *Optimal Delaunay triangulations.* J. Computational Mathematics 22(2):299–308, 2004. No DOI.
    Publisher page: https://ojs.global-sci.org/index.php/JCM/article/view/11631. Author page with PDF link:
    https://www.math.uci.edu/~chenlong/CX2004.html
  - L. Chen, M. Holst. *Efficient mesh optimization schemes based on Optimal Delaunay Triangulations.* CMAME
    200(9–12):967–984, 2011. DOI [10.1016/j.cma.2010.11.007](https://doi.org/10.1016/j.cma.2010.11.007).
- **Method.** The Delaunay triangulation minimizes the interpolation error of ‖x‖² over all triangulations of a point
  set. Chen–Holst give the exact and approximate vertex-update rules (ODT smoothing), including a boundary-aware variant.
- **Complexity.** Low.

**A3.5 Tournois, Wormser, Alliez, Desbrun 2009, Delaunay refinement interleaved with optimization.** **Verdict: build on (the overall schedule).**
- **Citation.** *Interleaving Delaunay Refinement and Optimization for Practical Isotropic Tetrahedron Mesh
  Generation.* ACM TOG 28(3):75, 2009. DOI [10.1145/1531326.1531381](https://doi.org/10.1145/1531326.1531381).
- **Method.** Alternate Delaunay refinement and ODT/NODT smoothing. The refinement uses the radius–edge and boundary
  criteria. NODT is a Natural ODT variant for boundary vertices. Finish with sliver perturbation and exudation.
- **Results.** Minimum dihedral angles are substantially better than refinement alone.
- **Problem solved.** The practical recipe used by CGAL Mesh_3-style generators.
- **Complexity.** Moderate.
- **Code.** CGAL is GPL: FORBIDDEN. This paper is the readable source of the scheme.

- **Implemented** in `fieldcraft-mesh::tet::smooth` as the density-weighted ODT update (Chen–Xu, Alliez et al.; weight
  `h⁻³`), interleaved with refinement in three rounds as in A3.5; boundary vertices move in their face's tangent plane
  and are projected back onto the B-rep face, curve vertices move along their curve, corners stay. A move is kept when
  every tetrahedron stays positive and its star's worst quality stays above min(old worst, sin 17.5°). Tetrahedra whose
  circumcentre is further than their longest edge from their centroid contribute their centroid instead.

**A3.6 Klingner and Shewchuk 2007, Aggressive Tetrahedral Mesh Improvement.** **Verdict: BUILD ON (final sliver-removal stage).**
- **Citation.** B. M. Klingner, J. R. Shewchuk. IMR 16, 3–23, 2007.
- **Links.** DOI [10.1007/978-3-540-75103-8_1](https://doi.org/10.1007/978-3-540-75103-8_1).
  PDF: https://people.eecs.berkeley.edu/~jrs/papers/aggress.pdf
- **Method.**
  - §2: quality measures. Use volume–length, or the minimum sine of dihedral angles.
  - §3, operations:
    - §3.1 optimization-based smoothing with an active-set method on the worst quality;
    - §3.2 edge removal;
    - §3.3 multi-face removal;
    - §3.4 vertex insertion;
    - §3.5 composites.
  - §4.2: hill-climbing schedule that accepts only operations improving the local minimum quality.
- **Results.** Dihedral angles are usually between **31° and 149°**, or 23°–136° with another objective.
- **Problem solved.** Kills slivers that ODT and refinement leave behind. Works on any valid tet mesh.
- **Robustness.** Monotone in worst quality, so it never makes the mesh worse.
- **Limitations.** Slow. Boundary vertices are constrained.
- **Complexity.** Moderate: edge removal requires Shewchuk's dynamic programming over the edge ring.
  PDF: https://people.eecs.berkeley.edu/~jrs/papers/edge.pdf, not separately verified.
- **Code.** The Stellar code's licence was not verified, so it is not read.
- **Implemented** in `fieldcraft-mesh::tet::{ops, improve}` from the two papers: Klincsek's dynamic programming for
  edge removal (interior edges; boundary edges of one flat B-rep face end in a 2-2 flip; the triangulation never takes
  an edge the mesh has elsewhere, else a boundary vertex with three faces round it is flipped onto a face the mesh
  has), multi-face removal grown greedily over the sandwiched faces, nonsmooth optimization smoothing (active set,
  smallest vector of the convex hull of the active gradients, line search; gradients by central differences; a
  boundary vertex moves only while none of its faces turns by more than 30° or crosses another boundary face), and the
  composite vertex insertion with rollback (barycentres of boundary faces, of the tetrahedron, midpoints of boundary
  edges, and a point inside off a flat boundary edge). Deviations: the inserted vertex takes the Delaunay cavity
  rather than the optimal cut of §3.4; each pass smooths and flips the tetrahedra below sin 40° and inserts only when
  the worst did not improve; passes and the insertion list are capped. Objective: the smaller of the minimum dihedral
  sine and the scaled Jacobian at the worst corner weighted so that a Jacobian of 0.3 counts as a 30° dihedral angle:
  the sine misses needles (every dihedral angle fine, one corner nearly flat; Shewchuk 2002); the scaled Jacobian
  analysts check (Knupp 2001, B6.6) catches them, and with the weight it decides only near that limit. Every operation
  is a checked replacement of a region (positive by exact `orient3d`, faces matched with opposite orientation, removed
  vertices fully inside, no face or edge it makes held by a tetrahedron outside it, its new boundary faces meeting the
  rest of the boundary only in shared vertices by the exact triangle test of the surface checks, found through a grid
  of the boundary faces), so the mesh is valid after each step: two sheets of a thin run-out, hundredths of a
  millimetre apart, are not pushed through each other.

**A3.7 Cheng, Dey, Edelsbrunner, Facello, Teng 2000, Sliver exudation.** **Verdict: reference.**
- **Citation.** J. ACM 47(5):883–904, 2000. DOI [10.1145/355483.355487](https://doi.org/10.1145/355483.355487).
- **Method.** Assign weights to vertices so the weighted Delaunay triangulation contains no slivers. This gives a
  provable but tiny dihedral bound.
- **Use.** Klingner–Shewchuk dominates it in practice. Useful as a cheap pass inside the A3.5 schedule.

**A3.8 Labelle and Shewchuk 2007, Isosurface stuffing.** **Verdict: reference (useful only for organic, smooth domains).**
- **Citation.** ACM TOG 26(3):57, 2007.
- **Links.** DOI [10.1145/1276377.1276448](https://doi.org/10.1145/1276377.1276448).
  PDF: https://people.eecs.berkeley.edu/~jrs/papers/stuffing.pdf
- **Method.** BCC lattice plus snapping and warping with stencils (§3). It has **guaranteed** dihedral bounds (§4).
- **Limitations.** It does **not preserve sharp features**, so it is wrong for CAD.
- **Code.** `crawforddoran/quartet` is ISC: READ OK. Do not read `alecjacobson/isostuffer` (LGPL-3).

---

### A4. Voronoi-based: CVT, restricted Voronoi diagrams, Lp-CVT

**A4.1 Du, Faber, Gunzburger 1999; Du and Wang 2003.** **Verdict: reference.**
- **Citations.**
  - *Centroidal Voronoi Tessellations: Applications and Algorithms.* SIAM Review 41(4):637–676, 1999.
    DOI [10.1137/S0036144599352836](https://doi.org/10.1137/S0036144599352836).
  - Q. Du, D. Wang. *Tetrahedral mesh generation and optimization based on centroidal Voronoi tessellations.*
    IJNME 56(9):1355–1373, 2003. DOI [10.1002/nme.616](https://doi.org/10.1002/nme.616).
- **Method.** Lloyd relaxation toward a CVT. The dual Delaunay triangulation gives well-shaped tets apart from slivers.
- **Assessment.** Alliez 2005 shows ODT beats CVT for tet quality, so CVT is not our tet smoother.

**A4.2 Liu, Wang, Lévy, Sun, Yan, Lu, Yang 2009, CVT energy smoothness and quasi-Newton.** **Verdict: build on (the CVT solver).**
- **Citation.** ACM TOG 28(4):101, 2009. DOI [10.1145/1559755.1559758](https://doi.org/10.1145/1559755.1559758).
- **Method.** The CVT energy is C². Minimize it with L-BFGS instead of Lloyd for an order-of-magnitude speedup.

**A4.3 Yan, Lévy, Liu, Sun, Wang 2009 (restricted Voronoi diagram); Yan et al. 2013 (clipped Voronoi diagram).** **Verdict: build on (the RVD and clipped-VD kernel).**
- **Citations.**
  - Yan et al. *Isotropic Remeshing with Fast and Exact Computation of Restricted Voronoi Diagram.* CGF 28(5):1445–1454,
    2009. DOI [10.1111/j.1467-8659.2009.01521.x](https://doi.org/10.1111/j.1467-8659.2009.01521.x).
  - D.-M. Yan, W. Wang, B. Lévy, Y. Liu. *Efficient computation of clipped Voronoi diagram for mesh generation.*
    CAD 45(4):843–852, 2013. DOI [10.1016/j.cad.2011.09.004](https://doi.org/10.1016/j.cad.2011.09.004).
- **Method.** Compute the Voronoi diagram restricted to a surface (or clipped by a volume). Propagate over triangles
  and clip each against the bisector half-spaces of the nearest seeds (Sutherland–Hodgman).
- **Problem solved.** The integration domain for surface remeshing and for CVT/Lp-CVT inside a closed volume.
- **Robustness.** Exact with PCK-style predicates (geogram).
- **Complexity.** Moderate.
- **Code.** geogram `voronoi/RVD*.cpp` (BSD): READ OK.

**A4.4 Lévy and Liu 2010, Lp Centroidal Voronoi Tessellation.** **Verdict: reference (hex-dominant point placement, superseded for that purpose by Sokolov 2016).**
- **Citation.** B. Lévy, Y. Liu. *Lp Centroidal Voronoi Tessellation and its Applications.* ACM TOG 29(4):119, 2010.
  DOI [10.1145/1778765.1778856](https://doi.org/10.1145/1778765.1778856).
- **Method.** Replace the L2 CVT energy with an Lp norm in a local frame (anisotropy tensor). As p grows, Voronoi cells
  become cubes aligned to a field. Feature constraints come from integrating over the RVD of the surface and its
  sharp edges.
- **Use.** Produces cube-like point sets for hex-dominant meshing.
- **Limitations.** Non-convex energy. Needs a smooth frame field.

---

## Part B: Hexahedral meshing

### B0. Survey and benchmark (read first)

**B0.1 Pietroni et al. 2022, survey.** **Verdict: reference (required reading).**
- **Citation.** N. Pietroni, M. Campen, A. Sheffer, G. Cherchi, D. Bommes, X. Gao, R. Scateni, F. Ledoux,
  J.-F. Remacle, M. Livesu. *Hex-Mesh Generation and Processing: a Survey.* ACM TOG 42(2):16, 2022.
- **Links.** DOI [10.1145/3554920](https://doi.org/10.1145/3554920). arXiv: https://arxiv.org/abs/2202.12670
- **Structure.**
  - §2.5 integer-grid maps;
  - §4.6 grid-based;
  - §4.7 polycube;
  - §4.8 frame fields;
  - §4.9 hex-dominant;
  - §5 topological operators;
  - §6 geometric optimization;
  - §8 open problems;
  - §9 resources.
- **Key quote (§4.6).** Grid-based methods "really stand out in terms of robustness. To date, they are the only fully
  automatic methods capable of successfully hex-meshing any input shape... they are the only automatic methods
  currently implemented in professional software." It also notes their cost: orientation dependence, intricate
  high-valence singular structure, and no coarse block structure.
- **§8.1.** Frame-field methods are "limited to rather simple shapes" fully automatically. No characterization of
  hex-meshable fields exists. Polycube graph characterization is limited to 3-connected graphs and, in the Eppstein–Mumford
  form, genus 0.

**B0.2 Beaufort, Reberol, Kalmykov, Liu, Ledoux, Bommes 2022, Hex Me If You Can (HexMe).** **Verdict: reference (benchmark, with a provenance caveat in C).**
- **Citation.** CGF 41(5):125–134 (SGP), 2022.
- **Links.**
  - DOI [10.1111/cgf.14608](https://doi.org/10.1111/cgf.14608).
  - arXiv: https://arxiv.org/abs/2111.10295
  - Dataset DOI [10.5281/zenodo.6642020](https://doi.org/10.5281/zenodo.6642020), **CC-BY-4.0**.
  - Catalog: https://hexme.algohex.eu
- **Contents.** 189 tet meshes from 63 CAD models in three categories: simple, nasty and industrial. Each comes with
  feature tags.
- **Baseline result.** A state-of-the-art frame-field → IGM → HexEx pipeline meshes all features correctly on only
  **19 of 189**. On nasty and industrial shapes about half the meshes preserve under 1–8% of surfaces.
- **Why it matters.** This is the evidence that frame-field methods are not production-ready for CAD.

---

### B1. Grid / octree-based (robust family)

**B1.1 Schneiders 1996.** **Verdict: reference (historical).**
- **Citation.** R. Schneiders. *A grid-based algorithm for the generation of hexahedral element meshes.* Engineering
  with Computers 12:168–177, 1996. DOI [10.1007/BF01198732](https://doi.org/10.1007/BF01198732).
- **Method.** An octree with templates resolves hanging nodes (3-refinement). The boundary layer is projected onto the
  surface.

**B1.2 Zhang and Bajaj 2006.** **Verdict: reference.**
- **Citation.** *Adaptive and quality quadrilateral/hexahedral meshing from volumetric data.* CMAME 195:942–960, 2006.
  DOI [10.1016/j.cma.2005.02.016](https://doi.org/10.1016/j.cma.2005.02.016).
- **Method.** Dual contouring on an octree with templates, followed by quality improvement.

**B1.3 Maréchal 2009, octree-based all-hex with sharp features.** **Verdict: reference (feature handling ideas).**
- **Citation.** L. Maréchal. *Advances in Octree-Based All-Hexahedral Mesh Generation: Handling Sharp Features.*
  IMR 18, 65–84, 2009. DOI [10.1007/978-3-642-04319-2_5](https://doi.org/10.1007/978-3-642-04319-2_5).
- **Method.** A 2-refinement octree. Hanging nodes are removed by **dualization** with flat transition schemes. Sharp
  features are captured by snapping dual vertices and by padding.
- **Significance.** The basis of commercial octree hex meshers.
- **Limitation.** The scheme set is incomplete: Livesu 2021 shows it misses some transitions.

**B1.4 Gao, Shen, Panozzo 2019, Feature Preserving Octree-Based Hexahedral Meshing.** **Verdict: reference (pipeline design).**
- **Citation.** CGF 38(5):135–149 (SGP), 2019. DOI [10.1111/cgf.13795](https://doi.org/10.1111/cgf.13795).
- **Pipeline.**
  1. Octree refinement driven by surface approximation.
  2. Dual-based hex extraction.
  3. Feature-preserving surface projection.
  4. Padding.
  5. Global optimization that maximizes the minimum scaled Jacobian while fitting features.
- **Results.** The authors report high robustness on a large dataset.
- **Code.** `gaoxifeng/Feature-Preserving-Octree-Hex-Meshing` has **NO LICENCE**: do not read.

**B1.5 Livesu, Pitzalis, Cherchi 2021, Optimal Dual Schemes for Adaptive Grid Based Hexmeshing.** **Verdict: BUILD ON (the core of the first hex mesher).**
- **Citation.** ACM TOG 41(2):15, 2022 (online 2021).
- **Links.** DOI [10.1145/3494456](https://doi.org/10.1145/3494456). arXiv: https://arxiv.org/abs/2103.07745
- **Method, part 1.**
  - Enumerates **all** transitions dual methods must handle (§4: §4.1 convex, §4.2 concave, §4.3 around corners).
  - Shows prior schemes miss some of them.
  - Gives the minimal scheme set with the simplest singular structure (§5).
- **Method, part 2.**
  - Extends support from strongly to **weakly balanced** grids (§6; discussion in §7.1). Coarser meshes for the same accuracy.
  - Compares with Gao 2019 and Maréchal 2009 (§7.2).
  - Spells out the implementation details that earlier papers left implicit.
- **Robustness.** Combinatorially complete for balanced grids, so the dual is always all-hex. Geometry
  (projection/padding/untangling) is outside the paper's scope.
- **Limitations.**
  - Most grid methods do not guarantee the output keeps the input's genus or component count (survey §4.6).
  - Sharp features need a separate projection and optimization stage.
- **Complexity.** Moderate: the schemes are tables plus a balancing pass.
- **Code.** Integrated into `mlivesu/cinolib`, **MIT** (paper p. 2 and §8): READ OK. Do not open cinolib's optional
  TetGen/Triangle wrappers if present.

**B1.6 Pitzalis, Livesu, Cherchi, Gobbetti, Scateni 2021, Generalized Adaptive Refinement.** **Verdict: build on (refinement and balancing stage).**
- **Citation.** ACM TOG 40(6):257, 2021. DOI [10.1145/3478513.3480508](https://doi.org/10.1145/3478513.3480508).
- **Method.** Generalizes balanced adaptive refinement so it works on any input hex mesh, not just a regular grid,
  with conforming templates and no dualization.
- **Use.** Lets us refine a coarse aligned hex mesh, for example polycube output, as well as an octree.
- **Code.** `cg3hci/Gen-Adapt-Ref-for-Hexmeshing` is **MIT**: READ OK.
- **Finding (2026-10).** The reference pipeline chooses transition vertices with an integer linear
  program (Gurobi) and installs CinoLib's schemes by heuristics that are retried with other
  parameters when they fail; neither was adopted (see B1.10).

**B1.7 Schneiders, Algorithms for quadrilateral and hexahedral mesh generation (VKI lecture notes,
2000), §5.** **Verdict: built on (the transition templates of `fieldcraft-mesh::octree`).**
- **Link.** https://www.robertschneiders.de/papers/vki.pdf (author's copy).
- **Method.** 3-refinement driven by nodal marks: an element is split by the template of its
  marked corners, and faces split by their own marks only, so the refinement is conforming.
- **Limit (stated there).** The 2D templates of the original set make some 3D mark patterns
  impossible: a hexahedron with seven marked corners has 51 boundary quads, an odd number, so no
  hexahedral mesh of it exists (Mitchell's parity condition). Concave refined regions therefore
  need their concave cells refined. The fix is a modified edge rule (an edge with a marked end
  splits in three) with buffer layers in the face templates.
- **Finding (#85, 2026-10).** Read in full (§5.3, Fig. 57–61). Under the modified rule the 256 mark
  patterns fall into 14 classes of split edges; the faces are plain, one-corner (5 quads),
  three-edge (template 2a, 10 quads) or full 3 × 3, so parity never fails and seven marks are the
  27-split. The six-mark cell (a concave edge) is given only as a sketch ("the two trapezoidal
  faces become a dangling face"). Working it through: its two 2a faces both fold their three
  columns into the one unsplit edge, the cell corner there has five surface quads, and in the
  sheet structure the chord through the two trapezoids must cross a transverse sheet twice, so
  every hexahedral filling folds a sheet back on itself near that edge. No construction with
  usable quality was found, so the rule was not adopted.

**B1.8 Ito, Shih, Soni 2009, Octree-based reasonable-quality hexahedral mesh generation using a
new set of refinement templates.** **Verdict: reference.**
- **Citation.** IJNME 77(13):1809–1833, 2009. Templates for concave edges of the refined region
  (paper not obtained; cited after Tong and Zhang 2025).

**B1.9 Tong, Halilaj, Zhang 2024, HybridOctree_Hex.** **Verdict: reference (buffer-zone rule adopted).**
- **Citation.** H. Tong, E. Halilaj, Y. J. Zhang, *HybridOctree_Hex: Hybrid octree-based adaptive
  all-hexahedral mesh generation with Jacobian control*, arXiv 2401.05984, 2024.
- **Method.** Dual of a strongly balanced, paired octree with five templates (a 13-hexahedron
  frustum per transition face plus four transition-edge templates); a core cleared at the
  boundary so that the normals of the quads around each core boundary vertex lie in one
  half-space, then a buffer layer projected onto the surface and optimised.
- **Code.** `CMU-CBML/HybridOctree_Hex` is **BSD-3-Clause**: READ OK (hard-coded template
  geometry, not adopted).

**B1.10 Tong, Zhang 2025, MCHex.** **Verdict: reference (next boundary stage to evaluate).**
- **Citation.** H. Tong, Y. J. Zhang, *MCHex: Marching cubes based adaptive hexahedral mesh
  generation with guaranteed positive Jacobian*, arXiv 2511.02064, 2025.
- **Method.** The same 3-refinement templates (B1.7); the boundary cells are cut by marching
  cubes on a signed distance field, every cut cell's pieces are 3-regular polyhedra, and their
  midpoint subdivision gives hexahedra with a proven positive Jacobian (Bernstein coefficients) and
  a boundary that is pillowed by construction; 17 of the 31 cube cases are refined locally first.
- **Limits (stated there).** Sharp edges and corners are lost (marching cubes smooths them), and
  the convex-only transitions over-refine. FieldCraft keeps the core-padding-pillow boundary with
  feature capture for CAD parts; MCHex's cut cells are the candidate for organic shapes.

**Implemented (`fieldcraft-mesh::octree`, `mesh.hex.octree`, 2026-10).** B1.7 templates on a
27-tree refined by curvature, thickness and narrow faces; a core cleared by the B1.9 rule; the
feature map of B1.4 (patches first); padding plus a pillow per patch (B1.4 §3.4); the B6.4 energy
with a fitting term (B1.4 §3.5); the certified check of B6.6 with local repair. The dual schemes
of B1.5 were not adopted (B1.6 finding); concave refined regions are closed by refinement (B1.7
limit), so curvature refines one level only. Since #85 the levels and root size are chosen before
the grid is built, from the thinnest wall (measured to a facing wall only: a ray that leaves
through the face beyond a convex edge measures how near the edge it started, not a wall) and the
narrowest patch; a part beyond the levels or cells allowed fails at once with the numbers.

---

### B2. Polycube maps

**B2.1 Tarini, Hormann, Cignoni, Montani 2004, PolyCube-Maps.** **Verdict: reference (the concept).**
- **Citation.** ACM TOG 23(3):853–860, 2004. DOI [10.1145/1015706.1015810](https://doi.org/10.1145/1015706.1015810).

**B2.2 Gregson, Sheffer, Zhang 2011, All-hex mesh generation via volumetric polycube deformation.** **Verdict: reference.**
- **Citation.** CGF 30(5):1407–1416, 2011.
  DOI [10.1111/j.1467-8659.2011.02015.x](https://doi.org/10.1111/j.1467-8659.2011.02015.x).
- **Method.** Rotation-driven volumetric deformation toward axis-aligned normals, then labelling, then hex sampling.
- **Limitation.** Heuristic cleanup of invalid labellings.

**B2.3 Livesu, Vining, Sheffer, Gregson, Scateni 2013, PolyCut.** **Verdict: reference.**
- **Citation.** ACM TOG 32(6):171, 2013. DOI [10.1145/2508363.2508388](https://doi.org/10.1145/2508363.2508388).
- **Links.** Project: https://www.cs.ubc.ca/labs/imager/tr/2013/polycut/
- **Method.** Graph-cut labelling with monotonicity. Applies the Eppstein–Mumford validity conditions.
- **Limitation.** Not guaranteed to converge to a valid polycube (survey §4.7).

**B2.4 Huang et al. 2014 (ℓ1 polycube); Fang et al. 2016 (closed-form induced polycube); Fu et al. 2016.** **Verdict: reference.**
- **Citations.**
  - Huang, Jiang, Shi, Tong, Bao, Desbrun. *ℓ1-Based Construction of Polycube Maps from Complex Shapes.* ACM TOG
    33(3):25, 2014. DOI [10.1145/2602141](https://doi.org/10.1145/2602141).
  - Fang, Xu, Bao, Huang. *All-hex meshing using closed-form induced polycube.* ACM TOG 35(4):124, 2016.
    DOI [10.1145/2897824.2925957](https://doi.org/10.1145/2897824.2925957).
  - Fu, Bai, Liu. *Efficient volumetric PolyCube-map construction.* CGF 35(7):97–106, 2016.
    DOI [10.1111/cgf.13007](https://doi.org/10.1111/cgf.13007).
- **Method.** Deformation-based polycube construction with different energies.
- **Limitation.** They "seldom converge to a topologically consistent polyhedron" without a valid guiding labelling
  (survey §4.7).

**B2.5 Guo, Liu, Yan, Liu 2020, Cut-enhanced PolyCube-maps.** **Verdict: reference.**
- **Citation.** *Cut-enhanced PolyCube-maps for feature-aware all-hex meshing.* ACM TOG 39(4):106, 2020.
  DOI [10.1145/3386569.3392378](https://doi.org/10.1145/3386569.3392378).
- **Method.** Interior cuts let polycubes align to features that a plain polycube cannot.

**B2.6 Dumery, Protais, Mestrallet, Bourcier, Ledoux 2022, Evocube.** **Verdict: reference.**
- **Citation.** *Evocube: a Genetic Labelling Framework for Polycube-Maps.* CGF 41(6):467–479, 2022.
- **Links.** DOI [10.1111/cgf.14649](https://doi.org/10.1111/cgf.14649). arXiv: https://arxiv.org/abs/2205.00738
- **Method.** Genetic search over labellings, scored by validity and distortion. Slow, and still with no guarantee.
- **Code.** `LIHPC-Computational-Geometry/evocube` is **GPL-3**: FORBIDDEN.

**B2.7 Protais, Reberol, Ray, Corman, Ledoux, Sokolov 2022, Robust Quantization for Polycube Maps.** **Verdict: build on (when polycube is built).**
- **Citation.** CAD 150:103321, 2022. DOI [10.1016/j.cad.2022.103321](https://doi.org/10.1016/j.cad.2022.103321).
- **Method.** Robust integer quantization of a polycube map's chart sizes. Guarantees a valid hex layout from a valid
  labelling.
- **Code.** `fprotais/robustPolycube` is **AGPL-3**: FORBIDDEN.

**B2.8 Mandad, Chen, Bommes, Campen 2022, Intrinsic Mixed-Integer Polycubes.** **Verdict: reference.**
- **Citation.** CAGD 94:102078, 2022. DOI [10.1016/j.cagd.2022.102078](https://doi.org/10.1016/j.cagd.2022.102078).

**B2.9 Snoep, Speckmann, Verbeek 2025, DualCube.** **Verdict: BUILD ON (when polycube is built). The first labelling method with a validity guarantee.**
- **Citation.** *Robust Construction of Polycube Segmentations via Dual Loops.* CGF 44(5):e70195 (SGP), 2025.
- **Links.** DOI [10.1111/cgf.70195](https://doi.org/10.1111/cgf.70195). arXiv: https://arxiv.org/abs/2402.00652
- **Method.**
  - Work in the **dual** loop structure of a polycube. §3 gives the dual structure, its characterization and the
    valid modifications.
  - §4: start from a trivial polycube of the right genus, then iteratively add loops (§4.1) while keeping validity.
    Primalize into a segmentation (§4.2) and optimize the loop structure (§4.3).
  - §5.2 covers hex meshing.
- **Guarantee.** "Guaranteed to return a valid polycube segmentation for 3D models of any genus." It is valid at every
  iteration, so it can be stopped at any time.
- **Limitations.** Quality and feature alignment trade against complexity. Feature preservation is weaker than with
  feature-driven labellings.
- **Companion.** *Polycubes via Dual Loops*, arXiv https://arxiv.org/abs/2410.16865 (2024): the dual-loop
  characterization for arbitrary genus.
- **Code.** The repo named in search results (`maximsnoep/polycube-layouts-via-iterative-dual-loops`) returned 404.
  Licence unknown.

---

### B3. Frame fields and integer-grid maps (IGM)

**B3.1 Nieser, Reitebuch, Polthier 2011, CubeCover.** **Verdict: reference (defines the IGM pipeline).**
- **Citation.** *CubeCover – Parameterization of 3D Volumes.* CGF 30(5):1397–1406, 2011.
  DOI [10.1111/j.1467-8659.2011.02014.x](https://doi.org/10.1111/j.1467-8659.2011.02014.x).
- **Method.** Frame field → seamless volumetric parametrization whose integer iso-surfaces form hexes. It is the 3D
  analogue of Bommes, Zimmer, Kobbelt, *Mixed-integer quadrangulation*, ACM TOG 28(3):77, 2009,
  DOI [10.1145/1531326.1531383](https://doi.org/10.1145/1531326.1531383).

**B3.2 Huang, Tong, Wei, Bao 2011, Boundary aligned smooth 3D cross-frame field.** **Verdict: reference.**
- **Citation.** ACM TOG 30(6):143, 2011. DOI [10.1145/2070781.2024177](https://doi.org/10.1145/2070781.2024177).
- **Method.** Spherical-harmonic (SH) frame representation, which removes the octahedral-symmetry ambiguity.

**B3.3 Li, Liu, Xu, Wang, Guo 2012, All-hex meshing using singularity-restricted field.** **Verdict: reference.**
- **Citation.** ACM TOG 31(6):177, 2012. DOI [10.1145/2366145.2366196](https://doi.org/10.1145/2366145.2366196).
- **Method.** Local operations repair singularity graphs to be hex-meshable.

**B3.4 Ray, Sokolov, Lévy 2016, Practical 3D frame field generation.** **Verdict: build on (only if the frame-field family is pursued).**
- **Citation.** ACM TOG 35(6):233, 2016. DOI [10.1145/2980179.2982408](https://doi.org/10.1145/2980179.2982408).
- **Method.** SH-based smoothing with local Euler-angle optimization and boundary alignment. Fast and simple.

**B3.5 Solomon, Vaxman, Bommes 2017, Boundary element octahedral fields.** **Verdict: reference.**
- **Citation.** ACM TOG 36(3):28, 2017. DOI [10.1145/3065254](https://doi.org/10.1145/3065254).

**B3.6 Liu, Zhang, Chien, Solomon, Bommes 2018, Singularity-constrained octahedral fields.** **Verdict: reference.**
- **Citation.** ACM TOG 37(4):93, 2018. DOI [10.1145/3197517.3201344](https://doi.org/10.1145/3197517.3201344).
- **Method.** Enumerates local meshability conditions for valence-3 and valence-5 singular edges. Computes a field
  conforming to a given singularity graph.
- **Limitation.** The singular graph is user- or heuristic-supplied.

**B3.7 Corman and Crane 2019, Symmetric Moving Frames.** **Verdict: reference.**
- **Citation.** ACM TOG 38(4):87, 2019. DOI [10.1145/3306346.3323029](https://doi.org/10.1145/3306346.3323029).
- **Links.** PDF: https://www.cs.cmu.edu/~kmcrane/Projects/SymmetricMovingFrames/paper.pdf
- **Method.** Cartan's moving frames. Frame fields with *prescribed* singularity graph and non-orthogonal frames,
  obtained as connection 1-forms with integrability constraints.

**B3.8 Palmer, Bommes, Solomon 2020, Algebraic representations for volumetric frame fields.** **Verdict: reference.**
- **Citation.** ACM TOG 39(2):16, 2020.
- **Links.** DOI [10.1145/3366786](https://doi.org/10.1145/3366786). arXiv: https://arxiv.org/abs/1908.05411
- **Method.** Octahedral and odeco frames as varieties. Manifold-optimization (MBO-style) field computation. Supports
  non-orthogonal frames.

**B3.9 Lyon, Bommes, Kobbelt 2016, HexEx; Kohler, Heistermann, Bommes 2025, HexHex.** **Verdict: reference (needed only for the frame-field family).**
- **Citations.**
  - *HexEx: Robust Hexahedral Mesh Extraction.* ACM TOG 35(4):123, 2016.
    DOI [10.1145/2897824.2925976](https://doi.org/10.1145/2897824.2925976).
  - *HexHex: Highspeed Extraction of Hexahedral Meshes.* ACM TOG 44(4), 2025.
    DOI [10.1145/3730940](https://doi.org/10.1145/3730940). Project: https://cgg.unibe.ch/publications/hexhex/
- **Method.** Extract the hex mesh from an IGM, tolerating some defects (HexEx). HexHex assumes a *locally injective*
  map and is unconditionally robust, with about 30× speedup through a propeller data structure and conservative
  rasterization.
- **Code.**
  - libHexEx is **GPLv3** (RWTH page): FORBIDDEN.
  - `cgg-bern/libHexHex` LICENSE is **GPL** (header plus `GPL.txt`): FORBIDDEN.

**B3.10 Brückler, Bommes, Campen: 3D motorcycle complex (CGF 2022), QGP3D (TOG 2022), cell-complex collapsing (TOG 2023), Integer-Sheet-Pump (CGF 2024).** **Verdict: reference.**
- **Citations.**
  - DOI [10.1111/cgf.14470](https://doi.org/10.1111/cgf.14470): *The 3D Motorcycle Complex for Structured Volume
    Decomposition*, CGF 41(2):221–235.
  - DOI [10.1145/3528223.3530123](https://doi.org/10.1145/3528223.3530123): *Volume parametrization quantization for
    hexahedral meshing*, TOG 41(4):60.
  - DOI [10.1145/3618384](https://doi.org/10.1145/3618384): *Collapsing Embedded Cell Complexes for Safer Hexahedral
    Meshing*, TOG 42(6):205.
  - DOI [10.1111/cgf.15131](https://doi.org/10.1111/cgf.15131): *Integer-Sheet-Pump Quantization*, CGF 43(5).
- **Method.** Robust quantization of seamless volume maps into valid IGMs through the motorcycle complex.
- **Code.** QGP3D, MC3D, TrulySeamless3D and C4HexMeshing are all **GPL-3**: FORBIDDEN.

**B3.11 Liu and Bommes 2023, Locally Meshable Frame Fields.** **Verdict: reference (best current evidence of the family's state).**
- **Citation.** ACM TOG 42(4), 2023. DOI [10.1145/3592457](https://doi.org/10.1145/3592457).
- **Result.** Raises the valid-IGM success rate on HexMe "**from 2% to 58%**" (project page abstract). Local
  meshability is necessary but not sufficient.
- **Code.** No code link on the project page.

**B3.12 Nigolian, Campen, Bommes 2023, Expansion Cones.** **Verdict: reference (injective volumetric maps).**
- **Citation.** ACM TOG 42(4), 2023. DOI [10.1145/3592421](https://doi.org/10.1145/3592421).
- **Use.** Relevant for polycube maps too.

---

### B4. Loop- and sweep-based block decomposition

**B4.1 Livesu, Pietroni, Puppo, Sheffer, Cignoni 2020, LoopyCuts.** **Verdict: reference.**
- **Citation.** *LoopyCuts: Practical Feature-Preserving Block Decomposition for Strongly Hex-Dominant Meshing.*
  ACM TOG 39(4):121, 2020. DOI [10.1145/3386569.3392472](https://doi.org/10.1145/3386569.3392472).
- **Method.** Trace field-aligned surface loops, cut the volume along them into simple blocks, then mesh the blocks.
- **Guarantee.** Strongly hex-dominant output, *not* all-hex.
- **Fit.** Good for CAD with features.
- **Complexity.** High.
- **Code.** `mlivesu/LoopyCuts` has **NO LICENCE**: do not read.

**B4.2 Reberol, Verhetsel, Henrotte, Bommes, Remacle 2023, all-hex boundary layers.** **Verdict: reference (CFD boundary layers later).**
- **Citation.** *Robust Topological Construction of All-hexahedral Boundary Layer Meshes.* ACM TOMS 49(1), 2023.
  DOI [10.1145/3577196](https://doi.org/10.1145/3577196).

---

### B5. Hex-dominant via tet recombination or aligned points

**B5.1 Meshkat and Talmor 2000.** **Verdict: reference.**
- **Citation.** *Generating a mixed mesh of hexahedra, pentahedra and tetrahedra from an underlying tetrahedral mesh.*
  IJNME 49(1–2):17–30, 2000.
  DOI [10.1002/1097-0207(20000910/20)49:1/2<17::AID-NME920>3.0.CO;2-U](https://doi.org/10.1002/1097-0207(20000910/20)49:1/2%3C17::AID-NME920%3E3.0.CO;2-U).
- **Method.** Graph search for tet groups forming hexes or prisms.

**B5.2 Yamakawa and Shimada 2003, rectangular-solid cell packing.** **Verdict: reference.**
- **Citation.** IJNME 57(15):2099–2129, 2003. DOI [10.1002/nme.754](https://doi.org/10.1002/nme.754).
- **Method.** Pack cells to place hex-friendly points, then recombine.

**B5.3 Baudouin, Remacle, Marchandise, Henrotte, Geuzaine 2014, frontal hex-dominant.** **Verdict: reference.**
- **Citation.** Adv. Model. Simul. Eng. Sci. 1:8, 2014. DOI [10.1186/2213-7467-1-8](https://doi.org/10.1186/2213-7467-1-8).
- **Method.** Frontal point insertion along a cross field, Delaunay of the points, then recombination.

**B5.4 Botella, Lévy, Caumon 2016, tet recombination.** **Verdict: reference.**
- **Citation.** *Indirect unstructured hex-dominant mesh generation using tetrahedra recombination.* Comput. Geosci.
  20(3):437–451, 2016. DOI [10.1007/s10596-015-9484-9](https://doi.org/10.1007/s10596-015-9484-9).

**B5.5 Sokolov, Ray, Untereiner, Lévy 2016, Hexahedral-Dominant Meshing.** **Verdict: reference (strong candidate for the hex-dominant product).**
- **Citation.** ACM TOG 35(5):157, 2016. DOI [10.1145/2930662](https://doi.org/10.1145/2930662).
- **Method.**
  1. Frame field.
  2. Global parametrization, tolerating defects.
  3. Point set from the integer iso-values.
  4. Delaunay of the points.
  5. Recombine tets into hexes.
  6. Fill the rest with tets, prisms and pyramids.
- **Result.** Robust hex-dominant output even where the parametrization is defective.
- **Code.** Graphite is GPL-2: FORBIDDEN.

**B5.6 Gao, Jakob, Tarini, Panozzo 2017, field-guided polyhedral agglomeration.** **Verdict: reference.**
- **Citation.** ACM TOG 36(4):114, 2017. DOI [10.1145/3072959.3073676](https://doi.org/10.1145/3072959.3073676).
- **Method.** Instant-Meshes-style position field in 3D, giving a polyhedral, mostly-hex mesh.
- **Code.** `gaoxifeng/robust_hex_dominant_meshing` is **MIT**: READ OK.

**B5.7 Pellerin, Johnen, Verhetsel, Remacle 2018, vertex-based tet-to-hex identification.** **Verdict: build on (the recombination step).**
- **Citation.** *Identifying combinations of tetrahedra into hexahedra: a vertex based strategy.* CAD 105:1–10, 2018.
  DOI [10.1016/j.cad.2018.05.004](https://doi.org/10.1016/j.cad.2018.05.004).
- **Method.** Enumerate **all** potential hexes in a tet mesh exhaustively from vertex neighbourhoods. This is faster
  and complete compared with pattern-matching on tet groups. A selection step (maximum compatible set) then follows.
- **Fit.** Reuses our tet pipeline directly.
- **Code.** Implemented in Gmsh (GPL): FORBIDDEN.

**B5.8 Verhetsel, Pellerin, Remacle 2019, hexahedrizations of small quadrangulations.** **Verdict: skip for now (research).**
- **Citation.** *Finding hexahedrizations for small quadrangulations of the sphere.* ACM TOG 38(4):53, 2019.
  DOI [10.1145/3306346.3323017](https://doi.org/10.1145/3306346.3323017).

---

### B6. Hex quality improvement and untangling

**B6.1 Livesu, Sheffer, Vining, Tarini 2015, edge-cone rectification.** **Verdict: reference.**
- **Citation.** *Practical hex-mesh optimization via edge-cone rectification.* ACM TOG 34(4):141, 2015.
  DOI [10.1145/2766905](https://doi.org/10.1145/2766905).
- **Method.** Local-global optimization that makes each corner's edge cone well-shaped. It fixes inverted elements,
  with surface constraints.

**B6.2 Xu, Gao, Chen 2018, edge-angle optimization.** **Verdict: reference.**
- **Citation.** *Hexahedral mesh quality improvement via edge-angle optimization.* Computers & Graphics 70:17–27, 2018.
  DOI [10.1016/j.cag.2017.07.002](https://doi.org/10.1016/j.cag.2017.07.002).
- **Note.** This is the "Xu et al. 2018" hex improvement and untangling paper.

**B6.3 Marschner, Palmer, Zhang, Solomon 2020, sum-of-squares hex repair.** **Verdict: reference.**
- **Citation.** *Hexahedral Mesh Repair via Sum-of-Squares Relaxation.* CGF 39(5):133–147, 2020.
  DOI [10.1111/cgf.14074](https://doi.org/10.1111/cgf.14074).
- **Method.** Polynomial (SOS) relaxation of the trilinear-hex Jacobian positivity. It certifies validity of the
  *trilinear* element, not just its corners.
- **Assessment.** Heavy machinery, since it needs an SDP solver. Reference for the **validity check**: corner Jacobians
  are insufficient for trilinear hexes.

**B6.4 Garanzha, Kaporin, Kudryavtseva, Protais, Ray, Sokolov 2021, Foldover-free maps in 50 lines of code.** **Verdict: BUILD ON (the untangler for every hex family and for polycube/IGM maps).**
- **Citation.** ACM TOG 40(4):102, 2021.
- **Links.** DOI [10.1145/3450626.3459847](https://doi.org/10.1145/3450626.3459847). arXiv: https://arxiv.org/abs/2102.03069
- **Method.**
  - §2.1: a variational barrier energy, with χ-regularization of det J.
  - §2.2: penalty method. ε_k is driven to zero (§4.2).
  - §2.3: L-BFGS resolution scheme.
- **Result.** Untangles meshes reliably and is simple to implement. The paper reports it on a large benchmark (§3.1).
- **Code.** `ssloy/invertible-maps` is **AGPL-3**: FORBIDDEN. The paper itself is enough, by design.

**B6.5 Visualization and metrics.** **Verdict: reference.**
- **HexaLab** (MIT, https://github.com/mlivesu/HexaLab, https://www.hexalab.net) is READ OK. It is a hex quality
  viewer whose metrics match the literature's.

**B6.6 Element validity and quality metrics (implemented in `fieldcraft-mesh::quality`).** **Verdict: BUILD ON.**
- **Johnen, Remacle, Geuzaine 2013.** *Geometrical validity of curvilinear finite elements.* J. Comput. Phys.
  233:359–372, 2013. The Jacobian
  determinant in Bernstein form is bounded by its coefficients, equals them at the vertices and converges under
  subdivision: a proof of `det J > 0` on the whole element. We use it for every solid kind, including the trilinear
  hex whose corner Jacobians do not certify it (B6.3). (The paper's authors also wrote Gmsh code, which stays
  forbidden; the paper is enough.)
- **Shewchuk 2002.** *What is a good linear finite element? Interpolation, conditioning, anisotropy, and quality
  measures.* Unpublished preprint, https://people.eecs.berkeley.edu/~jrs/papers/elemj.pdf. Tet dihedral angles,
  radius ratio and aspect ratio; with Klingner–Shewchuk's volume–length measure (A3.6).
- **Knupp 2001.** *Algebraic mesh quality metrics.* SIAM J. Sci. Comput. 23(1):193–218. Corner-Jacobian scaled Jacobian and
  shape metrics for hexes, wedges and quads.


### B7. Planar quadrilateral meshing

> Shell (curved-face) quad meshing — advancing front in (u,v) with a cross field, recombination and cleanup — is surveyed separately in [shell-quad-literature.md](shell-quad-literature.md).

 (implemented in `fieldcraft-mesh::quad`)

The face meshers under mapped and swept hex meshing and shell meshing. All quads on prescribed
boundary nodes, so faces sharing an edge share its nodes. DOIs below were checked through
`doi.org/api/handles` and Crossref on 2026-10-13; papers without a DOI are cited by venue only.

**B7.1 Gordon and Hall 1973, transfinite interpolation.** **Verdict: BUILD ON (mapped meshing).**
- **Citation.** W. J. Gordon, C. A. Hall. *Construction of curvilinear co-ordinate systems and applications to mesh
  generation.* IJNME 7(4):461–477, 1973. DOI [10.1002/nme.1620070405](https://doi.org/10.1002/nme.1620070405).
- **Use.** The bilinearly blended Coons map of a four-sided region, dimension-generic (`quad::transfinite_grid`),
  with arc-length blending parameters for unequal spacing on opposite sides.

**B7.2 Free all-quad: the methods compared.**
- **Paving.** T. D. Blacker, M. B. Stephenson. *Paving: a new approach to automated quadrilateral mesh generation.*
  IJNME 32(4):811–847, 1991. DOI [10.1002/nme.1620320410](https://doi.org/10.1002/nme.1620320410). **Verdict:
  reference.** Best boundary rows; front collisions are its failure mode.
- **Q-Morph.** S. J. Owen, M. L. Staten, S. A. Canann, S. Saigal. *Q-Morph: an indirect approach to advancing front
  quad meshing.* IJNME 44(9):1317–1340, 1999. DOI
  [10.1002/(SICI)1097-0207(19990330)44:9<1317::AID-NME532>3.0.CO;2-N](https://doi.org/10.1002/(SICI)1097-0207(19990330)44:9%3C1317::AID-NME532%3E3.0.CO;2-N).
  **Verdict: reference (upgrade path for boundary alignment).** Front states, side-edge recovery by swaps, seam
  closing.
- **Blossom-Quad.** J.-F. Remacle, J. Lambrechts, B. Seny, E. Marchandise, A. Johnen, C. Geuzaine. *Blossom-Quad: a
  non-uniform quadrilateral mesh generator using a minimum-cost perfect-matching algorithm.* IJNME 89(9):1102–1119,
  2012. DOI [10.1002/nme.3279](https://doi.org/10.1002/nme.3279). **Verdict: reference (upgrade path for the
  matching).** Paper only: the implementation lives in Gmsh, which is GPL and forbidden.
- **Recombination and split.** C. K. Lee, S. H. Lo. *A new scheme for the generation of a graded quadrilateral
  mesh.* Computers & Structures 52(5):847–857, 1994. DOI
  [10.1016/0045-7949(94)90070-1](https://doi.org/10.1016/0045-7949(94)90070-1). **Verdict: BUILD ON.** Triangles
  merged into quads, then every element split once, which makes the mesh all-quad whatever the merge found.
  The split topology is that of E. Catmull, J. Clark, CAD 10(6):350–355, 1978, DOI
  [10.1016/0010-4485(78)90110-0](https://doi.org/10.1016/0010-4485(78)90110-0).
- **Why the split method.** All-quad by construction given an even point count per loop; robustness is that of the
  constrained Delaunay triangulation (exists for every valid input; exact predicates); low implementation risk.
  Its cost is irregular nodes and 120° corners in triangle-born quads, reduced by smoothing. The reasoning in full
  is in the `quad` module docs.

**B7.3 Planar constrained Delaunay triangulation and refinement.** **Verdict: BUILD ON.**
- A. Bowyer, *Computing Dirichlet tessellations*, Comput. J. 24(2):162–166, 1981,
  DOI [10.1093/comjnl/24.2.162](https://doi.org/10.1093/comjnl/24.2.162); D. F. Watson, *Computing the n-dimensional
  Delaunay tessellation…*, Comput. J. 24(2):167–172, 1981, DOI
  [10.1093/comjnl/24.2.167](https://doi.org/10.1093/comjnl/24.2.167).
- S. W. Sloan. *A fast algorithm for generating constrained Delaunay triangulations.* Computers & Structures
  47(3):441–450, 1993. DOI [10.1016/0045-7949(93)90239-A](https://doi.org/10.1016/0045-7949(93)90239-A). Edge
  recovery by flipping.
- J. Ruppert. *A Delaunay refinement algorithm for quality 2-dimensional mesh generation.* J. Algorithms
  18(3):548–585, 1995. DOI [10.1006/jagm.1995.1021](https://doi.org/10.1006/jagm.1995.1021); J. R. Shewchuk,
  *Delaunay refinement algorithms for triangular mesh generation*, CGTA 22:21–74, 2002, DOI
  [10.1016/S0925-7721(01)00047-5](https://doi.org/10.1016/S0925-7721(01)00047-5). Used without segment splitting
  (boundary nodes are prescribed).

**B7.3a Free triangle meshes on prescribed boundary nodes** (`free_tri_mesh`, used by
`mesh.planar` for tri3/tri6 2D models): the same constrained Delaunay triangulation of the full
resolution loops, refined to the size by circumcentre insertion that never splits a boundary
segment (J. Ruppert, *A Delaunay refinement algorithm for quality 2-dimensional mesh generation*,
J. Algorithms 18(3):548–585, 1995, without segment splitting, so the boundary nodes stay the
input points), Steiner points kept half the local spacing from the boundary.

**B7.4 Smoothing and untangling of the quads.** **Verdict: BUILD ON.**
- L. A. Freitag. *On combining Laplacian and optimization-based mesh smoothing techniques.* AMD-Vol. 220, ASME 1997;
  S. A. Canann, J. R. Tristano, M. L. Staten, *An approach to combined Laplacian and optimization-based smoothing for
  triangular, quadrilateral, and quad-dominant meshes*, IMR 1998. Laplacian moves kept only when the worst element
  does not get worse.
- L. A. Freitag, P. Plassmann. *Local optimization-based simplicial mesh untangling and improvement.* IJNME
  49(1–2):109–125, 2000. DOI
  [10.1002/1097-0207(20000910/20)49:1/2<109::AID-NME925>3.0.CO;2-U](https://doi.org/10.1002/1097-0207(20000910/20)49:1/2%3C109::AID-NME925%3E3.0.CO;2-U).
  Untangling by maximizing the smallest (linear) corner area.

**B7.5 Quadrilateral quality.** **Verdict: BUILD ON.**
- Knupp 2001 (B6.6) for the scaled Jacobian. J. Robinson, *CRE method of element testing and the Jacobian shape
  parameters*, Eng. Comput. 4(2):113–118, 1987, for aspect ratio and skew from the principal axes. Planar validity
  is exact: a planar bilinear quad's Jacobian is affine in each reference coordinate, so four exact `orient2d` corner
  tests decide it.

### B7a. Node relaxation of a finished mesh (implemented in `fieldcraft-mesh::relax`, `mesh.node.relax`)

- **Method.** Smart Laplacian (S. Canann, J. Tristano, M. Staten, *An approach to combined Laplacian and
  optimization-based smoothing for triangular, quadrilateral, and quad-dominant meshes*, IMR 1998; L. Freitag 1997
  above): a node's weighted step toward the centroid of its edge neighbours is kept only when the worst element
  around it does not get worse and no element inverts, halved otherwise. Then one step up the gradient of the worst
  element's quality (the active max–min element of Freitag and Ollivier-Gooch, IJNME 40, 1997; P. Knupp, IJNME 48,
  2000), kept only when it raises the worst. Quality is the scaled Jacobian of `crate::quality`; validity its
  certified check.
- **Surface nodes in UV.** A node on a B-rep face moves in the face's own parameters (the displacement resolved on
  the surface's coordinate directions, then the surface evaluated there), so it stays exactly on an analytic face;
  free-form faces fall back to projection onto the face's tessellation (SolidCraft's kernel exposes the analytic
  surfaces only). Nodes on a B-rep edge slide along it; vertex nodes are held; partition cut faces are not
  geometry, so their nodes move freely.
- **Relaxation fields.** What a relaxation did is kept as the relaxed region's cells before and after; a later
  mesh takes the displacement at each node's parametric coordinates in the old cell holding it (barycentric or
  trilinear shape functions), slides it onto its face or edge, undoes any move that would invert an element, and
  smooths once (`mesh.node.relax.apply`, or `reapply` after every remesh).
- **Strong mode** (`"mode": "strong"`, `fieldcraft-mesh::relax::morph`). A mesh morph round the bad elements
  (scaled Jacobian below `below`): their movable corner nodes are grouped into clusters (nodes within the reach of
  each other, from a fixed-radius grid query); per cluster, the bad nodes first move to better places by a few
  sweeps of their own optimisation with the rest held (a central-difference gradient step with backtracking on
  Σ SJ⁻¹⁶ + ½ ln²(V/V₀) over their elements, an invalid element costing far more than any valid one), then that
  displacement is carried to every node within the reach by harmonic interpolation: the graph Laplacian of the
  element edges solved (Gauss–Seidel) with the bad nodes' displacements and zero at held nodes and beyond the reach.
  Face nodes go back on their face (UV), edge and vertex nodes are held, nodes outside the reach do not move. A
  cluster's morph is kept, at full or a halved scale, only when no more elements are invalid, the worst element it
  touches is no worse, the smooth objective Σ (1/SJ − 1)² there drops, and no element gets longer than the longest
  was; rounds repeat while any cluster is kept. The result never has a lower smallest SJ than the start.

---

### B8. Mapped and swept hexahedral meshing (implemented in `fieldcraft-mesh::hex`)

`mesh.hex.generate` (and `mesh.hex.sweep`, `mesh.hex.mapped`, `mesh.hex.blocks`). Every method below is implemented from
its paper; no code was read. DOIs marked "checked" were confirmed on 2026-10-13 through a repository
record or the publisher; the others are cited by venue only.

**B8.1 Blacker 1996, the Cooper tool.** **Verdict: BUILD ON (sweep recognition).**
- **Citation.** T. D. Blacker. *The Cooper tool.* Proc. 5th International Meshing Roundtable, 1996, pp. 13–29.
- **Use.** Source, target and linking faces from B-rep topology: rings of logically rectangular
  linking faces followed from the source loops until they close on one target (`hex::sweep`).
  A closed (periodic) face the kernel gives whole, such as a cylinder's side in one piece, takes part
  through its seam: the seam bounds the face on both sides, so its loop is the bottom chain, the seam,
  the top chain and the seam back, a logically rectangular linking face whose two sides are one edge
  (`hex::brep`, the kernel's seams kept as edges from a face to itself); a torus in one face gets the
  level of its revolution from its meridian seam (`hex::layered`).

**B8.2 Lai, Benzley, White 2000, many-to-many sweeping.** **Verdict: BUILD ON (multi-level extrusions).**
- **Citation.** M. Lai, S. Benzley, D. White. *Automated hexahedral mesh generation by generalized multiple
  source to multiple target sweeping.* IJNME 49(1–2):261–275, 2000.
- **Use.** Imprinting every level's outline into one projection plane so each region of the arrangement
  is a one-to-one column (`hex::layered`, `hex::arrangement`); here done in a straight or rotating
  frame, with the regions' fill decided by the generalized winding number (Jacobson et al. 2013).

**B8.3 Knupp 1998, next-generation sweep tool.** **Verdict: BUILD ON (residual correction).**
- **Citation.** P. M. Knupp. *Next-generation sweep tool: a method for generating all-hex meshes on
  two-and-one-half dimensional geometries.* Proc. 7th International Meshing Roundtable, 1998, pp. 505–513.
- **Use.** The boundary residual of each layer's map interpolated into the interior.

**B8.4 Roca and Sarrate 2010, least-squares sweeping.** **Verdict: BUILD ON (inner nodes).**
- **Citation.** X. Roca, J. Sarrate. *An automatic and general least-squares projection procedure for sweep
  meshing.* Engineering with Computers 26(4):391–406, 2010. DOI 10.1007/s00366-009-0172-z (checked, repository
  record). Conference version: Proc. 15th IMR, 2006, pp. 487–506.
- **Use.** Each layer's interior nodes from the source mesh by the least-squares affine map of the source
  ring onto the layer's ring (`hex::sweep::Affine`). Related: Roca, Sarrate, Huerta, *A new least-squares
  approximation of affine mappings for sweep algorithms*, Eng. Comput. 26(3):327–337, 2010.

**B8.5 Interval assignment: Tam and Armstrong 1993; Mitchell 1997/2021.** **Verdict: BUILD ON.**
- **Citations.** T. K. H. Tam, C. G. Armstrong. *Finite element mesh control by integer programming.* IJNME
  36:2581–2605, 1993. S. A. Mitchell. *High fidelity interval assignment.* Proc. 6th IMR, 1997 (IJCGA
  10(4):399–415, 2000). S. A. Mitchell. *Incremental interval assignment by integer linear algebra.*
  International Meshing Roundtable 2021 (Zenodo 5559025); extended in Computer-Aided Design 158:103485, 2023.
- **Use.** Division counts of all curves together: side equalities of mapped faces and linking faces,
  even (or multiple-of-four) loop sums, minimum loop sums, and dominance (`Σ_A ≥ Σ_B + d`: a side of a
  triangle primitive shorter than the other two together). `hex::interval` uses incremental repair from the
  rounded goals with increments only, chosen by squared relative error; it is bounded and reports an
  infeasible system instead of looping. Mitchell's integer-linear-algebra solver (Hermite normal form plus
  nullspace improvement) is the upgrade if the increment-only scheme proves too coarse. A curve's goal is
  its length in units of the size field, raised to a floor from the curve's chordal deviation: the
  geometric size map of a curve (P. J. Frey, P.-L. George, *Mesh Generation: Application to Finite
  Elements*, 2nd ed., Wiley–ISTE, 2008), where a piece of radius `ρ` whose chord stays within the
  tolerance `δ` turns by at most `2 acos(1 − δ / ρ)` (the angle whose sagitta is `δ`), so a curve gets
  `∫ dθ / (2 acos(1 − δ / ρ))` divisions (`size::curve_divisions`), with `δ` the chordal tolerance the
  boundary was tessellated with: the mesh strays from a curve no further than its tessellation, a
  large circle gets many divisions and a small fillet few. No piece turns by more than a quarter turn,
  so a whole circle in one edge keeps at least four (the sides of an O-grid's block). The floor is a
  minimum of the assignment (`hex::curve_goal`).

**B8.6 Whiteley, White, Benzley, Blacker 1996, submapping.** **Verdict: BUILD ON (vertex classification).**
- **Citation.** M. Whiteley, D. White, S. Benzley, T. Blacker. *Two and three-quarter dimensional meshing
  facilitators.* Engineering with Computers 12:144–154, 1996.
- **Use.** Corners by interior angle (≤ 135° end, ≥ 225° reflex); cross-section regions are cut from
  their reflex corners to the boundary (B8.6b), which decomposes rectilinear sections into rectangles (an
  explicit form of submapping's virtual decomposition), then mapped.

**B8.6a Midpoint subdivision: Tam and Armstrong 1991; Li, Armstrong, McKeag 1997.** **Verdict: BUILD ON
(triangle primitive).**
- **Citations.** T. K. H. Tam, C. G. Armstrong. *2D finite element mesh generation by medial axis
  subdivision.* Advances in Engineering Software 13(5–6):313–324, 1991. T. S. Li, C. G. Armstrong,
  R. M. McKeag. *Quad mesh generation for k-sided faces and hex mesh generation for trivalent polyhedra.*
  Finite Elements in Analysis and Design 26(4):279–301, 1997.
- **Use.** A cross-section region with three corners (a fillet's quarter disk, a pocket corner cut off by
  a partition) is split by lines from a centre to one node on each side into three mapped blocks
  (`hex::region`, scheme `Triangle`): side `k` of `n_k` divisions splits into the counts of the two lines
  that end beside it, so each line has `m_k = (n_{k+1} + n_{k+2} − n_k)/2` (an even loop, each side
  shorter than the other two together). The centre is the mean of the split nodes; the interior is then
  smoothed (`quad::smooth`) and every quadrilateral checked; a fold falls back to the free mesher.

**B8.6b Automatic cross-section decomposition (medial-axis subdivision).** **Verdict: BUILD ON (blocks
of 2.5D bodies).**
- **Citations.** Tam and Armstrong 1991 (B8.6a): a region is subdivided at its concave features and the
  branches of its medial axis into simple primitives. Li, Armstrong, McKeag 1997 (B8.6a): the primitives.
  Whiteley et al. 1996 (B8.6): vertex classification.
- **Use.** `hex::decompose` cuts the regions of a multi-level sweep's cross-section that no structured
  scheme fits into blocks, before interval assignment (B8.5), so the cuts are curves of the arrangement
  and their counts are assigned with the rest; nothing reaches the model. Per region: several holes are
  separated by lines across the gaps between neighbouring holes; a piece with one hole gets an O-grid,
  the layout (lines from the outer corners towards the hole or along their bisectors, from the hole's
  corners, or along four directions from the hole) chosen by the score of its pieces; pieces without
  holes are cut by a two-ply search over lines from concave corners and fillets (bisector, edge
  extensions, normals of strongly turning curves) and convex corners' bisectors. In place of the
  medial axis's geometric rules each candidate is scored by meshing its pieces with their schemes (the
  smallest corner scaled Jacobian), preferring fewer free pieces, then less shortfall from 0.9, then
  fewer pieces. A region's blocks are kept when fair (≥ 0.5) or no worse than its free mesh; cuts never
  end inside a seeded curve.

**B8.7 Gordon and Thiel 1982, trivariate transfinite interpolation.** **Verdict: BUILD ON (mapped blocks).**
- **Citation.** W. J. Gordon, L. C. Thiel. *Transfinite mappings and their application to grid generation.*
  Applied Mathematics and Computation 10–11:171–233, 1982.
- **Use.** Interior nodes of a logical cube from its six face grids (`hex::mapped`).

**B8.8 Multi-block meshing of partitioned parts.** **Verdict: BUILD ON (conforming cells).**
- **Citations.** J. F. Thompson, B. K. Soni, N. P. Weatherill (eds.). *Handbook of Grid Generation.* CRC Press,
  1999, ch. 3 (transfinite and algebraic grids) and ch. 13 (multiblock grids). A. Sheffer, M. Bercovier,
  T. Blacker, J. Clements. *Virtual topology operators for meshing.* IJCGA 10(3):309–331, 2000.
- **Use.** `hex::cells` (`mesh.hex.blocks`): cells that meet face to face, each a logical cube after smoothly
  meeting faces are joined into composite faces (virtual topology); vertices, edges and faces shared by
  identity; interval assignment over all cells (B8.5); each edge divided once, each face's transfinite grid
  made once and used by both cells; edges that leave a curved face graded from the spacing along it with a
  growth of 1.2 per element; interiors by B8.7. The merged mesh is checked for free internal faces.

**B8.9 Support.** Generalized winding number for the inside test (A. Jacobson, L. Kavan, O. Sorkine-Hornung,
*Robust inside-outside segmentation using generalized winding numbers*, ACM TOG 32(4):33, 2013; solid angle
of a triangle: Van Oosterom and Strackee 1983); closest point on a triangle (C. Ericson, *Real-Time
Collision Detection*, 2005, §5.1.5).

### B9. Projection blocking (implemented in `fieldcraft-mesh::project`)

A block topology is laid over the part and projected onto its B-rep instead of cutting the
geometry into sweepable cells (the block-structured projection method as publicly described:
blocks in index space, carved, outer faces projected onto surfaces, edges onto curves, vertices
onto points, then smoothed). No code, manual or text of any commercial mesher was read.

- **Transfinite interpolation** places face and interior nodes from projected edges and faces:
  W. J. Gordon, C. A. Hall, *Construction of curvilinear co-ordinate systems and applications to
  mesh generation*, IJNME 7(4):461–477, 1973; W. J. Gordon, L. C. Thiel, *Transfinite mappings and
  their application to grid generation*, Appl. Math. Comput. 10–11:171–233, 1982.
- **Elliptic smoothing**: A. M. Winslow, *Numerical solution of the quasilinear Poisson equation
  in a nonuniform triangle mesh*, J. Comput. Phys. 2:149–172, 1967; J. F. Thompson, Z. U. A. Warsi,
  C. W. Mastin, *Numerical Grid Generation: Foundations and Applications*, North-Holland 1985
  (ch. 6, the three-dimensional equations and their point iteration).
- **Control functions**: P. D. Thomas, J. F. Middlecoff, *Direct control of the grid point
  distribution in meshes generated by elliptic equations*, AIAA J. 18(6):652–656, 1980 (one-
  dimensional form, evaluated on the algebraic grid).
- **Curve and corner projection**: Newton's method on the tangent planes of the two (three)
  surfaces, falling back to alternating projections (von Neumann) where they meet tangentially.
- **Worst-element repair**: compass search on the local smallest scaled Jacobian (optimization-
  based smoothing, L. A. Freitag, C. Ollivier-Gooch, IJNME 40:3979–4002, 1997; metric: P. M. Knupp,
  SIAM J. Sci. Comput. 23(1):193–218, 2001).
- **Ring-and-pillow layout of a mesh plan** (`mesh.plan`, #455): a swept body whose section is at
  an angle to the axes (a crank arm) is blocked as an annulus round each boss's hole and a carved
  grid over the rest, with a layer of blocks inserted between the grid's boundary and the walls
  it faces: pillowing, S. A. Mitchell, T. J. Tautges, *Pillowing doublets: refining a mesh to
  ensure that faces share at most one edge*, Proc. 4th International Meshing Roundtable, 231–240,
  1995, here between a carved grid and the geometry, as the grid-based method projects its
  boundary layer onto the surface (B1.1, Schneiders 1996).

### B10. Extrusion up to a surface (implemented in `fieldcraft-mesh::hex::extrude_to`, `mesh.hex.extrude`)

A quad section is carried along straight per-node paths and its last layer ends on target B-rep
faces, so a region need not be sweepable between two parallel faces.

- **Path ends**: Newton's method on the tangent plane of the target's closest point (the trimmed
  projection of B9), one step exact on a plane; side-face corners by alternating projections
  (J. von Neumann, *Functional Operators, vol. II*, Princeton 1950, ch. 13).
- **Interior nodes**: straight paths plus the boundary's displacement onto the side faces,
  interpolated by inverse-distance weights: the boundary-residual correction of P. M. Knupp,
  *Next-generation sweep tool: a method for generating all-hex meshes on two-and-one-half
  dimensional geometries*, Proc. 7th International Meshing Roundtable, 1998 (as in B8's sweep);
  D. Shepard, *A two-dimensional interpolation function for irregularly-spaced data*, Proc. ACM
  National Conference 1968.
- **Layer grading**: geometric progression of thicknesses (last over first = bias).

### B11. Hex–tet transitions and decomposition proposals (implemented in `fieldcraft-mesh::transition`, `mesh.cells`, `mesh.plan`)

- **Pyramid transition layer**: a pyramid on each quadrilateral of a hex mesh's interface with a
  tetrahedral region, its apex inside the region, so tetrahedra share its triangular sides:
  S. J. Owen, S. A. Canann, S. Saigal, *Pyramid elements for maintaining tetrahedra to hexahedra
  conformability*, AMD Vol. 220, Trends in Unstructured Mesh Generation, ASME 1997, 123–129. The
  tet mesher keeps the pyramids' sides as frozen boundary faces (A3's refinement and improvement
  refuse to change them).
- **Pyramid element**: L. Bedrosian, *Shape functions and integration formulas for
  three-dimensional finite element analysis*, IJNME 35(1):95–108, 1992 (rational functions, a
  collapsed Gauss rule).
- **Decomposition into sweepable cells**: cutting a part so most of it sweeps and a remainder goes
  to a free mesher, as in T. D. Blacker, *The Cooper Tool*, Proc. 5th International Meshing
  Roundtable, 13–29, 1996 (B4, B8); here the remainder is the box round the features the
  provenance plan does not map.

## Part C: Datasets and benchmarks

| Dataset | Content | Licence | Usable for FieldCraft? |
|---|---|---|---|
| **Thingi10K**: Zhou & Jacobson 2016, arXiv https://arxiv.org/abs/1605.04797; site https://ten-thousand-models.appspot.com; tooling https://github.com/Thingi10K/Thingi10K | 10,000 STL meshes. 45% self-intersecting, 50% non-solid. The 4,408 "valid" subset is the standard tet-meshing robustness benchmark (fTetWild, Diazzi 2023). | Tooling Apache-2.0. **Each model has its own licence** (10 licences, including CC-BY-NC and GPL). Filter by the `license` field. | **Conditionally.** Use CC-BY / CC0 / public-domain models for robustness testing. **Provenance risk:** Thingiverse records no authoring tool, and many uploads were modelled in commercial CAD. Under the "no test data produced with commercial CAD products" rule this is an **owner decision**. |
| **HexMe**: DOI 10.5281/zenodo.6642020, https://hexme.algohex.eu | 189 tet meshes from 63 CAD models, with feature tags. | Dataset **CC-BY-4.0**. The source models come from the **ABC dataset (Onshape terms), GrabCAD (GrabCAD terms)**, MAMBO (Apache-2.0), and models made with **Gmsh and Siemens NX** (the paper's own statement). | **Only a filtered subset.** The ABC/Onshape, GrabCAD and Siemens NX models are produced with commercial CAD: exclude them. Gmsh-authored ones are acceptable as data. |
| **MAMBO**: https://gitlab.com/franck.ledoux/mambo | STEP models: basic, simple, medium. Designed for hex block-meshing evaluation. | **Apache-2.0** | **Risky.** The STEP headers I sampled say `'Spatial InterOp 3D'` as originating system and authors `'Cubit 15.2'` / `'MGX'`. These are commercial or closed tools. Treat as commercial-CAD-produced unless the owner rules otherwise. |
| **ABC**: Koch et al., CVPR 2019, DOI [10.1109/CVPR.2019.00983](https://doi.org/10.1109/CVPR.2019.00983), arXiv https://arxiv.org/abs/1812.06216, https://deep-geometry.github.io/abc-dataset/ | 1M CAD models (B-rep and meshes). | "Copyright of the CAD models is owned by their creators", under **Onshape Terms of Use**. The models come from Onshape. | **No.** Produced with a commercial CAD product, and the licence is not open. Skip. |

**Recommendation for test data.** Build FieldCraft's own corpus:
- **Generated cases.** Parametric primitives plus booleans from our own kernel:
  - cubes and cylinders with holes;
  - fillets;
  - thin walls;
  - multi-material interfaces;
  - genus-n plates;
  - a "nasty" series modelled on HexMe's categories by description, not by copying its models.
- **Analytic cases.** Sphere and torus volumes, with exact volumes for checks.
- **Robustness stress.** The CC-licensed Thingi10K subset, if the owner approves.
- **Hex-meshing claims.** Report results on our own corpus. Report the HexMe Gmsh-authored subset separately.

---

## Part D: Recommended pipeline (opinionated)

### Step 1: Robust predicates (weeks 1–3)
1. Shewchuk adaptive `orient2d`, `orient3d`, `incircle` and `insphere` (A0.1, §2–§4).
   - Port the expansion primitives from the paper.
   - Use `predicates.c` (public domain) and `georust/robust` (Apache-2.0) only as differential oracles in tests.
2. Devillers–Teillaud symbolic perturbation for `insphere` (A0.2).
3. Indirect predicates for LPI/TPI implicit points (A0.3, §4–§5). Use a small build-time code generator that emits a
   floating-point filter, an interval filter and an expansion fallback. Model it on PCK (A0.4, BSD-readable). Do not
   use Attene's LGPL code.
4. **Tests.** Exhaustive near-degenerate fuzzing against an exact rational reference (`num-rational` / BigInt in tests only).

**If wrong.** Every later stage inherits silent topological corruption. This is the one stage where "almost right" is fatal.

### Step 2: Delaunay plus boundary recovery (weeks 3–12)
1. Bowyer–Watson with BRIO and Hilbert order, straight walking, and a compact adjacency array (A1.1).
2. **Boundary recovery: implement Diazzi et al. 2023 (A2.3).**
   - Si–Gärtner segment recovery (§3.3) and face recovery (§3.4).
   - Implicit Steiner points (§4.3).
   - The gift-wrapping fallback (§4.4), which is mandatory: it closes the theoretical gap in TetGen.
   - Inside/outside labelling (§4.5).
   - **Why CDT over the envelope approach.** The B-rep tessellator gives a watertight, non-self-intersecting PLC. CAE
     needs exact face conformity: B-rep face IDs on boundary triangles, shared interfaces in assemblies, contact
     surfaces. CDT preserves this by construction. fTetWild (A2.5) does not preserve all input triangles.
3. Keep **fTetWild's design (MPL-2.0, readable) as the documented fallback** for imported dirty geometry. Do not build
   it now. Convex polyhedral meshing (A2.4) is the alternative if exact conformity on self-intersecting input is ever
   required.

**If wrong.** If the tessellator emits self-intersections (common at tolerant B-rep edges), the CDT rejects the input.
- **Mitigation.** Run an exact self-intersection check before CDT.
- **Then route to repair.** Use the arrangement method (A2.6), or the fTetWild path once it exists.

### Step 3: Quality (weeks 12–18)
1. Constrained Delaunay refinement with a radius–edge bound of about 2 and a sizing field (A3.1, A2.2/Si 2010 for small
   input angles).
2. Interleave with ODT vertex smoothing following Tournois 2009 (A3.5). Use the Chen–Holst update formulas (A3.4) and
   keep boundary vertices on their B-rep faces.
3. **Final pass: the Klingner–Shewchuk schedule (A3.6, §3–§4.2).** Smoothing, edge removal, multi-face removal and
   vertex insertion, monotone in worst quality.
4. **Targets.** Report min/max dihedral angles and an FE quality metric. Use Shewchuk's *What is a good linear
   element?*, PDF https://people.eecs.berkeley.edu/~jrs/papers/elemj.pdf, for the metric choice.
5. **CVT/RVD (A4.2–A4.3).** Build only for **surface remeshing** and for Lp-CVT if needed later. Not for tet
   smoothing, where ODT is better.

**If wrong.** ODT alone leaves boundary slivers, which can make FE results wrong. Step 3.3 is the guard. Do not ship without it.

### Step 4: Hex meshing. Build the grid/octree dual family first.

**Order.**
1. **(a) Octree dual all-hex.** Livesu 2021 schemes (B1.5) with Pitzalis 2021 refinement (B1.6). Then the Gao 2019
   stage list (B1.4): feature projection, padding layer, and global min-scaled-Jacobian optimization. Use the Garanzha
   2021 untangler (B6.4) for any inverted elements. Validate with the Marschner 2020 insight (B6.3) that corner
   Jacobians do not certify trilinear hexes; use a subdivided or Bernstein Jacobian check.
2. **(b) Hex-dominant by recombination.** Frame-field-aligned point placement:
   - Ray 2016 field (B3.4);
   - frontal or Sokolov-style points (B5.3/B5.5);
   - our own Delaunay;
   - Pellerin 2018 vertex-based hex identification (B5.7) and selection;
   - pyramids, prisms and tets for the rest.
   This reuses all of steps 1–3.
3. **(c) Polycube for blocky CAD parts.**
   - DualCube labelling (B2.9): guaranteed valid.
   - Volumetric map with an injective solve (B6.4).
   - Robust quantization (B2.7, from the paper only).
   - Pitzalis refinement (B1.6) for grading.
4. **(d) Frame-field/IGM all-hex last, as R&D.** All mature implementations are GPL. HexMe shows 19/189 fully correct
   outputs, and the best 2023 method reaches 58% valid IGMs.

**Why octree first.**
- It is the **only family with a robustness guarantee on arbitrary input** (survey §4.6, quoted in B0.1).
- It is the only family shipped in professional automatic hex meshers.
- Every component has a permissive reference: CinoLib MIT for the dual schemes, Gen-Adapt-Ref MIT for refinement,
  geogram BSD.
- A CAE customer gets *a* valid all-hex mesh for any part on day one, which polycube and frame-field methods cannot promise.

**Risks per family.**

| Family | Main risk | Mitigation |
|---|---|---|
| Octree dual | Poor element quality at the boundary (low scaled Jacobian near concave features); orientation dependence; high element count; edges where more than 3 faces meet; genus not guaranteed (survey §4.6). | Padding layer. Feature-aware projection plus global optimization (B1.4). Align the octree to the part's principal axes. Weakly balanced grids (B1.5 §6) reduce element count. A topology check rejects genus changes. |
| Hex-dominant recombination | Solvers must support pyramids and prisms; hex ratio is often only 60–90%; pyramids hurt accuracy. | Check the solver targets (Abaqus/Nastran element libraries support mixed meshes). Report the hex ratio honestly. |
| Polycube | High distortion on organic or highly non-axis-aligned parts; features at oblique angles are lost; the injective volumetric map can fail. | Use only when the labelling distortion score is acceptable; otherwise fall back to octree. The Garanzha untangler handles the map. |
| Frame-field IGM | Low success rate on CAD (HexMe, Liu & Bommes 2023); mixed-integer solves; all reference code is GPL. | Treat as R&D. Keep it out of the critical path. |

**If the octree-first call is wrong.** Customers wanting high-quality structured hex for thin-walled or extruded parts
would find octree meshes too unstructured. The answer then is sweeping and block decomposition (LoopyCuts-style, B4.1),
not frame fields. That family is not covered in depth here and would need a follow-up survey on sweeping and
multi-sweeping.

---

## Licence summary of associated code (metadata only; no source read)

| Repo / code | Licence | Status |
|---|---|---|
| Shewchuk `predicates.c` | Public domain | READ OK |
| Shewchuk Triangle | Proprietary (non-commercial) | FORBIDDEN |
| `georust/robust` | Apache-2.0 (GitHub metadata) | READ OK |
| `BrunoLevy/geogram` | BSD-3, **except** `third_party/tetgen` (AGPL), `third_party/triangle` and the `delaunay_tetgen/triangle` wrappers | READ OK, with exclusions |
| `wildmeshing/fTetWild` | MPL-2.0 | READ OK |
| `wildmeshing/wildmeshing-toolkit` | MIT | READ OK |
| `Yixin-Hu/TetWild` | GPL-3.0 | FORBIDDEN |
| TetGen | AGPL | FORBIDDEN |
| CGAL, Gmsh (incl. HXT, recombination), Netgen | GPL / LGPL | FORBIDDEN |
| `MarcoAttene/CDT`, `Delaunay3D` | GPL, or LGPL by option | FORBIDDEN |
| `MarcoAttene/Indirect_Predicates` | LGPL-2.1 | FORBIDDEN |
| `MarcoAttene/VolumeMesher` | GPL-3 | FORBIDDEN |
| `libigl/libigl` | MPL-2.0 core; GPL `copyleft/` modules (GitHub reports GPL-3 because both licence files exist) | Core READ OK; `copyleft/` FORBIDDEN |
| `mlivesu/cinolib` | MIT | READ OK |
| `cg3hci/Gen-Adapt-Ref-for-Hexmeshing` | MIT | READ OK |
| `mlivesu/HexaLab` | MIT | READ OK |
| `crawforddoran/quartet` (isosurface stuffing) | ISC | READ OK |
| `alecjacobson/isostuffer` | LGPL-3 | FORBIDDEN |
| `gaoxifeng/robust_hex_dominant_meshing` | MIT | READ OK |
| `gaoxifeng/Feature-Preserving-Octree-Hex-Meshing` | none | NO LICENCE |
| `mlivesu/LoopyCuts` | none | NO LICENCE |
| `LIHPC-Computational-Geometry/evocube` | GPL-3 | FORBIDDEN |
| `LIHPC-Computational-Geometry/dds-hexmeshing` (HexMe-related workflow) | MIT | READ OK |
| `fprotais/robustPolycube`, `ultimaille`, `marchinghex` | AGPL-3 | FORBIDDEN |
| `ssloy/invertible-maps` (foldover-free) | AGPL-3 | FORBIDDEN |
| libHexEx (RWTH) | GPLv3 | FORBIDDEN |
| `cgg-bern/libHexHex` | GPL | FORBIDDEN |
| `HendrikBrueckler/QGP3D`, `MC3D`, `TrulySeamless3D`, `C4HexMeshing` | GPL-3 | FORBIDDEN |
| `mlivesu/CoMISo` (mixed-integer solver fork) | GPL-3 | FORBIDDEN |
| Graphite (`BrunoLevy/GraphiteThree`) | GPL-2 | FORBIDDEN |
| OpenVolumeMesh, CGoGN | Not checked in this pass. Believed LGPL. | Treat as FORBIDDEN until verified |

## Unverified items
- Stellar (Klingner's mesh improvement code): licence not checked.
- `https://people.eecs.berkeley.edu/~jrs/papers/edge.pdf`: listed on Shewchuk's page but not fetched.
- The DualCube code repo: 404 at the name given in search results.
- The HexMe workflow repo URL: not stated in the paper. `dds-hexmeshing` is a related LIHPC repo, not confirmed as the
  HexMe workflow.
