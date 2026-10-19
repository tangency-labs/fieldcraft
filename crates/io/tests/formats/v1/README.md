# Format 1 fixtures

Real `.fieldcraft` files of format `fieldcraft/1`, one folder per program that wrote them. Format
1 was written by FieldCraft 0.4.0 to 0.7.0, and what it holds grew release by release (inclined
symmetry planes, thermal fluid domains, ties, contact, datums, view state, nonlinear and
buckling steps, …), so every release's files are here:

| Folder | Written by | Commit |
|---|---|---|
| `0.7.0/` | FieldCraft 0.7.0 (release tag `v0.7.0`) | `e9ccf106` |
| `0.6.0/` | FieldCraft 0.6.0 (`v0.6.0`) | `c012a9e4` |
| `0.5.0/` | FieldCraft 0.5.0 (`v0.5.0`) | `e3bdf963` |
| `0.4.0/` | FieldCraft 0.4.0 (`v0.4.0`) | `0c9952ba` |
| `develop-2eb3515c/` | an unreleased `develop` build, 35 commits after `v0.7.0` (its generator still says 0.7.0) | `2eb3515c` |

The `develop-2eb3515c` files hold what no release wrote (a fluid step's `steady` stop rule), so
no user has files like them; they are kept as a check on what format 1 could hold at its end.

`crates/io/tests/formats.rs` opens every one with the current reader and checks that it gives
what its writer's own reader gave (after the migration steps' expected changes), that the step
1 → 2 changes only the tag, and that saving the migrated model and reading it again round-trips
bit for bit. The engine's tests open them through `file.open` and check that opening writes
nothing and the first save keeps the original.

**Never edit, regenerate or replace these files.** They stand for the files users have on disk;
a reader change that breaks one of them breaks those files too. A new case is a new file.

| Fixture | What it holds |
|---|---|
| `hex-static` | a hex8 block (`mesh.block`) with named faces, mm-t-s-N-MPa, a static step with a clamp and a pressure |
| `tet-frequency` | a tet10 block, m-kg-s-N-Pa, a frequency step |
| `tet-design` (0.7.0 only) | a tet4 mesh of the design `hollow-box.solidcraft` (a copy of `crates/geometry/tests/designs/hollow-box.solidcraft`, referenced by a relative path), with geometry tags on the B-rep faces |
| `fluid-box` | a fluid domain with every box face condition and a fluid step (thermal settings from 0.7.0; the `steady` rule in the develop build only); no mesh |
| `unspecified-units` | a wedge6 block in the `unspecified` unit system, a transient heat transfer step |
| `every-kind` | two hex8 blocks with every boundary condition and load kind, reference points, both coupling kinds, point masses, a spring and a dashpot, plastic, damaged and orthotropic materials and every step procedure the release had; from 0.6.0 also a tie, contact pairs and general contact, a rigid wall, datum planes, axes and points and hidden view items; from 0.5.0 a hyperelastic material (neo-Hookean in 0.5.0, Mooney–Rivlin later) |
| `empty` | an empty model |
| `design-metres` (0.7.0 only) | a tet4 mesh of `hollow-box.solidcraft` with an edge seed and two size controls, in m-kg-s-N-Pa: 0.7.0 solved the design's millimetre coordinates as metres; format 3 reads the mesh in millimetres and says so |
| `design-inches` (0.7.0 only) | the same in in-lbf-s-psi |
| `mixed-metres` (0.7.0 only) | a hex8 block typed in metres beside a tet4 mesh of the design, m-kg-s-N-Pa: each body migrates by its origin |
| `thermal-inches` (0.7.0 only) | a hex8 block in in-lbf-s-psi with a thermal material that leaves out its reference temperature, a Johnson–Cook material without a reference strain rate, a plane-strain section without a thickness, a transient heat step without an initial temperature, a convection load and a fluid domain whose heat settings leave out the initial temperature (format 2's defaults, in this system) |
| `thermal-millisecond-kn`, `thermal-millisecond-gram` (0.7.0 only) | the same in mm-kg-ms-kN-GPa and mm-g-ms-N-MPa, beside a tet4 mesh of the design with edge seeds (count with bias, size), size controls (faces, box), a size growth rate, a hex-meshing plan and a relaxation field (`mesh.node.relax`) |

| `fluid-metres` (0.7.0 only) | a heated duct in m-kg-s-N-Pa (an inlet, a moving wall and held temperatures, in kelvin) with a fluid step: its lattice counts from 0 K, as 0.7.0 ran it |

The 0.7.0 release had no porous zones or fans, and stitching needs an imported 3MF surface, so
no fixture holds those.

## How they were made

Each `<name>.recipe.json` is the script that wrote `<name>.fieldcraft`, run in the folder with
that release's `fieldcraft-cli run <name>.recipe.json` (its last step saves the file). The
recipes are the same model trimmed to what each release had: a command or field a release did
not know was left out, nothing was added.

Each `<name>.read.json` is what that release's reader (`fieldcraft_io::read` at the same commit)
gave for the file, as one line of JSON: the whole model (`serde_json::to_value(&model)`, mesh
included), the geometry, the journal and the view (`null` for 0.4.0 and 0.5.0, whose reader had
no view state). It was written by this test, run once at that commit with `FX` set to the folder:

```rust
//! Records what this release's reader gives for every `.fieldcraft` file in `$FX`.
#[test]
fn dump() {
    let dir = std::path::PathBuf::from(std::env::var("FX").unwrap());
    for e in std::fs::read_dir(&dir).unwrap() {
        let p = e.unwrap().path();
        if p.extension().and_then(|x| x.to_str()) != Some("fieldcraft") { continue; }
        let f = fieldcraft_io::read(&std::fs::read(&p).unwrap()).unwrap();
        let v = serde_json::json!({
            "model": serde_json::to_value(&f.model).unwrap(),
            "geometry": serde_json::to_value(&f.geometry).unwrap(),
            "journal": f.journal.iter().map(|(c, p)| serde_json::json!({"command": c, "params": p})).collect::<Vec<_>>(),
            "view": serde_json::to_value(&f.view).unwrap(),
        });
        let mut text = serde_json::to_vec(&v).unwrap();
        text.push(b'\n');
        std::fs::write(p.with_extension("read.json"), text).unwrap();
    }
}
```
