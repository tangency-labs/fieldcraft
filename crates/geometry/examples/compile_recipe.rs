//! Compile a SolidCraft recipe script (`{"commands": [...]}`) into a `.solidcraft` design that
//! `geometry.open` (and so `fieldcraft-cli run`/`snapshot`) can open. Used by the cooling-examples
//! runner (`tools/cooling/make.py`) to build the committed example geometry from its recipe.
//!
//! `cargo run --example compile_recipe -p fieldcraft-geometry --features scripts -- recipe.json out.solidcraft`
fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (recipe, out) = match args.as_slice() {
        [recipe, out] => (recipe, out),
        _ => return Err("usage: compile_recipe <recipe.json> <out.solidcraft>".into()),
    };
    let script = std::fs::read_to_string(recipe).map_err(|e| format!("read {recipe}: {e}"))?;
    let bytes = fieldcraft_geometry::script_design(&script).map_err(|e| format!("compile: {e}"))?;
    std::fs::write(out, &bytes).map_err(|e| format!("write {out}: {e}"))?;
    eprintln!("compile_recipe: wrote {out} ({} bytes)", bytes.len());
    Ok(())
}
