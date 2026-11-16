//! `fieldcraft-cli --version`, which the release packaging runs to check each build: the version
//! the build was packaged as (`FIELDCRAFT_VERSION`, a nightly's suffix included), else the crate's.

use std::process::Command;

/// The version the build was given, written out here rather than taken from the code under test.
fn packaged() -> &'static str {
    match option_env!("FIELDCRAFT_VERSION") {
        Some(v) if !v.is_empty() => v,
        _ => env!("CARGO_PKG_VERSION"),
    }
}

#[test]
fn version_prints_the_packaged_version() {
    for flag in ["--version", "-V"] {
        let out = Command::new(env!("CARGO_BIN_EXE_fieldcraft-cli")).arg(flag).output().unwrap();
        assert!(out.status.success(), "{flag}: {out:?}");
        assert_eq!(String::from_utf8_lossy(&out.stdout).trim(), format!("fieldcraft-cli {}", packaged()));
    }
}
