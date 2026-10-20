//! `oracle [--case NAME]... [--bless NAME|all]... [--quick] [--nightly] [--no-write]`: see
//! [`fieldcraft_testkit::oracle`]; run it as `cargo xtask oracle`.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

use std::process::ExitCode;

use fieldcraft_testkit::oracle::{Options, run};

const USAGE: &str = "usage: oracle [--case NAME]... [--bless NAME|all]... [--quick] [--nightly] [--no-write]";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut options = Options { write: true, ..Options::default() };
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let value = |it: &mut std::slice::Iter<'_, String>| it.next().cloned().ok_or_else(|| format!("{a} needs a value\n{USAGE}"));
        let step = match a.as_str() {
            "--case" => value(&mut it).map(|v| options.cases.push(v)),
            "--bless" => value(&mut it).map(|v| options.bless.push(v)),
            "--quick" => {
                options.quick = true;
                Ok(())
            }
            "--nightly" => {
                options.nightly = true;
                Ok(())
            }
            "--no-write" => {
                options.write = false;
                Ok(())
            }
            other => Err(format!("unknown argument `{other}`\n{USAGE}")),
        };
        if let Err(e) = step {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    }
    match run(&fieldcraft_testkit::corpus::workspace_root(), &options) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
