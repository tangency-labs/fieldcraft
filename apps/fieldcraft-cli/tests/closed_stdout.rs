//! The CLI ends quietly when whoever reads its output stops early (`fieldcraft-cli commands | head`):
//! a closed stdout is not a crash.

use std::process::{Command, Stdio};

/// A folder removed when the test ends, passed or not.
struct TempDir(std::path::PathBuf);

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn output_ends_without_a_panic_when_stdout_is_closed() {
    // Every subcommand that writes to stdout without a network or an input file; `crashes` reads
    // an empty folder of the test's own, removed when the test ends.
    let crash_dir = TempDir(std::env::temp_dir().join(format!("fieldcraft-cli-closed-stdout-{}", std::process::id())));
    let all: [&[&str]; 9] = [
        &["commands"],
        &["commands", "--json"],
        &["features"],
        &["--help"],
        &["--version"],
        &["crashes", "list"],
        &["crashes", "list", "--json"],
        &["crashes", "dir"],
        &["crashes", "clear"],
    ];
    for args in all {
        // The reading end is gone before the first line is written, as when `head` has had enough.
        let (reader, writer) = std::io::pipe().unwrap();
        drop(reader);
        let mut cli = Command::new(env!("CARGO_BIN_EXE_fieldcraft-cli"));
        cli.args(args).env("FIELDCRAFT_CRASH_DIR", &crash_dir.0);
        let out = cli.stdout(writer).stderr(Stdio::piped()).output().unwrap();
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!stderr.contains("panicked"), "{args:?} panicked: {stderr}");
        assert!(out.status.success(), "{args:?}: {:?} {stderr}", out.status);
    }
}

/// `mcp` answers requests from stdin: a reader that left before the reply (`… | fieldcraft-cli mcp
/// | head -n 1` with more requests behind) ends it quietly too.
#[test]
fn mcp_ends_without_an_error_when_stdout_is_closed() {
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    let mut cli = Command::new(env!("CARGO_BIN_EXE_fieldcraft-cli"))
        .args(["mcp", "--headless"])
        .stdin(Stdio::piped())
        .stdout(writer)
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let requests = b"{\"jsonrpc\":\"2.0\",\"id\":1,\"method\":\"ping\"}\n{\"jsonrpc\":\"2.0\",\"id\":2,\"method\":\"ping\"}\n";
    // The CLI may be gone before it has read everything; that is not what is tested.
    let _ = std::io::Write::write_all(cli.stdin.as_mut().unwrap(), requests);
    drop(cli.stdin.take());
    let out = cli.wait_with_output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!stderr.contains("panicked"), "mcp panicked: {stderr}");
    assert!(!stderr.contains("error"), "mcp: {stderr}");
    assert!(out.status.success(), "mcp: {:?} {stderr}", out.status);
}

/// Only a closed reader ends the output quietly: another write error (here a full disk) is reported
/// and fails the command, as the output is incomplete.
#[cfg(target_os = "linux")]
#[test]
fn a_full_disk_on_stdout_is_an_error() {
    for args in [&["--version"][..], &["commands", "--json"]] {
        let full = std::fs::OpenOptions::new().write(true).open("/dev/full").unwrap();
        let out = Command::new(env!("CARGO_BIN_EXE_fieldcraft-cli")).args(args).stdout(full).stderr(Stdio::piped()).output().unwrap();
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!stderr.contains("panicked"), "{args:?} panicked: {stderr}");
        assert!(!out.status.success(), "{args:?} succeeded on a full disk: {stderr}");
        assert!(stderr.contains("writing to stdout"), "{args:?}: {stderr}");
    }
}

/// The same for `mcp`'s replies: a full disk fails it as a stdout write, not as an MCP error.
#[cfg(target_os = "linux")]
#[test]
fn a_full_disk_under_mcp_is_a_stdout_error() {
    let full = std::fs::OpenOptions::new().write(true).open("/dev/full").unwrap();
    let mut cli = Command::new(env!("CARGO_BIN_EXE_fieldcraft-cli"))
        .args(["mcp", "--headless"])
        .stdin(Stdio::piped())
        .stdout(full)
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let _ = std::io::Write::write_all(cli.stdin.as_mut().unwrap(), b"{\"jsonrpc\":\"2.0\",\"id\":1,\"method\":\"ping\"}\n");
    drop(cli.stdin.take());
    let out = cli.wait_with_output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!stderr.contains("panicked"), "mcp panicked: {stderr}");
    assert!(!out.status.success(), "mcp succeeded on a full disk: {stderr}");
    assert!(stderr.contains("writing to stdout"), "mcp: {stderr}");
}
