FieldCraft is open-source CAE pre-processing in Rust: it turns geometry into analysis-ready
finite element models, writes Abaqus, Nastran and LS-DYNA decks, and solves linear static and
modal problems with its own solver. AI agents can drive it over MCP (`fieldcraft-cli mcp`).
[README](https://github.com/tangency-labs/fieldcraft#readme) ·
[roadmap](https://github.com/tangency-labs/fieldcraft/blob/main/ROADMAP.md) ·
[issues](https://github.com/tangency-labs/fieldcraft/issues)

**Which file?**
- **macOS 11+** (Apple silicon): open `fieldcraft-{version}-macos-aarch64.dmg` and
  drag `FieldCraft.app` onto Applications. The command-line tool is included at
  `/Applications/FieldCraft.app/Contents/MacOS/fieldcraft-cli` and also available separately in
  `fieldcraft-cli-{version}-macos-aarch64.zip`.
- **Linux** (x86_64, glibc 2.35+: Ubuntu 22.04, Debian 12, Fedora 36 or newer):
  `fieldcraft-{version}-linux-x86_64.AppImage` (`chmod +x`, then run it), or unpack
  `fieldcraft-{version}-linux-x86_64.tar.gz` and run `bin/fieldcraft` (`bin/fieldcraft-cli` is the
  command-line tool).
- **Windows 10+** (x64): unzip `fieldcraft-{version}-windows-x64.zip` and run `fieldcraft.exe`
  (`fieldcraft-cli.exe` is the command-line tool).
- `SHA256SUMS.txt`: `sha256sum -c SHA256SUMS.txt --ignore-missing` (macOS:
  `shasum -a 256 -c SHA256SUMS.txt --ignore-missing`).

To start with a sample model, run the app with `--sample cantilever` (macOS:
`/Applications/FieldCraft.app/Contents/MacOS/FieldCraft --sample cantilever`).

**Updates.** {update_signing}

**Signing.** {macos_signing} The Linux and Windows builds are not signed; on Windows, SmartScreen
may say "Windows protected your PC": choose **More info**, then **Run anyway**.

Every package carries `LICENSE-MIT`, `LICENSE-APACHE`, `NOTICE` and `THIRD-PARTY-LICENSES.txt`.
