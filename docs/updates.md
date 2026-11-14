# Updates

FieldCraft checks for new releases, verifies them and installs them itself, with your consent at
each step. The code is `crates/update` (`fieldcraft-update`); the window's side is
`crates/ui-egui/src/updates.rs`, the desktop app's `apps/fieldcraft/src/update.rs`.

## Checking

The app asks the GitHub Releases API for the latest published release
(`GET https://api.github.com/repos/tangency-labs/fieldcraft/releases/latest`; drafts and pre-releases
are never offered) and compares its version with the running one (Semantic Versioning order).

- **At start-up**, at most once a day, when **Check for updates automatically** is on (the
  default). It is shown with the first run's question and in **Help › Settings… › Updates**.
  Off means no request at all. A failed automatic check is silent (it goes to the log). Runs
  driven by automation (`--control`) and CI never check automatically.
- **Help › Check for Updates** asks now, and says when FieldCraft is up to date or why the check
  failed.
- **Agents and scripts**: the `app.update.check` command (in the app, and through MCP connected
  to it) and `fieldcraft-cli update --check` reply with `current`, `latest`, `newer`, `url` (the
  release page), `published`, `signed` and `asset` (the file this install would download; `null`
  when it cannot update itself). They only report.

The request carries the FieldCraft version in its `User-Agent`; GitHub sees the request's IP
address, as with any download. Nothing else is sent ([privacy.md](privacy.md)).

A newer release shows a notice in the status bar: **FieldCraft 0.5.9 is available — What's new ·
Install · Skip**. What's new opens the release page; Skip stops offering that version (Settings
can offer it again).

## Verifying

Every release carries `SHA256SUMS.txt` and, signed with the FieldCraft release key,
`SHA256SUMS.txt.minisig` (Ed25519 in minisign's format; the public key is built into the app,
`fieldcraft_update::PUBLIC_KEY`, and printed in [releasing.md](releasing.md)). Installing:

1. downloads the list and its signature, and refuses unless the signature is the release key's;
2. downloads the asset for this install and refuses unless its SHA-256 is the one the list gives;
3. on macOS, refuses an app that fails `codesign --verify --deep --strict`, is not
   FieldCraft (`io.github.tangencylabs.fieldcraft`, or `io.github.jyviko.fieldcraft` up to 0.7.x),
   or is rejected by Gatekeeper (`spctl`) when the installed app is accepted.

A release without a signature cannot be installed from the app: Install opens its page instead.
A development build (run from `target/`) checks but never installs.

## Installing

Install downloads and verifies in the background, unpacks the release beside the installed copy
(`.fieldcraft-update-<version>`, so the swap is a rename on one file system), then asks to restart.
**Restart Now** saves or asks about unsaved changes first, as opening another model does; the
window closes, the autosave is closed, the swap is made and the new version starts. **Later**
installs when you quit. A failed swap puts the previous version back; a staging folder left
behind is removed at the next start.

| Installed as | Asset | Swap |
|---|---|---|
| `FieldCraft.app` (macOS, Apple silicon) | `fieldcraft-<v>-macos-aarch64.dmg` | the bundle is replaced, then opened |
| `fieldcraft-cli` (macOS, Apple silicon) | `fieldcraft-cli-<v>-macos-aarch64.zip` | the binary is replaced |
| A Windows folder | `fieldcraft-<v>-windows-x64.zip` | each file is replaced (a running `.exe` is renamed away, never overwritten); the app restarts |
| An AppImage | `fieldcraft-<v>-linux-x86_64.AppImage` | the file is replaced and started |
| A tar.gz folder | `fieldcraft-<v>-linux-x86_64.tar.gz` | the folder is replaced; when its parent is not writable the app says so and links the download |
| `.deb` / `.rpm` | `fieldcraft_<v>_amd64.deb` / `fieldcraft-<v>.x86_64.rpm` | never replaced by the app: the verified package is downloaded to the app's data folder (`updates/`), and the app offers **Open with Installer** and the `sudo apt install …` / `sudo dnf install …` command |
| Flatpak | none | never replaced by the app: `flatpak update` installs new releases from the [Flatpak repository](releasing.md#flatpak-repository) |

**Help › Reinstall FieldCraft…** downloads and verifies the running version again and swaps it
in the same way, to repair a damaged install.

`fieldcraft-cli update --install` does the same from the command line for the install the CLI
belongs to (the tar.gz folder or Windows folder it is in, or the macOS CLI binary), without
restarting anything; for `.deb`/`.rpm` it prints the package manager command.

## The web app

The browser build has no downloads to install: every half hour the page fetches its own address
past the browser's cache, and when the deployed page loads another build than the running one
(Trunk names each build's files by their content hash), the status bar offers **Reload to
update**. Reloading with unsaved changes asks first, as closing the page does.

## Testing

`FIELDCRAFT_UPDATE_HOST` replaces `https://api.github.com` with another server, and allows plain
`http` for it. It is for tests against a local server only (`fieldcraft_update::mock::Server`,
behind the crate's `mock` feature, plays the Releases API); with it set, the app also checks
automatically in automation runs, so the notice can be looked at:

```sh
FIELDCRAFT_UPDATE_HOST=http://127.0.0.1:8899 fieldcraft-cli update --check
```

The tests cover version order, the daily schedule and no request when checks are off, asset
selection and the swap plan for every install kind, a checksum mismatch and a missing or foreign
signature refused (without downloading the asset), archives that try to write outside their
folder, undoing a failed swap, and the CLI against a local server. They never replace a running
program.
