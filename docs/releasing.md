# Releasing FieldCraft

Releases are deliberate. Development lands on `develop`, the trunk; `main` holds released code
only. Fast-forwarding `main` to the chosen, verified `develop` commit builds every package and
creates a **draft** GitHub Release; a maintainer verifies the draft and publishes it. Pushes to
`develop` and tags never start a release, and nothing bumps the version automatically. The workflow is `.github/workflows/release.yml`; the packaging scripts
are in `packaging/` and run the same way on a developer machine.

## What a release contains

| Platform | Files | Script |
|---|---|---|
| macOS 11+, universal (arm64 + x86_64) | `fieldcraft-<v>-macos-universal.dmg` (FieldCraft.app, including `Contents/MacOS/fieldcraft-cli`, on a drag-to-Applications image), `fieldcraft-cli-<v>-macos-universal.zip` | `packaging/macos/package.sh` |
| Linux x86_64, glibc 2.35+ | `fieldcraft-<v>-linux-x86_64.tar.gz` (`bin/`, `share/`), `fieldcraft-<v>-linux-x86_64.AppImage`, `fieldcraft_<v>_amd64.deb` (Debian, Ubuntu 22.04+), `fieldcraft-<v>.x86_64.rpm` (Fedora, RHEL, openSUSE) | `packaging/linux/package.sh` |
| Windows x64 | `fieldcraft-<v>-windows-x64.zip` (`fieldcraft.exe`, `fieldcraft-cli.exe`) | `packaging/windows/package.ps1` |

Plus `SHA256SUMS.txt`, its signature `SHA256SUMS.txt.minisig` (below), and the notes from `.github/release-notes.md` followed by GitHub's
generated list of changes. Every package carries `README.md`, `LICENSE-MIT`, `LICENSE-APACHE`,
`NOTICE` and `THIRD-PARTY-LICENSES.txt` (on macOS in `FieldCraft.app/Contents/Resources/Licenses`).
`THIRD-PARTY-LICENSES.txt` is not committed: each packaging job runs `cargo xtask licences` first,
which fails when a shipped crate names no licence and carries no licence file, and the release
attaches the file on its own too. Which licences a dependency may have is `[licenses]` in
`deny.toml`, which CI checks on every pull request (`cargo deny check licenses`).
The application id is `io.github.jyviko.fieldcraft`: the macOS bundle id, the Linux `.desktop`
entry and icon name, and the window's Wayland app_id / X11 class. The icons come from
`assets/app-icon`: `icon.icns` and `Assets.car` (macOS 26 and later) in the app bundle,
`icon.ico` embedded in `fieldcraft.exe` (`apps/fieldcraft/build.rs`), the iconset PNGs as the
Linux hicolor icons.

The `.deb` and `.rpm` are built by `cargo-deb` and `cargo-generate-rpm` (versions pinned in
`release.yml`) from the same stripped release binaries as the tar.gz, with the same files under
`/usr` (`/usr/bin`, the desktop entry, the hicolor icons, the licences in
`/usr/share/doc/fieldcraft`); their metadata and file lists are in `apps/fieldcraft/Cargo.toml`.
The binaries link only glibc; the windowing, GPU and file-dialog libraries they open at run time
are listed as dependencies by hand (Debian package names, RPM soname capabilities). A
pre-release `0.2.0-rc.1` is versioned `0.2.0~rc.1` inside the packages so it sorts before
`0.2.0`; the file names keep the release version. The Linux job compares both file lists with
the tar.gz, runs lintian, and installs the packages in Ubuntu 22.04 and Fedora containers to run
`fieldcraft-cli --version`.

Native rendering and GPU solvers require a Vulkan loader and a compatible Vulkan driver. D3D/DX12,
Metal and native OpenGL are unsupported; no DXC DLL is needed. Linux packages require the Vulkan
loader, not an EGL alternative. On macOS this requires a separately installed Vulkan portability
implementation such as MoltenVK; the package does not bundle one. A stock macOS installation
without Vulkan cannot open the desktop viewport; CLI CPU solves and CPU snapshots remain usable.
A release must verify its intended macOS Vulkan setup before claiming desktop support there.
Browser rendering remains WebGPU with WebGL2 fallback.

## Signing

**macOS.** When the repository has the secrets below, the app, the dmg and the CLI are signed
with a Developer ID (hardened runtime, `packaging/macos/entitlements.plist`) and notarized, and
the tickets are stapled. Without them the build is **ad-hoc signed**, which runs on Apple silicon
but which Gatekeeper blocks on first launch until the user right-clicks › Open; a run warns
about it, and the release notes say which of the two the draft contains.

| Secret | What it is |
|---|---|
| `APPLE_CERTIFICATE` | The "Developer ID Application" certificate with its private key, exported as .p12 and base64-encoded (`base64 -i cert.p12 \| pbcopy`) |
| `APPLE_CERTIFICATE_PASSWORD` | The .p12's password |
| `KEYCHAIN_PASSWORD` | Optional: the temporary CI keychain's password (random when unset) |
| `APPLE_ID` | The Apple ID used for notarization |
| `APPLE_PASSWORD` | An app-specific password for that Apple ID (appleid.apple.com › Sign-In and Security) |
| `APPLE_TEAM_ID` | The 10-character team id of the Developer ID |

Add them to the `release` environment ([Where the secrets live](#where-the-secrets-live)) with
`gh secret set <NAME> --env release -R tangency-labs/fieldcraft`. With the certificate but without the
three notarization secrets the app is signed and not notarized. A dmg made locally with
`packaging/macos/package.sh` is ad-hoc signed unless `MACOS_SIGN_IDENTITY` names an identity
in the keychain.

**The release key (updates).** The app installs a release only when `SHA256SUMS.txt` is signed
with the FieldCraft release key ([updates.md](updates.md)): an Ed25519 key in minisign's format.
Its public half is built into the app (`fieldcraft_update::PUBLIC_KEY` in
`crates/update/src/lib.rs`):

```
untrusted comment: minisign public key: FEF410EFEF5BAB21
RWQhq1vv7xD0/uAOt8fqERP1hP+1UCM9U9QHTDYTYytWKiuedTkT7Lgr
```

The secret half is the `release` environment secret `FIELDCRAFT_SIGNING_KEY`: the whole secret key file
(two lines), made without a password. With it, the release job signs `SHA256SUMS.txt` with
[rsign2](https://github.com/jedisct1/rsign2) (pinned in `release.yml`), checks the signature
against the public key in the source, and uploads `SHA256SUMS.txt.minisig`. Without it the
draft is built unsigned, the run warns, and the notes say the app cannot install that release
itself (users download it from the page). Anyone can check a release by hand with `minisign -Vm
SHA256SUMS.txt -P RWQhq1vv7xD0/uAOt8fqERP1hP+1UCM9U9QHTDYTYytWKiuedTkT7Lgr`.

Setting it up, once:

```sh
cargo install --locked rsign2 --version 0.6.7
rsign generate -W -p fieldcraft-release.pub -s fieldcraft-release.key   # -W: no password
gh secret set FIELDCRAFT_SIGNING_KEY --env release -R tangency-labs/fieldcraft < fieldcraft-release.key
```

Keep an offline copy of the secret key file (a password manager, an encrypted backup) and
delete it from the working machine. Never commit it.

Rotating it (a key that leaked, or a planned change): generate a new pair, put the new public
key into `PUBLIC_KEY` and the block above, re-sign `crates/update/tests/data/SHA256SUMS.txt`
with it (`rsign sign -W -s <new key> -x crates/update/tests/data/SHA256SUMS.txt.minisig
crates/update/tests/data/SHA256SUMS.txt`, which a test checks), release that version, and only
then replace the secret. Copies that still carry the old public key cannot verify releases
signed with the new one: they still show the update, but Install refuses it ("not signed by the
FieldCraft release key"), and a download from the release page brings the new key. If the old key leaked, also delete the secret at once (releases go out
unsigned until the new one is in place) and say so in the next release notes.

**Linux and Windows** builds are not signed; the checksums are the integrity check. On Windows
SmartScreen warns about an unsigned download (More info › Run anyway) until a code-signing
certificate is added.

## Usage analytics key

The release builds (macOS, Linux and Windows jobs) build in the PostHog project key from the
`release` environment secret `POSTHOG_KEY` as `FIELDCRAFT_POSTHOG_KEY` ([privacy.md](privacy.md)).
Without the secret the builds ask nothing and send nothing, as local, fork and dry-run builds do.
The web build on Pages gets the key only when it runs from `main`: its build job then runs in the
`release` environment (the deploy job already uses `github-pages`, and a job has one environment).

## Where the secrets live

Every release secret (the `APPLE_*` secrets, `KEYCHAIN_PASSWORD`, `FIELDCRAFT_SIGNING_KEY`,
`POSTHOG_KEY`) is a secret of the GitHub Environment `release`, whose deployment branches allow
`main` only. Repository secrets are readable by any workflow a writer can run on any branch;
environment secrets reach only jobs that name the environment and run from an allowed branch, so a
contributor with write access cannot read the signing keys from a branch of their own. The jobs
that read them (`release.yml`'s macos, linux, windows and release jobs, `pages.yml`'s build job)
ask for `release` only on `main`
(`environment: ${{ github.ref == 'refs/heads/main' && 'release' || '' }}`); elsewhere the name is
empty, the job runs in no environment, sees no secrets and builds unsigned. `.github/CODEOWNERS`
keeps the workflows, packaging and updater under the owner's review.

Moving them from repository secrets, once (the owner; each value comes from its own source file,
never the shell history):

```sh
for n in APPLE_CERTIFICATE APPLE_CERTIFICATE_PASSWORD KEYCHAIN_PASSWORD APPLE_ID APPLE_PASSWORD \
         APPLE_TEAM_ID FIELDCRAFT_SIGNING_KEY POSTHOG_KEY; do
  gh secret set "$n" --env release -R tangency-labs/fieldcraft < "$n.txt"
done
gh secret list --env release -R tangency-labs/fieldcraft     # all eight listed
# After a release run on main has signed with them:
for n in APPLE_CERTIFICATE APPLE_CERTIFICATE_PASSWORD KEYCHAIN_PASSWORD APPLE_ID APPLE_PASSWORD \
         APPLE_TEAM_ID FIELDCRAFT_SIGNING_KEY POSTHOG_KEY; do
  gh secret delete "$n" -R tangency-labs/fieldcraft
done
```

## Choose and prepare the version

1. Look at the previous published release (`gh release list -R tangency-labs/fieldcraft`) and the
   changes since it, and choose one version for the batch: a patch for compatible fixes; a minor
   for new capabilities, and while on `0.x` also for breaking changes (say what breaks in the
   notes); `-rc.N` for a release candidate. `1.0` is the owner's decision.
2. Set `[workspace.package] version` in the root `Cargo.toml` and run
   `cargo update --workspace --offline`, which changes only the workspace crates' versions in
   `Cargo.lock`. Land it on `develop` through a pull request like any change, and let CI on
   `develop` pass. Note the commit to release.

## Build the draft

```sh
git fetch origin
git push origin <commit>:main      # fast-forward main to the verified develop commit
```

`main` holds released code only: it moves only by fast-forwarding to a verified `develop`
commit (a maintainer push; the branch is protected), never by merges, rewrites or force-pushes. The run
first verifies the commit with everything CI has (see below); then it builds the three platforms; when
all succeed, the last job checks that every file is there, writes `SHA256SUMS.txt` and the notes,
and creates the draft `FieldCraft v<version>` targeting the pushed commit. Pushing
`main` again with the same version rebuilds and replaces the draft's files (its notes are
kept). The job refuses a version whose tag exists or whose release is already published: bump
the version instead.

## What CI runs when

| When | What runs |
|---|---|
| Push to `develop` | `plan` takes the newest `develop` commit whose CI run passed as the base and runs `cargo xtask affected --base <it>`. `cargo xtask ci --affected <base>`: fmt, assets, layers, unit-factors (always, whole tree), then clippy and the tests of the changed crates and every crate that depends on them, and the wasm check when one of them builds for wasm. `heavy` (`cargo xtask test --heavy`) and `oracle` run only when the change reaches the solvers or the oracle (`solver`), or everything. A newer push cancels the run before it; the next one tests everything since the last green run, so nothing is skipped. |
| Change to `Cargo.toml`, `Cargo.lock`, the toolchain, `xtask/` or `.github/` | Everything: the full gate, heavy and oracle (`everything`). |
| Pull request (ours, a fork's once the maintainer approves the run, or Dependabot's; not a draft) | The same plan against the pull request's merge base: the gate, and `heavy` and `oracle` when the change reaches the solvers or everything, on Blacksmith runners only (never the self-hosted runner). A new push cancels the pull request's previous run; the merge runs its own plan on `develop`. |
| Nightly (03:23 UTC) and `gh workflow run CI --ref <branch>` | Everything: the full `cargo xtask ci`, heavy and oracle on Linux, the oracle with its `nightly` cases (manual runs add the macOS oracle). |
| Weekly (Sunday) | The oracle on macOS. |
| Daily (06:07 UTC) and `gh workflow run solidcraft-compat.yml` | `SolidCraft compatibility`: SolidCraft's latest release saves the design test recipes and `fieldcraft-geometry`'s design tests open them (docs/architecture.md, "SolidCraft pin"). |
| Release (`release.yml`, push to `main` and dry runs) | `verify-ci` (full `cargo xtask ci`), `verify-heavy`, and `verify-oracle` on Linux and macOS, in parallel with the packaging jobs (macOS universal, Linux, Windows), which only upload workflow artifacts. The `release` job, the only one that signs `SHA256SUMS.txt` or creates/updates the draft, needs every verify and packaging job: a failing verify means no draft and no signature. The verify jobs restore develop's CI caches (shared keys); only `main` saves release caches. A run takes about as long as its slowest job, ~30 min (the macOS universal package; it was ~55–60 min when packaging waited for verify). |

`cargo xtask affected --base <rev> [--head <rev>]` prints the plan for any range, e.g. what a
commit needs: `cargo xtask affected --base <sha>^ --head <sha>`. A file outside every crate
(`tests/oracle/…`, `assets/…`) belongs to the crates whose sources name its top-level directory;
`solver` is set when the affected crates include the oracle's (`fieldcraft-testkit`) or one with
`#[ignore]`d heavy tests.

CI has no hardware GPU. `cargo xtask ci` runs the GPU suites (`fieldcraft-lbm-gpu`'s
verification, `fieldcraft-mesh-gpu`'s agreement, the drawing tests) on lavapipe, Mesa's software
Vulkan, with `FIELDCRAFT_GPU=software FIELDCRAFT_REQUIRE_GPU=1`, so a missing adapter fails
instead of skipping (lavapipe's own crash, #638, showed up this way). The heavy job runs with
`FIELDCRAFT_GPU=cpu` on Blacksmith: the GPU throughput benchmark needs a 131 MB buffer and
lavapipe allows 128 MB. Bugs of real Vulkan GPUs and their drivers (such as NVIDIA device loss
#636) show up on none of these, so a change to `crates/gpu` or `crates/lbm-gpu` says
in its pull request which real GPU it ran on (as #641 and #651 did).

## Verify and publish

Publishing is a maintainer's call; an agent preparing a release stops at the draft.

1. Check the run succeeded for the intended commit and the draft targets it
   (`gh release view v<version> -R tangency-labs/fieldcraft`).
2. Download the assets into an empty folder (`gh release download v<version> -R tangency-labs/fieldcraft`),
   run `shasum -a 256 -c SHA256SUMS.txt`, check the signature
   (`minisign -Vm SHA256SUMS.txt -P <the public key above>`, or `rsign verify -P …`), and check
   every platform's files are there.
3. Smoke-test what you can run: open the dmg, drag the app to Applications, launch it with
   `--sample cantilever`; run `fieldcraft-cli --version` from the zip or tar.gz; on Linux run the
   AppImage, or install the `.deb` (`sudo apt install ./fieldcraft_<v>_amd64.deb`) or `.rpm`
   (`sudo dnf install ./fieldcraft-<v>.x86_64.rpm`); on Windows unzip and start `fieldcraft.exe`. Look at the run's logs for the
   platforms you could not try (each job checks its own package).
4. Edit the notes: what is new, fixes, breaking changes, known limitations. Keep the install and
   signing paragraphs; never claim signing the notes don't state.
5. Publish: `gh release edit v<version> -R tangency-labs/fieldcraft --draft=false --latest`
   (a release candidate: `--draft=false --prerelease --latest=false`). Publishing creates the tag
   at the draft's commit; never push release tags by hand.

## Dry runs

`gh workflow run release.yml -R tangency-labs/fieldcraft --ref <branch>` runs the verify jobs and builds every package in parallel, and
uploads them with the notes as the `release-assets` artifact of the run, without creating a
release. `-f draft=true` creates or updates the draft, and is accepted on `main` only.
Locally, `packaging/macos/package.sh --arch aarch64` (or `x86_64`, `universal`) and
`packaging/linux/package.sh` write to `dist/release/` (the `.deb` and `.rpm` need
`cargo install --locked cargo-deb cargo-generate-rpm` at the versions in `release.yml`;
`--formats "tar appimage"` skips them).
