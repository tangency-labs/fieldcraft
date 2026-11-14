# Releasing FieldCraft

Releases are deliberate. Development lands on `develop`, the trunk; `main` holds released code
only. Fast-forwarding `main` to the chosen, verified `develop` commit builds every package and
creates a **draft** GitHub Release; a maintainer verifies the draft and publishes it. Pushes to
`develop` and tags never start a release, and nothing bumps the version automatically. The workflow is `.github/workflows/release.yml`; the packaging scripts
are in `packaging/` and run the same way on a developer machine.

## What a release contains

| Platform | Files | Script |
|---|---|---|
| macOS 11+, Apple silicon (aarch64; no Intel build since 0.8.0) | `fieldcraft-<v>-macos-aarch64.dmg` (FieldCraft.app, including `Contents/MacOS/fieldcraft-cli`, on a drag-to-Applications image), `fieldcraft-cli-<v>-macos-aarch64.zip` | `packaging/macos/package.sh` |
| Linux x86_64, glibc 2.35+ | `fieldcraft-<v>-linux-x86_64.tar.gz` (`bin/`, `share/`), `fieldcraft-<v>-linux-x86_64.AppImage`, `fieldcraft_<v>_amd64.deb` (Debian, Ubuntu 22.04+), `fieldcraft-<v>.x86_64.rpm` (Fedora, RHEL, openSUSE) | `packaging/linux/package.sh` |
| Linux x86_64 Flatpak | `fieldcraft-<v>-linux-x86_64.flatpak`: built from source with flatpak-builder on the freedesktop runtime (`org.freedesktop.Platform` 26.08, Rust from the `rust-stable` SDK extension), offline against the crates `cargo vendor` fetched; the job installs it and runs both programs. Not yet attached to the release ([Flatpak repository](#flatpak-repository)) | `packaging/linux/flatpak-bundle.sh`, manifest `packaging/linux/flatpak/io.github.tangencylabs.fieldcraft.yml` |
| Windows x64 | `fieldcraft-<v>-windows-x64.zip` (`fieldcraft.exe`, `fieldcraft-cli.exe`) | `packaging/windows/package.ps1` |

Plus `SHA256SUMS.txt`, its signature `SHA256SUMS.txt.minisig` (below), and the notes from `.github/release-notes.md` followed by GitHub's
generated list of changes. Every package carries `README.md`, `LICENSE-MIT`, `LICENSE-APACHE`,
`NOTICE` and `THIRD-PARTY-LICENSES.txt` (on macOS in `FieldCraft.app/Contents/Resources/Licenses`).
`THIRD-PARTY-LICENSES.txt` is not committed: each packaging job runs `cargo xtask licences` first,
which fails when a shipped crate names no licence and carries no licence file, and the release
attaches the file on its own too. Which licences a dependency may have is `[licenses]` in
`deny.toml`, which CI checks on every pull request (`cargo deny check licenses`).
The application id is `io.github.tangencylabs.fieldcraft`: the macOS bundle id, the Linux
`.desktop` entry and icon name, and the window's Wayland app_id / X11 class. The icons come from
`assets/app-icon`: `icon.icns` and `Assets.car` (macOS 26 and later) in the app bundle,
`icon.ico` embedded in `fieldcraft.exe` (`apps/fieldcraft/build.rs`), the iconset PNGs as the
Linux hicolor icons.

Up to 0.7.x the id was `io.github.jyviko.fieldcraft`. An installed 0.7.x macOS app updates only
to an app signed with that id, so it refuses the first release with the new one: **that
release's notes must tell macOS users to download the dmg once by hand**; from there the app
updates itself again. The updater accepts both ids (`ACCEPTED_BUNDLE_IDS` in
`crates/update/src/install.rs`). On Linux a `.deb` or `.rpm` upgrade replaces the old desktop
entry and icons; pins or favourites made with the old desktop id are lost.

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
`POSTHOG_KEY`, `FLATPAK_GPG_KEY`) is a secret of the GitHub Environment `release`, whose deployment branches allow
`main` only. Repository secrets are readable by any workflow a writer can run on any branch;
environment secrets reach only jobs that name the environment and run from an allowed branch, so a
contributor with write access cannot read the signing keys from a branch of their own. The jobs
that read them (`release.yml`'s macos, linux, flatpak, windows and release jobs, `pages.yml`'s build job)
ask for `release` only on `main`
(`environment: ${{ github.ref == 'refs/heads/main' && 'release' || '' }}`); elsewhere the name is
empty, the job runs in no environment, sees no secrets and builds unsigned. `FLATPAK_GPG_KEY`
reaches only `pages.yml`'s `flatpak` job, which runs on `main` only and builds no code
([Flatpak repository](#flatpak-repository)). `.github/CODEOWNERS`
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
| `gh workflow run corpus-bless.yml --ref <branch> -f parts=<names>` | `Corpus bless`: blesses the named meshing corpus parts on Linux x86_64 (the heavy job's runner and toolchain) and uploads the fixtures that changed as an artifact; it pushes nothing (docs/meshing/corpus.md, "Blessing"). |
| Release (`release.yml`, push to `main` and dry runs) | `verify-ci` (full `cargo xtask ci`), `verify-heavy`, and `verify-oracle` on Linux and macOS, in parallel with the packaging jobs (macOS aarch64, Linux, Flatpak, Windows), which only upload workflow artifacts. The `release` job, the only one that signs `SHA256SUMS.txt` or creates/updates the draft, needs every verify and packaging job but Flatpak, which is not yet part of the release: a failing verify means no draft and no signature. The verify jobs restore develop's CI caches (shared keys); only `main` saves release caches. A run takes about as long as its slowest job, ~30 min with warm caches (the macOS package; it was ~55–60 min when packaging waited for verify). |

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
lavapipe allows 128 MB. Bugs of real GPUs and their drivers (NVIDIA device loss #636, Windows
DX12 #588, Metal) show up on none of these, so a change to `crates/gpu` or `crates/lbm-gpu` says
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
6. Once releases carry the Flatpak bundle, check the Flatpak repository took it: the `Flatpak
   repository` run starts a `Pages` run on `main`, whose `Install from the Flatpak repository` job
   installs the new version from Pages and runs it ([Flatpak repository](#flatpak-repository)).

## Flatpak repository

Pages serves a Flatpak repository with the latest published release (not pre-releases) at
<https://tangency-labs.github.io/fieldcraft/flatpak/>, next to the book and the web app, once
releases carry the bundle (below; until then the site has no repository). Users will install and
update FieldCraft with:

```sh
flatpak install --user https://tangency-labs.github.io/fieldcraft/flatpak/fieldcraft.flatpakref
flatpak update --user
```

(`fieldcraft.flatpakref` adds the remote `fieldcraft`; the freedesktop runtime comes from
Flathub. `fieldcraft.flatpakrepo` adds the remote alone.) How it is made:

- `release.yml` builds `fieldcraft-<v>-linux-x86_64.flatpak` (`packaging/linux/flatpak-bundle.sh`).
  **Not yet in the release:** the job runs on every release run and dry run and uploads the bundle
  as the `flatpak-x86_64` workflow artifact, but the `release` job neither waits for it nor
  attaches it. Once a dry run (`gh workflow run release.yml --ref <branch>` on a branch of this
  repository) has built it, a follow-up change adds `flatpak` to the `release` job's `needs`,
  renames its artifact into `release-*`, adds the file to the set (so it is in `SHA256SUMS.txt`)
  and adds the install lines to the release notes and the README.
- Every Pages deploy (`pages.yml` on `main`) runs `packaging/flatpak-site.sh` in its own `flatpak`
  job, the only one that gets `FLATPAK_GPG_KEY`: a fresh runner with Ubuntu's flatpak, ostree and
  minisign, no Rust toolchain and nothing from the book or web build (which runs crates.io build
  scripts). `deploy` adds its output to the `site` artifact of the `build` job. When the `flatpak`
  job fails, the book and the web app still deploy, without `/flatpak/` (the run warns), so the
  repository is missing until the next good run. The script downloads the latest published
  release's bundle, checks `SHA256SUMS.txt.minisig` with the release key the app is built with (so
  a release without `FIELDCRAFT_SIGNING_KEY` gets no repository) and the bundle against
  `SHA256SUMS.txt`, imports it into a new OSTree repository under `flatpak/repo` and writes the
  `.flatpakrepo`, the `.flatpakref` and an index page. Each deploy holds exactly one release, so a
  docs-only deploy keeps the current one, and the size stays that of one release. Releases from
  before the Flatpak have no bundle, and the site then has no repository.
- Publishing a release starts `flatpak-pages.yml`, which runs `pages.yml` on `main` again when
  the release is the latest. A release event runs on its tag, and only `main` may deploy Pages.
  A release published by a workflow's `GITHUB_TOKEN` starts no workflows; a maintainer publishing
  it (on the page, or `gh release edit` as above) does.
- After a deploy that publishing (or a manual run) started, `pages.yml`'s `flatpak-check` job
  waits until Pages serves the new release (its cache keeps files up to 10 minutes), installs
  FieldCraft from the `.flatpakref` on a clean runner and runs both programs.

To put a chosen published release in (going back after a bad release, say):
`gh workflow run pages.yml -R tangency-labs/fieldcraft --ref main -f release_tag=v<version>`.

**Signing (optional).** With the `release` environment secret `FLATPAK_GPG_KEY` (an
ASCII-armored GPG secret key without a passphrase), the commits and the repository summary are
signed and the `.flatpakrepo` and `.flatpakref` carry the public key, so clients verify what they
download. Without it the repository is unsigned (the run warns): clients add it with GPG
verification off and rely on HTTPS. Setting it up, once (the owner):

```sh
export GNUPGHOME="$(mktemp -d)"
gpg --batch --passphrase '' --quick-gen-key "FieldCraft Flatpak repository" rsa4096 sign never
gpg --batch --armor --export-secret-keys | gh secret set FLATPAK_GPG_KEY --env release -R tangency-labs/fieldcraft
```

Keep an offline copy of the key and delete `$GNUPGHOME`. Remotes added before the key existed
stay unsigned until users add the repository again (`flatpak remote-delete --user fieldcraft`,
then install from the `.flatpakref`); the same goes for a new key after a rotation.

## Dry runs

`gh workflow run release.yml -R tangency-labs/fieldcraft --ref <branch>` runs the verify jobs and builds every package in parallel, and
uploads them with the notes as the `release-assets` artifact of the run, without creating a
release. `-f draft=true` creates or updates the draft, and is accepted on `main` only.
Locally, `packaging/macos/package.sh --arch aarch64` (or `x86_64`, `universal`) and
`packaging/linux/package.sh` write to `dist/release/` (the `.deb` and `.rpm` need
`cargo install --locked cargo-deb cargo-generate-rpm` at the versions in `release.yml`;
`--formats "tar appimage"` skips them).
