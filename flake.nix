{
  description = "FieldCraft — native app development shell and package for Linux and macOS";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  # Official Rust builds (bundled LLVM) at the repo's pinned version.
  inputs.rust-overlay = {
    url = "github:oxalica/rust-overlay";
    inputs.nixpkgs.follows = "nixpkgs";
  };

  # Adapted from SolidCraft's flake.nix (MIT OR Apache-2.0; see ATTRIBUTION.md).
  outputs =
    {
      self,
      nixpkgs,
      rust-overlay,
    }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "aarch64-darwin"
      ];
      forAllSystems =
        f:
        nixpkgs.lib.genAttrs systems (
          system:
          f (
            import nixpkgs {
              inherit system;
              overlays = [ rust-overlay.overlays.default ];
            }
          )
        );

      # The Rust of rust-toolchain.toml (channel, components, targets), as rustup and CI use it;
      # plus rust-src and rust-analyzer for editors. Not nixpkgs' rustc: its 1.99.0 (system
      # LLVM 21.1.8) segfaults in the SLP vectorizer on fieldcraft-solver.
      rustToolchain = (nixpkgs.lib.importTOML ./rust-toolchain.toml).toolchain;
      rustVersion = rustToolchain.channel;
      toolchain =
        pkgs:
        pkgs.rust-bin.stable.${rustVersion}.default.override {
          extensions = rustToolchain.components or [ ] ++ [
            "rust-src"
            "rust-analyzer"
          ];
          targets = rustToolchain.targets or [ ];
        };
      rustPlatform =
        pkgs:
        pkgs.makeRustPlatform {
          cargo = toolchain pkgs;
          rustc = toolchain pkgs;
        };

      cargoToml = nixpkgs.lib.importTOML ./Cargo.toml;
      appId = "io.github.tangencylabs.fieldcraft";
      # SolidCraft checkout at the rev Cargo.lock pins (the same fetch importCargoLock makes).
      solidcraftSrc =
        let
          pkg = nixpkgs.lib.findFirst (
            p: p.name == "solidcraft-kernel"
          ) null (nixpkgs.lib.importTOML ./Cargo.lock).package;
          rev = nixpkgs.lib.last (nixpkgs.lib.splitString "#" pkg.source);
        in
        builtins.fetchGit {
          url = "https://github.com/tangency-labs/solidcraft";
          inherit rev;
          allRefs = true;
        };
      cargoPackages = [
        "-p"
        "fieldcraft"
        "-p"
        "fieldcraft-cli"
      ];

      # Linux only: libraries dlopened at runtime, absent from LD_LIBRARY_PATH on NixOS.
      # - wgpu (app and CLI snapshots): EGL loader (libglvnd), Vulkan loader
      # - winit (app): Wayland client / X11 libraries, libxkbcommon
      # - rfd (app): xdg-desktop-portal through libdbus-1
      # GPU drivers (and their EGL vendor / Vulkan ICD files) come from the system:
      # /run/opengl-driver (NixOS hardware.graphics.enable), or nixGL elsewhere.
      # macOS: AppKit + Metal, nothing needed.
      gpuLibs =
        pkgs:
        pkgs.lib.optionals pkgs.stdenv.hostPlatform.isLinux (
          with pkgs;
          [
            libglvnd
            vulkan-loader
          ]
        );
      runtimeLibs =
        pkgs:
        gpuLibs pkgs
        ++ pkgs.lib.optionals pkgs.stdenv.hostPlatform.isLinux (
          with pkgs;
          [
            wayland
            libxkbcommon
            dbus
            libx11
            libxcursor
            libxi
            libxrandr
            libxrender
            libxext
            libxcb
          ]
        );
      driverLibDir = pkgs: if pkgs.stdenv.hostPlatform.isLinux then "/run/opengl-driver/lib" else "";
      driverShareDir = pkgs: if pkgs.stdenv.hostPlatform.isLinux then "/run/opengl-driver/share" else "";
      mkLibraryPath =
        pkgs: libs:
        if pkgs.stdenv.hostPlatform.isLinux then
          "${pkgs.lib.makeLibraryPath libs}:${driverLibDir pkgs}"
        else
          "";
      libraryPath = pkgs: mkLibraryPath pkgs (runtimeLibs pkgs);

      shellHook = pkgs: ''
        ${pkgs.lib.optionalString pkgs.stdenv.hostPlatform.isLinux ''
          if [ ! -d ${driverLibDir pkgs} ]; then
            echo "fieldcraft: ${driverLibDir pkgs} is missing — on NixOS enable hardware.graphics; elsewhere use a matching nixGL wrapper" >&2
          fi
          export XDG_DATA_DIRS="${driverShareDir pkgs}:''${XDG_DATA_DIRS:-/usr/local/share:/usr/share}"
        ''}
        command -v cc >/dev/null || echo "fieldcraft: no C linker on PATH, cargo cannot link — use nix develop .#rust" >&2
        echo "fieldcraft: cargo run --release -p fieldcraft -- --sample cantilever"
      '';
    in
    {
      devShells = forAllSystems (pkgs: {
        # Default: no toolchain of its own (Rust or C); uses the one on PATH (rustup, system),
        # so its build cache and linker stay consistent. Only adds the runtime library path.
        default = pkgs.mkShellNoCC {
          name = "fieldcraft";
          LD_LIBRARY_PATH = libraryPath pkgs;
          shellHook = shellHook pkgs;
        };

        # No Rust on the machine: `nix develop .#rust`. Separate target dir under target/
        # (gitignored), so switching toolchains does not thrash the cache.
        rust = pkgs.mkShell {
          name = "fieldcraft-rust";
          packages = [ (toolchain pkgs) ];
          LD_LIBRARY_PATH = libraryPath pkgs;
          CARGO_TARGET_DIR = "target/nix";
          shellHook = ''
            ${shellHook pkgs}
            echo "fieldcraft: $(rustc --version), building into target/nix"
          '';
        };
      });

      packages = forAllSystems (pkgs: rec {
        default = fieldcraft;
        # The CLI alone: the `cli` output of the one build (the workspace compiles once).
        fieldcraft-cli = pkgs.symlinkJoin {
          name = "fieldcraft-cli-${fieldcraft.version}";
          paths = [ fieldcraft.cli ];
          meta = fieldcraft.meta // {
            description = "FieldCraft command line: scripts, snapshots, MCP server";
            mainProgram = "fieldcraft-cli";
          };
        };
        # The desktop app (out) and the CLI (cli), one build.
        fieldcraft =
          assert pkgs.lib.assertMsg
            (pkgs.lib.versionAtLeast rustVersion cargoToml.workspace.package.rust-version)
            "FieldCraft requires the Rust version in Cargo.toml (${cargoToml.workspace.package.rust-version}); update rustVersion.";
          (rustPlatform pkgs).buildRustPackage {
            pname = "fieldcraft";
            version = cargoToml.workspace.package.version;
            src = self;
            outputs = [
              "out"
              "cli"
            ];

            cargoLock = {
              lockFile = ./Cargo.lock;
              # SolidCraft, truck, egui-wgpu, winit: git dependencies pinned by rev in Cargo.toml.
              allowBuiltinFetchGit = true;
            };

            # SolidCraft crates embed files from their repo root (`../../assets`, LICENSE-*,
            # NOTICE, THIRD-PARTY-LICENSES.txt); a vendored crate dir has none of them. Vendored
            # crates sit in $NIX_BUILD_TOP/cargo-vendor-dir/<crate>, so `../..` is $NIX_BUILD_TOP.
            preBuild = ''
              for f in assets LICENSE-MIT LICENSE-APACHE NOTICE THIRD-PARTY-LICENSES.txt; do
                ln -s ${solidcraftSrc}/$f "$NIX_BUILD_TOP/$f"
              done
            '';

            cargoBuildFlags = cargoPackages;
            cargoInstallFlags = cargoPackages;
            # Test gate is `cargo xtask ci` on CI (wasm target, GPU); the sandbox has neither.
            doCheck = false;

            # Crash reports record the commit (crates/engine/src/crash.rs).
            FIELDCRAFT_BUILD_SHA = self.rev or self.dirtyRev or "unknown";

            nativeBuildInputs = [ pkgs.makeWrapper ];
            buildInputs = runtimeLibs pkgs;

            # cargo installs both binaries into $out; the CLI moves to its own output.
            # Linux: same share/ tree as packaging/env.sh install_share.
            postInstall = ''
              install -Dm755 $out/bin/fieldcraft-cli $cli/bin/fieldcraft-cli
              rm $out/bin/fieldcraft-cli
            ''
            + pkgs.lib.optionalString pkgs.stdenv.hostPlatform.isLinux ''
              install -Dm644 packaging/linux/${appId}.desktop $out/share/applications/${appId}.desktop
              install -Dm644 packaging/linux/${appId}.metainfo.xml $out/share/metainfo/${appId}.metainfo.xml
              for icon in 16:icon_16x16 32:icon_32x32 64:icon_32x32@2x 128:icon_128x128 \
                256:icon_256x256 512:icon_512x512 1024:icon_512x512@2x; do
                size=''${icon%%:*}
                install -Dm644 assets/app-icon/icon.iconset/''${icon#*:}.png \
                  $out/share/icons/hicolor/''${size}x''${size}/apps/${appId}.png
              done
            '';

            # Linux: the app dlopens its window/GPU/portal libraries, so the binary gets the
            # shell's LD_LIBRARY_PATH. zenity: rfd's file dialog fallback without a portal,
            # after the system's own on PATH. The CLI only renders offscreen: GPU libraries only.
            postFixup = pkgs.lib.optionalString pkgs.stdenv.hostPlatform.isLinux ''
              wrapProgram $out/bin/fieldcraft \
                --prefix LD_LIBRARY_PATH : "${libraryPath pkgs}" \
                --suffix PATH : "${pkgs.lib.makeBinPath [ pkgs.zenity ]}" \
                --set-default XDG_DATA_DIRS /usr/local/share:/usr/share \
                --prefix XDG_DATA_DIRS : "$out/share:${driverShareDir pkgs}"
              wrapProgram $cli/bin/fieldcraft-cli \
                --prefix LD_LIBRARY_PATH : "${mkLibraryPath pkgs (gpuLibs pkgs)}" \
                --set-default XDG_DATA_DIRS /usr/local/share:/usr/share \
                --prefix XDG_DATA_DIRS : "${driverShareDir pkgs}"
            '';

            meta = {
              description = "Open-source CAE modeler: meshing, model setup, decks and solvers";
              homepage = cargoToml.workspace.package.repository;
              license = with pkgs.lib.licenses; [
                mit
                asl20
              ];
              mainProgram = "fieldcraft";
            };
          };
      });

      apps = forAllSystems (
        pkgs:
        let
          packages = self.packages.${pkgs.stdenv.hostPlatform.system};
          desktop = {
            type = "app";
            program = "${packages.fieldcraft}/bin/fieldcraft";
          };
          cli = {
            type = "app";
            program = "${packages.fieldcraft-cli}/bin/fieldcraft-cli";
          };
        in
        {
          default = desktop;
          fieldcraft = desktop;
          inherit cli;
          fieldcraft-cli = cli;
        }
      );

      formatter = forAllSystems (pkgs: pkgs.nixfmt);
    };
}
