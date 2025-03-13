{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
  };

  outputs =
    {
      self,
      nixpkgs,
    }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};

      gccEnv = pkgs.gcc11Stdenv;

      latex =
        with pkgs;
        (texlive.combine {
          inherit (texlive)
            scheme-small
            amsmath
            hyperref
            mhchem
            latexmk
            enumitem
            titlesec
            texcount
            framed
            siunitx
            upquote
            cm-super
            ;
        });

      shellInputs = with pkgs; [
        # Build system
        ninja
        cmake

        # Dependencies
        eigen
        nlohmann_json
        openblas

        # Documentation
        doxygen
        graphviz

        # Coverage
        gcovr
      ];
    in
    {
      devShells.${system} = rec {
        dev = gccEnv.mkDerivation {
          name = "loki-b-dev";
          buildInputs =
            shellInputs
            ++ (with pkgs; [
              # Development tools
              clang-tools

              # Needed to make clang-tools not complain.
              llvmPackages.openmp

              # WebAssembly compilation
              emscripten

              # Plotting
              gnuplot

              # Workflow testing
              act

              # LaTeX
              latex
            ]);
          EM_CACHE = "/tmp/emscripten_cache";
        };
        ci = gccEnv.mkDerivation {
          name = "loki-b-ci";
          buildInputs = shellInputs;
        };
        default = dev;
      };

      packages.${system} = rec {
        loki-b = gccEnv.mkDerivation {
          pname = "loki-b";
          version = "0.0.1";

          src = ./.;

          nativeBuildInputs = with pkgs; [
            cmake
            ninja
            nlohmann_json
            eigen
          ];

          cmakeFlags = [
            "-DENABLE_INSTALL=ON"
            "-DCMAKE_SKIP_INSTALL_ALL_DEPENDENCY=ON"
          ];
          ninjaFlags = [
            "loki"
            "loki_legacytojson"
            "loki_offsidetojson"
          ];
        };
        loki-web = gccEnv.mkDerivation {
          pname = "loki-web";
          version = "0.0.1";

          src = ./.;

          nativeBuildInputs = with pkgs; [
            cmake
            ninja
            eigen
            nlohmann_json
            emscripten
            python3
          ];

          EM_CACHE = "/tmp/emscripten_cache";

          ninjaFlags = [
            "-C build"
            "loki_bindings"
          ];

          configurePhase = ''
            # Create the cache directory.
            mkdir -p $EM_CACHE

            # # Copy the prebuilt emscripten cache.
            cp -r ${pkgs.emscripten}/share/emscripten/cache/* $EM_CACHE

            # Set correct permissions for cache.
            chmod u+rwX -R $EM_CACHE

            # Configure using emscripten.
            emcmake cmake \
              -GNinja \
              -DEigen3_DIR=${pkgs.eigen}/share/eigen3/cmake \
              -Dnlohmann_json_DIR=${pkgs.nlohmann_json}/share/cmake/nlohmann_json \
              -DLOKIB_USE_OPENMP=OFF \
              -B build
          '';

          installPhase = ''
            mkdir -p $out/share/loki-web
            mv web/*.{html,js,wasm,data} $out/share/loki-web
          '';
        };
        coverage = gccEnv.mkDerivation {
          pname = "loki-b-coverage";
          version = "0.0.1";

          src = ./.;

          nativeBuildInputs = with pkgs; [
            cmake
            ninja
            nlohmann_json
            eigen
            gcovr
          ];

          cmakeFlags = [
            "-Dcoverage=ON"
            "-DCMAKE_BUILD_TYPE=Debug"
          ];
          ninjaFlags = [ "coverage" ];

          installPhase = ''
            mkdir $out
            mv coverage.xml $out/coverage.xml
          '';
        };
        notes = gccEnv.mkDerivation {
          pname = "loki-b-notes";
          version = "0.0.1";

          src = ./.;

          nativeBuildInputs = with pkgs; [
            cmake
            ninja
            nlohmann_json
            eigen
            latex
          ];

          cmakeFlags = [ "-Dwith-doc=ON" ];
          ninjaFlags = [ "lokib_doc_pdf" ];

          installPhase = ''
            mkdir $out
            mv doc/notes/cpp_notes.pdf $out
          '';
        };
        default = loki-b;
      };

      apps.${system} = rec {
        loki-b = {
          type = "app";
          program = "${self.packages.${system}.loki-b}/bin/loki";
        };
        default = loki-b;
      };
    };
}
