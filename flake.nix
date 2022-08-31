{
  inputs = { nixpkgs.url = "github:daanboer/nix-packages"; };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
    in {
      devShell.${system} = pkgs.mkShell {
        name = "lxcat-ng";
        buildInputs = with pkgs.unstable; [
          gcc
          gnumake
          cmake
          clang-tools
          emscripten
          eigen
          nodePackages.prettier
          nodePackages.typescript-language-server
        ];
      };
    };
}
