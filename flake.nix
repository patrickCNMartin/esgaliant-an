{
  description = "Esgaliant Analysis";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
      in {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            pkgs.python3
            pkgs.uv
            pkgs.git
            pkgs.nextflow
          ];

          shellHook = ''
            echo "Activating or creating Python .venv with uv..."

            if [ ! -d .venv ]; then
              echo "Creating virtual environment using uv..."
              uv venv
              echo "Created .venv with uv"
            fi

            source .venv/bin/activate

            if ! grep -q ".venv/" .gitignore 2>/dev/null; then
              echo ".venv/" >> .gitignore
              echo "Added .venv/ to .gitignore"
            fi

            if [ ! -f pyproject.toml ]; then
              echo "No pyproject.toml found. Initializing with uv..."
              uv init ./
            else
              echo "Found pyproject.toml. Installing dependencies with uv..."
              uv sync
            fi
          '';
        };
      }
    );
}