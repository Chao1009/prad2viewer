# Source this file to set up the PRad2 environment (bash/zsh).
#   source <prefix>/bin/setup.sh

PRAD2_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

export PATH="${PRAD2_DIR}/bin${PATH:+:$PATH}"
export LD_LIBRARY_PATH="${PRAD2_DIR}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PRAD2_DATABASE_DIR="${PRAD2_DIR}/database"
export PRAD2_RESOURCE_DIR="${PRAD2_DIR}/resources"
