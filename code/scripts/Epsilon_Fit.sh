#!/bin/bash

echo " "
echo " "

cat <<'EOF'

      ███████╗██████╗ ███████╗██╗██╗      ██████╗ ███╗   ██╗
      ██╔════╝██╔══██╗██╔════╝██║██║     ██╔═══██╗████╗  ██║
      █████╗  ██████╔╝███████╗██║██║     ██║   ██║██╔██╗ ██║
      ██╔══╝  ██╔═══╝ ╚════██║██║██║     ██║   ██║██║╚██╗██║
      ███████╗██║     ███████║██║███████╗╚██████╔╝██║ ╚████║
      ╚══════╝╚═╝     ╚══════╝╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═══╝

EOF

echo " "
echo "$(date)"
echo " "

# ==========================================
# Run glmmTMB error model fitting
# ==========================================

# Usage:
# bash run_error_model.sh <input_csv> <noise_level> <specific|generic> <output_model_name>
#
# Example:
# bash run_error_model.sh negative_controls.csv 0.01 specific WES_error_model


# ==========================================
# Command line arguments
# ==========================================

INPUT_CSV=""
NOISE_LEVEL=""
ALT_MODE=""
MODEL_NAME=""

while [[ $# -gt 0 ]]; do
    case "$1" in

        --input)
            INPUT_CSV="$2"
            shift 2
            ;;

        --noise_level)
            NOISE_LEVEL="$2"
            shift 2
            ;;

        --alt_mode)
            ALT_MODE="$2"
            shift 2
            ;;

        --output_model)
            MODEL_NAME="$2"
            shift 2
            ;;

        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 --input <negative_control_csv> --noise_level <value> --alt_mode <specific|generic> --output_model <model_name>"
            exit 1
            ;;

    esac
done


# Required arguments
if [ -z "$INPUT_CSV" ]; then
    echo "ERROR: Missing --input"
    exit 1
fi

if [ -z "$NOISE_LEVEL" ]; then
    echo "ERROR: Missing --noise_level"
    exit 1
fi

if [ -z "$ALT_MODE" ]; then
    echo "ERROR: Missing --alt_mode"
    exit 1
fi

if [ -z "$MODEL_NAME" ]; then
    echo "ERROR: Missing --output_model"
    exit 1
fi

# Check input file
if [ ! -f "$INPUT_CSV" ]; then
    echo "ERROR: Input file does not exist: $INPUT_CSV"
    exit 1
fi

# Check alt mode
if [[ "$ALT_MODE" != "specific" && "$ALT_MODE" != "generic" ]]; then
    echo "ERROR: alt_mode must be 'specific' or 'generic'"
    exit 1
fi

echo " "
echo " "
echo "=========================================="
echo "$(date)"
echo "Running error model fitting"
echo "Input:       $INPUT_CSV"
echo "Noise level: $NOISE_LEVEL"
echo "Alt mode:    $ALT_MODE"
echo "Output:      $MODEL_NAME"
echo "=========================================="
echo " "
echo " "

echo " "
echo " "


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "$SCRIPT_DIR/fit.R" \
    "$INPUT_CSV" \
    "$NOISE_LEVEL" \
    "$ALT_MODE" \
    "$MODEL_NAME"


if [ $? -ne 0 ]; then
    echo "ERROR: R model fitting failed."
    exit 1
fi

echo "$(date)"
echo " "
echo " "