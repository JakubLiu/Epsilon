#!/usr/bin/bash

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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"


# defaults
MODE="binomial_test"
METHOD="BH"
ALPHA="0.05"
PRIOR="0.3"
POSTERIOR="0.5"
ALT_MODE="generic"

MODEL_RDS=""
INPUT_DATA=""
OUTPUT_VCF=""


while [[ $# -gt 0 ]]; do
    case "$1" in

        --model)
            MODEL_RDS="$2"
            shift 2
            ;;

        --input)
            INPUT_DATA="$2"
            shift 2
            ;;

        --output)
            OUTPUT_VCF="$2"
            shift 2
            ;;

        --mode)
            MODE="$2"
            shift 2
            ;;

        --fdr_method)
            METHOD="$2"
            shift 2
            ;;

        --alpha)
            ALPHA="$2"
            shift 2
            ;;

        --prior)
            PRIOR="$2"
            shift 2
            ;;

        --posterior_cutoff)
            POSTERIOR="$2"
            shift 2
            ;;

        --alt_mode)
            ALT_MODE="$2"
            shift 2
            ;;

        *)
            echo "Unknown argument: $1"
            exit 1
            ;;

    esac
done


# required arguments
if [[ -z "$MODEL_RDS" ]]; then
    echo "ERROR: --model is required"
    exit 1
fi

if [[ -z "$INPUT_DATA" ]]; then
    echo "ERROR: --input is required"
    exit 1
fi

if [[ -z "$OUTPUT_VCF" ]]; then
    echo "ERROR: --output is required"
    exit 1
fi


if [[ "$MODE" != "binomial_test" && "$MODE" != "bayes_posterior" ]]; then
    echo "ERROR: --mode must be binomial_test or bayes_posterior"
    exit 1
fi


if [[ "$ALT_MODE" != "specific" && "$ALT_MODE" != "generic" ]]; then
    echo "ERROR: --alt_mode must be specific or generic"
    exit 1
fi


echo "=========================================="
echo "Running Epsilon variant caller"
echo "Model:        $MODEL_RDS"
echo "Input:        $INPUT_DATA"
echo "Output VCF:   $OUTPUT_VCF"
echo "Mode:         $MODE"
echo "Alt mode:     $ALT_MODE"
echo "=========================================="


if [[ "$MODE" == "binomial_test" ]]; then

    echo "Correction:   $METHOD"
    echo "Alpha:        $ALPHA"

    Rscript "$PROJECT_DIR/scripts/call.R" \
        "$MODEL_RDS" \
        "$INPUT_DATA" \
        "$OUTPUT_VCF" \
        "$MODE" \
        "$METHOD" \
        "$ALPHA" \
        "$ALT_MODE"


elif [[ "$MODE" == "bayes_posterior" ]]; then

    echo "Prior:             $PRIOR"
    echo "Posterior cutoff:  $POSTERIOR"

    Rscript "$PROJECT_DIR/scripts/call.R" \
        "$MODEL_RDS" \
        "$INPUT_DATA" \
        "$OUTPUT_VCF" \
        "$MODE" \
        "$PRIOR" \
        "$POSTERIOR" \
        "$ALT_MODE"

fi

echo "$(date)"
echo " "
echo " "