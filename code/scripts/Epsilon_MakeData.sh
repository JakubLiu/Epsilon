#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

# -------------------------------
# Defaults
# -------------------------------

OUTPUT_PREFIX="output"
NRANKS=1


# -------------------------------
# Help message
# -------------------------------

usage() {
    echo ""
    echo "Usage:"
    echo "  $0 --bamlist <file> --loci_list <file> --reference_genome <file> --alt_mode <specific|generic> --nranks <N> [options]"
    echo ""
    echo "Required arguments:"
    echo "  --bamlist             BAM list file"
    echo "  --loci_list           BED/loci file"
    echo "  --reference_genome    Reference FASTA file"
    echo "  --alt_mode            specific or generic"
    echo "  --nranks              Number of MPI ranks"
    echo ""
    echo "Optional arguments:"
    echo "  --output_prefix       Output file prefix (default: output)"
    echo ""
    exit 1
}


# -------------------------------
# Parse flags
# -------------------------------

while [[ $# -gt 0 ]]; do

    case $1 in

        --bamlist)
            BAMLIST="$2"
            shift 2
            ;;

        --loci_list)
            LOCI_LIST="$2"
            shift 2
            ;;

        --reference_genome)
            REFERENCE="$2"
            shift 2
            ;;

        --alt_mode)
            ALT_MODE="$2"
            shift 2
            ;;

        --nranks)
            NRANKS="$2"
            shift 2
            ;;

        --output_prefix)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;

        --help|-h)
            usage
            ;;

        *)
            echo "Unknown argument: $1"
            usage
            ;;

    esac

done


# -------------------------------
# Check required arguments
# -------------------------------

if [[ -z "${BAMLIST:-}" ||
      -z "${LOCI_LIST:-}" ||
      -z "${REFERENCE:-}" ||
      -z "${ALT_MODE:-}" ]]; then

    echo "ERROR: Missing required arguments"
    usage

fi


# -------------------------------
# Select python script
# -------------------------------

case "$ALT_MODE" in

    specific)
        # Check whether the loci file contains an alt column
        if head -n 1 "$LOCI_LIST" | grep -qw "alt"; then

            PYTHON_SCRIPT="$SCRIPT_DIR/make_data_alt_specific.py"

        else
            echo " "
            echo " "
            echo "$(date)"
            echo "WARNING: loci file does not contain an 'alt' column."
            echo "Cannot run alt-specific mode."
            echo "Switching to alt-generic mode."
            echo " "
            echo " "

            PYTHON_SCRIPT="$SCRIPT_DIR/make_data_alt_generic.py"
            ALT_MODE="generic"

        fi
        ;;

    generic)
        PYTHON_SCRIPT="$SCRIPT_DIR/make_data_alt_generic.py"
        ;;

    *)
        echo "ERROR: --alt_mode must be 'specific' or 'generic'"
        exit 1
        ;;

esac


# -------------------------------
# Run MPI
# -------------------------------

echo " "
echo " "
echo "--------------------------------"
echo "$(date)"
echo "Running alt extraction"
echo "Mode: $ALT_MODE"
echo "Ranks: $NRANKS"
echo "Script: $PYTHON_SCRIPT"
echo "--------------------------------"
echo " "
echo " "


mpirun -np "$NRANKS" python "$PYTHON_SCRIPT" \
    --bamlist "$BAMLIST" \
    --loci_list "$LOCI_LIST" \
    --reference_genome "$REFERENCE" \
    --output_file_prefix "$OUTPUT_PREFIX"


# ========================================== usage ================================================


# --------------------------------- alt specific mode -------------------------------------------
#./Epsilon_MakeData.sh \
#    --bamlist bamlist.txt \
#    --loci_list loci.bed \
#    --reference_genome hg38.fa \
#    --alt_mode specific \
#    --nranks 16 \
#    --output_prefix sample_specific



# -------------------------------------- alt generic mode --------------------------------------
#./Epsilon_MakeData.sh \
#    --bamlist bamlist.txt \
#    --loci_list loci.bed \
#    --reference_genome hg38.fa \
#    --alt_mode generic \
#    --nranks 32 \
#    --output_prefix sample_generic