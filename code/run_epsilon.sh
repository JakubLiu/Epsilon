#!/usr/bin/bash




config=""
jobs=""


while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            config="$2"
            shift 2
            ;;
        --jobs)
            jobs="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

snakemake --snakefile Epsilon.smk --configfile $config --cores $jobs