#!/usr/bin/python


import pysam
import pandas as pd
import numpy as np
from mpi4py import MPI
import os



def gc_fraction(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def split_chunks(lst, n):
    return np.array_split(lst, n)


def merge_rank_outputs(outdir, output_prefix, n_ranks, delete_temp=True):

    dfs = []

    for r in range(n_ranks):

        infile = os.path.join(
            outdir,
            f"{output_prefix}.rank{r}.txt"
        )

        if os.path.exists(infile):
            print(f"Reading {infile}", flush=True)
            dfs.append(pd.read_csv(infile))
        else:
            print(f"Warning: {infile} does not exist.", flush=True)

    if len(dfs) == 0:
        print("No output files found to merge.", flush=True)
        return

    merged_df = pd.concat(dfs, ignore_index=True)

    merged_file = os.path.join(
        outdir,
        f"{output_prefix}.txt"
    )

    merged_df.to_csv(merged_file, index=False)

    print(f"Merged output written to {merged_file}", flush=True)

    if delete_temp:
        for r in range(n_ranks):

            infile = os.path.join(
                outdir,
                f"{output_prefix}.rank{r}.txt"
            )

            if os.path.exists(infile):
                os.remove(infile)

        print("Temporary rank files deleted.", flush=True)



def extract_features_alt_generic(bamlist, reference_genome, loci_list, output_file_prefix="output", outdir="negative_control_data", window_size=100):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()


        if rank == 0:

            loci_df = pd.read_csv(loci_list)  # read the file with all the loci (doesnt have to be actionable)

            # REQUIRED COLUMNS:
            # chrom,pos,ref
            loci = loci_df[["chrom", "pos", "ref"]].values.tolist()   # turn that into a list

            with open(bamlist) as f:
                bamlist = [x.strip() for x in f]

            chunks = split_chunks(loci, size)  # use your function to split it into chunks, one chunk per worker

        else:
            bamlist = None
            chunks = None

        bamlist = comm.bcast(bamlist, root=0)  # broadcast the bamlist to all workers, because the iteration over the bams will be done sequentailly
        local_loci = comm.scatter(chunks, root=0)   # the iterations over the loci will be done in parallel, one chunk per worker



        # list of filehandles to open (each worker sees all), this is kinda already opening of the files
        bam_handles = [
            pysam.AlignmentFile(bam, "rb")
            for bam in bamlist
        ]

        fasta = pysam.FastaFile(reference_genome)


        results = []

        # this loop is done in parallel
        for idx, (chrom, pos, ref) in enumerate(local_loci):

            chrom = str(chrom)
            pos = int(pos)
            ref = ref.upper()

            print(f"[rank {rank}] processing {chrom}:{pos}", flush=True)

            

            # get the start, end indices for the current worker
            start = max(0, pos - window_size // 2)
            end = pos + window_size // 2

            context = fasta.fetch(chrom, start, end).upper()

            center = len(context) // 2

            GGC_upstream = int("GGC" in context[:center])
            CGG_downstream = int("CGG" in context[center:])
            GGT_upstream = int("GGT" in context[:center])
            TGG_downstream = int("TGG" in context[center:])

            GC_content = gc_fraction(context)

            

            # this is done serially (all worker sees all bams)
            for bam_path, bam in zip(bamlist, bam_handles):

                sample_id = os.path.basename(bam_path)

                # set the variables for the current bam
                ref_count = 0
                alt_count_generic = 0
                alt_A = 0
                alt_C = 0
                alt_G = 0
                alt_T = 0

                mapqs = []   # this is a list of all mapqs of all reads that have mapped to the current pos in the current bam
                baseqs = []   # this is a list of all baseqs at positions of reads that have mapped to the current pos in the current bam

                mismatch_values = [] # a list for the current position for the current bam

                ref_R1 = 0
                ref_R2 = 0
                alt_R1 = 0
                alt_R2 = 0


                for read in bam.fetch(chrom, pos - 1, pos):


                    if read.is_unmapped:
                        continue

                    if read.is_duplicate:
                        continue

                    if read.mapping_quality == 0:
                        continue

                    aligned_pairs = read.get_aligned_pairs(matches_only=False)

                    query_pos = None

                    for qpos, rpos in aligned_pairs:

                        if rpos == pos - 1:
                            query_pos = qpos
                            break

                    if query_pos is None:
                        continue

                    if query_pos >= len(read.query_sequence):
                        continue

                    base = read.query_sequence[query_pos].upper()

                    if base == "N":
                        continue

                    baseq = read.query_qualities[query_pos]
                    mapq = read.mapping_quality

                    mapqs.append(mapq)
                    baseqs.append(baseq)

                    try:
                        nm = read.get_tag("NM")
                    except KeyError:
                        nm = 0

                    mismatch_values.append(nm)

                    strand = "-" if read.is_reverse else "+"

                    

                    if base == ref:

                        ref_count += 1

                        if strand == "+":
                            ref_R1 += 1
                        else:
                            ref_R2 += 1

                    else:

                        alt_count_generic += 1

                        if base == "A":
                            alt_A += 1
                        elif base == "C":
                            alt_C += 1
                        elif base == "G":
                            alt_G += 1
                        elif base == "T":
                            alt_T += 1

            

                coverage = ref_count + alt_count_generic

                if coverage == 0:
                    continue

                # for the given position for the given bam compute the medians/means of the lists
                median_mapQ = np.median(mapqs) if mapqs else np.nan
                median_baseQ = np.median(baseqs) if baseqs else np.nan

                mismatches = np.mean(
                    np.array(mismatch_values) > 0
                )

            

                row = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "sampleID": sample_id,
                    "general_alt_counts" : alt_count_generic,
                    "ref_counts": ref_count,
                    "alt_A": alt_A,
                    "alt_C": alt_C,
                    "alt_G": alt_G,
                    "alt_T": alt_T,

                    "coverage": coverage,

                    "mismatches": mismatches,

                    "median_mapQ": median_mapQ,
                    "median_baseQ": median_baseQ,

                    "GGC_upstream": GGC_upstream,
                    "CGG_downstream": CGG_downstream,
                    "GGT_upstream": GGT_upstream,
                    "TGG_downstream": TGG_downstream,

                    "GC_content": GC_content,

                    "ref_R1": ref_R1,
                    "ref_R2": ref_R2,
                    "alt_R1": alt_R1,
                    "alt_R2": alt_R2
                }

                results.append(row)


        os.makedirs(outdir, exist_ok=True)

        local_df = pd.DataFrame(results)

        outfile = os.path.join(
            outdir,
            f"{output_file_prefix}.rank{rank}.txt"
        )

        local_df.to_csv(outfile, index = False)

        print(f"[rank {rank}] wrote {outfile}", flush=True)


        for bam in bam_handles:
            bam.close()

        fasta.close()

        comm.Barrier()

        if rank == 0:
            merge_rank_outputs(
                outdir=outdir,
                output_prefix=output_file_prefix,
                n_ranks=size,
                delete_temp=True
            )

            print("All ranks completed.", flush=True)
            print('\n', flush = True)
            print('\n', flush = True)