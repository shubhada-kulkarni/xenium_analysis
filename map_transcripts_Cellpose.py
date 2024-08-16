#!/usr/bin/env python3

import argparse
import csv
import os
import sys
import numpy as np
import pandas as pd
import re
import scipy.sparse as sparse
import scipy.io as sio
import subprocess


# Define constant.
# z-slices are 3 microns apart in morphology.ome.tif
Z_SLICE_MICRON = 3


def main():
    # Parse input arguments.
    args = parse_args()

    # Check for existence of input file.
    if (not os.path.exists(args.cellpose)):
        print("The specified Cellpose output (%s) does not exist!" % args.cellpose)
        sys.exit(0)
    if (not os.path.exists(args.transcript)):
        print("The specified transcripts.parquet file (%s) does not exist!" % args.transcript)
        sys.exit(0)

    # Check if output folder already exist.
    if (os.path.exists(args.out)):
        print("The specified output folder (%s) already exists!" % args.out)
        sys.exit(0)

    # Define additional constants
    NUC_EXP_PIXEL = args.nuc_exp / args.pix_size
    NUC_EXP_SLICE = args.nuc_exp / Z_SLICE_MICRON

    # Read Cellpose segmentation mask
    seg_data = np.load(args.cellpose, allow_pickle=True).item()
    mask_array = seg_data['masks']
    # Use regular expression to extract dimensions from mask_array.shape
    m = re.match("\((?P<z_size>\d+), (?P<y_size>\d+), (?P<x_size>\d+)", str(mask_array.shape))
    mask_dims = { key:int(m.groupdict()[key]) for key in m.groupdict() }

    # Read 5 columns from transcripts Parquet file
    transcripts_df = pd.read_parquet(args.transcript,
                                     columns=["feature_name",
                                              "x_location",
                                              "y_location",
                                              "z_location",
                                              "qv"])


    # Find distinct set of features.
    features = np.unique(transcripts_df["feature_name"])

    # Create lookup dictionary
    feature_to_index = dict()
    for index, val in enumerate(features):
        feature_to_index[str(val)] = index

    # Find distinct set of cells. Discard the first entry which is 0 (non-cell)
    cells = np.unique(mask_array)[1:]

    # Create a cells x features data frame, initialized with 0
    matrix = pd.DataFrame(0, index=range(len(features)), columns=cells, dtype=np.int32)


    # Iterate through all transcripts
    for index, row in transcripts_df.iterrows():
        if index % args.rep_int == 0:
            print(index, "transcripts processed.")

        feature = str(row['feature_name'])
        x = row['x_location']
        y = row['y_location']
        z = row['z_location']
        qv = row['qv']

        # Ignore transcript below user-specified cutoff
        if qv < args.qv_cutoff:
            continue

        # Convert transcript locations from physical space to image space
        x_pixel = x / args.pix_size
        y_pixel = y / args.pix_size
        z_slice = z / Z_SLICE_MICRON

        # Add guard rails to make sure lookup falls within image boundaries.
        x_pixel = min(max(0, x_pixel), mask_dims["x_size"] - 1)
        y_pixel = min(max(0, y_pixel), mask_dims["y_size"] - 1)
        z_slice = min(max(0, z_slice), mask_dims["z_size"] - 1)

        # Look up cell_id assigned by Cellpose. Array is in ZYX order.
        cell_id = mask_array[round(z_slice)] [round(y_pixel)] [round(x_pixel)]

        # If cell_id is 0, Cellpose did not assign the pixel to a cell. Need to perform
        # neighborhood search. See if nearest nucleus is within user-specified distance.
        if cell_id == 0:
            # Define neighborhood boundary for 3D ndarray slicing. Take image boundary into
            # consideration to avoid negative index.
            z_neighborhood_min_slice = max(0, round(z_slice-NUC_EXP_SLICE))
            z_neighborhood_max_slice = min(mask_dims["z_size"], round(z_slice+NUC_EXP_SLICE+1))
            y_neighborhood_min_pixel = max(0, round(y_pixel-NUC_EXP_PIXEL))
            y_neighborhood_max_pixel = min(mask_dims["y_size"], round(y_pixel+NUC_EXP_PIXEL+1))
            x_neighborhood_min_pixel = max(0, round(x_pixel-NUC_EXP_PIXEL))
            x_neighborhood_max_pixel = min(mask_dims["x_size"], round(x_pixel+NUC_EXP_PIXEL+1))

            # Call helper function to see if nearest nucleus is within user-specified distance.
            cell_id = nearest_cell(args, x_pixel, y_pixel, z_slice,
                                   x_neighborhood_min_pixel,
                                   y_neighborhood_min_pixel,
                                   z_neighborhood_min_slice,
                                   mask_array[z_neighborhood_min_slice : z_neighborhood_max_slice,
                                              y_neighborhood_min_pixel : y_neighborhood_max_pixel,
                                              x_neighborhood_min_pixel : x_neighborhood_max_pixel])

        # If cell_id is not 0 at this point, it means the transcript is associated with a cell
        if cell_id != 0:
            # Increment count in feature-cell matrix
            matrix.at[feature_to_index[feature], cell_id] += 1

    # Call a helper function to create Seurat and Scanpy compatible MTX output
    write_sparse_mtx(args, matrix, cells, features)



#--------------------------
# Helper functions

def parse_args():
    """Parses command-line options for main()."""
    summary = 'Map Xenium transcripts to Cellpose segmentation result. \
               Generate Seurat/Scanpy-compatible feature-cell matrix.'

    parser = argparse.ArgumentParser(description=summary)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-cellpose',
                               required = True,
                               help="The path to the *.ome_seg.npy file produced " +
                                    "by Cellpose.")
    requiredNamed.add_argument('-transcript',
                               required = True,
                               help="The path to the transcripts.parquet file produced " +
                                    "by Xenium.")
    requiredNamed.add_argument('-out',
                               required = True,
                               help="The name of output folder in which feature-cell " +
                                    "matrix is written.")
    requiredNamed.add_argument('-pix_size',
                               required = True,
                               type=float,
                               help="The size of each pixel, in microns, for the " +
                                    "image that is passed to Cellpose for nucleus " +
                                    "segmentation.")
    parser.add_argument('-nuc_exp',
                        default='10.0',
                        type=float,
                        help="The expansion distance from the nuclear boundary, " +
                             "in microns, for cell boundary. (default: 10.0)")
    parser.add_argument('-qv_cutoff',
                        default='20.0',
                        type=float,
                        help="Ignore transcripts with QV score below this " +
                             "threshold. (default: 20.0)")
    parser.add_argument('-rep_int',
                        default='10000',
                        type=int,
                        help="Reporting interval. Will print message to stdout " +
                             "whenever the specified # of transcripts is processed. " +
                             "(default: 10000)")

    try:
        opts = parser.parse_args()
    except:
        sys.exit(0)

    return opts


def nearest_cell(args, x_pixel, y_pixel, z_slice,
                 x_neighborhood_min_pixel, y_neighborhood_min_pixel,
                 z_neighborhood_min_slice, mask_array):
    """Check if nearest nucleus is within user-specified distance.
       If function returns 0, it means no suitable nucleus was found."""

    # Initialize constants
    NUC_EXP_PIXEL = args.nuc_exp / args.pix_size
    # For Euclidean distance, we need to convert z-slice to z-micron
    SLICE_TO_PIXEL = Z_SLICE_MICRON / args.pix_size
    # When we take a neighborhood slice of mask_array, all indices start at (0,0,0).
    # This INDEX_SHIFT is necessary to reconstruct coordinates from original mask_array.
    INDEX_SHIFT = np.array([z_neighborhood_min_slice,
                            y_neighborhood_min_pixel,
                            x_neighborhood_min_pixel])

    min_dist = NUC_EXP_PIXEL
    cell_id = 0

    # Enumerate through all points in the neighborhood
    for index, cell in np.ndenumerate(mask_array):
        # Current point is not assigned to a nucleus.
        if cell == 0:
            continue
        # Current point IS assigned to a nucleus. But is it within NUC_EXP_PIXEL?
        else:
            img_loc = np.asarray(index, dtype=float) + INDEX_SHIFT
            # Convert from z-slice to "z-pixel"
            img_loc[0] *= SLICE_TO_PIXEL

            transcript_loc = np.array([z_slice * SLICE_TO_PIXEL, y_pixel, x_pixel])
            # Calculate Euclidean distance between 2 points
            dist = np.linalg.norm(transcript_loc - img_loc)

            if dist < min_dist:
                min_dist = dist
                cell_id = cell

    return cell_id


def write_sparse_mtx(args, matrix, cells, features):
    """Write feature-cell matrix in Seurat/Scanpy-compatible MTX format"""

    # Create the matrix folder.
    os.mkdir(args.out)

    # Convert matrix to scipy's COO sparse matrix.
    sparse_mat = sparse.coo_matrix(matrix.values)

    # Write matrix in MTX format.
    sio.mmwrite(args.out + "/matrix.mtx", sparse_mat)

    # Write cells as barcodes.tsv. File name is chosen to ensure
    # compatibility with Seurat/Scanpy.
    with open(args.out + "/barcodes.tsv", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for cell in cells:
            writer.writerow(["cell_" + str(cell)])

    # Write features as features.tsv. Write 3 columns to ensure
    # compatibility with Seurat/Scanpy.
    with open(args.out + "/features.tsv", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for f in features:
            feature = str(f)
            if feature.startswith("NegControlProbe_") or feature.startswith("antisense_"):
                writer.writerow([feature, feature, "Negative Control Probe"])
            elif feature.startswith("NegControlCodeword_"):
                writer.writerow([feature, feature, "Negative Control Codeword"])
            elif feature.startswith("BLANK_"):
                writer.writerow([feature, feature, "Blank Codeword"])
            else:
                writer.writerow([feature, feature, "Gene Expression"])

    # Seurat expects all 3 files to be gzipped
    subprocess.run("gzip -f " + args.out + "/*", shell=True)



if __name__ == "__main__":
    main()

