#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from jakomics import colors, utilities

from skimage import exposure
import matplotlib.pyplot as plt
from skimage.filters import gaussian, threshold_local
from skimage.transform import rescale

from skimage import (
    color, feature, filters, io, measure, morphology, segmentation, util
)

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='XXX')

parser.add_argument('--in_dir',
                    help="Directory with fasta files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual fasta files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--crop_lrtb',
                    help="Comma-separated list of pixels to crop from ",
                    default="0,0,0,0",
                    required=False)

parser.add_argument('-s', '--scale',
                    help="Scale (microns per pixel)",
                    required=False,
                    default=0.0921,
                    type=float)

parser.add_argument('-md', '--min_distance',
                    help="Minimum distance between peaks",
                    required=False,
                    default=3,
                    type=int)

parser.add_argument('-bt', '--background_treshold',
                    help="Background threshold",
                    required=False,
                    default=0.8,
                    type=float)

parser.add_argument('-at', '--area_treshold',
                    help="Area threshold (in pixels)",
                    required=False,
                    default=15,
                    type=int)

parser.add_argument('-ba', '--blur_amount',
                    help="Blur amount",
                    required=False,
                    default=3,
                    type=int)

parser.add_argument('-bf', '--blur_factor',
                    help="Blur factor",
                    required=False,
                    default=5,
                    type=int)

parser.add_argument('-bs', '--block_size',
                    help="Block size",
                    required=False,
                    default=15,
                    type=int)

parser.add_argument('--visual',
                    '-v',
                    action='store_true',
                    help='Display images')

args = parser.parse_args()
args.crop_lrtb = args.crop_lrtb.split(",")

# FUNCTIONS ###################################################################


def blur(image, factor=4, sigma=1):
    g = rescale(image, factor, anti_aliasing=False)
    g = gaussian(g, sigma=sigma, mode='nearest')
    g = rescale(g, 1/factor, anti_aliasing=False)

    return(g)


def normalize_range(image):
    return (image - image.min()) / (image.max() - image.min())


def fix_contrast(image, low=0.2, high=99.8, background=0.8):
    v_min, v_max = np.percentile(image, (low, high))
    rescaled = exposure.rescale_intensity(image, in_range=(v_min, v_max))
    rescaled[rescaled <= background] = 0

    return rescaled


def subtract_background(image, radius=50, light_bg=False):
    # https://github.com/scikit-image/scikit-image/issues/3538
    from skimage.morphology import white_tophat, black_tophat, disk
    # you can also use 'ball' here to get a slightly smoother result at the cost of increased computing time
    str_el = disk(radius)
    if light_bg:
        return black_tophat(image, str_el)
    else:
        return white_tophat(image, str_el)

# MAIN ########################################################################


images = utilities.get_files(args.files, args.in_dir, ['tif', 'tiff'])

for file in images:

    file.view()

    image_orig = io.imread(file.file_path, as_gray=True)
    image_cropped = image_orig[int(args.crop_lrtb[2]):image_orig.shape[0]-int(args.crop_lrtb[3]),
                               int(args.crop_lrtb[0]):image_orig.shape[1]-int(args.crop_lrtb[1])]  # y axis, x-axis

    # print(image_orig.shape)
    # print(image_cropped.shape)

    if args.visual:
        fig, ax = plt.subplots(ncols=2, sharex=False, sharey=False)
        ax[0].imshow(io.imread(file.file_path))
        ax[0].set_title('original')
        ax[0].axis('off')
        ax[1].imshow(image_cropped)
        ax[1].set_title('cropped')
        ax[1].axis('off')
        plt.show()

    g = normalize_range(image_cropped)
    # g_tophat = subtract_background(g)
    g = fix_contrast(g, background=args.background_treshold)
    #
    # fig, ax = plt.subplots(ncols=3, sharex=True, sharey=True)
    # ax[0].imshow(image_cropped)
    # ax[0].set_title('original')
    # ax[0].axis('off')
    # ax[1].imshow(g)
    # ax[1].set_title('fixed contrast')
    # ax[1].axis('off')
    # ax[2].imshow(g_tophat)
    # ax[2].set_title('Tophat')
    # ax[2].axis('off')
    #
    # plt.show()

    better_contrast = blur(g, sigma=args.blur_amount, factor=args.blur_factor)

    # feature detection

    adaptive_thresh = threshold_local(better_contrast,
                                      args.block_size,
                                      method='mean',
                                      mode='reflect',
                                      cval=0)

    binary_adaptive = better_contrast > adaptive_thresh

    local_maxi = feature.peak_local_max(better_contrast,
                                        indices=False,
                                        min_distance=args.min_distance)

    markers = measure.label(local_maxi)
    segmented_cells_2 = segmentation.watershed(-better_contrast,
                                               markers,
                                               mask=binary_adaptive,
                                               watershed_line=False)

    #

    if args.visual:
        fig, ax = plt.subplots(ncols=3, sharex=True, sharey=True)
        ax[0].imshow(image_cropped)
        ax[0].set_title('original')
        ax[0].axis('off')
        ax[1].imshow(better_contrast)
        ax[1].set_title('background threshold')
        ax[1].axis('off')
        ax[2].imshow(color.label2rgb(segmented_cells_2, bg_label=0))
        ax[2].set_title('Watershed')
        ax[2].axis('off')

        plt.show()

    # get region properties

    props = measure.regionprops(segmented_cells_2)

    entries = []
    for p in props:
        entry = [p['label'], p['area'], p['perimeter'], *p['centroid'],
                 p['major_axis_length'], p['minor_axis_length']]
        entries.append(entry)

    df = pd.DataFrame(entries, columns=['label', 'area', 'perimeter', 'y', 'x', 'length', 'width'])
    df = df[df['area'] >= args.area_treshold]

    df['length'] = df['length'] * args.scale
    df['width'] = df['width'] * args.scale
    df['perimeter'] = df['perimeter'] * args.scale
    df['area'] = df['area'] * args.scale * args.scale

    # plotting

    color_labels = color.label2rgb(segmented_cells_2, better_contrast, alpha=0.5, bg_label=0)

    if args.visual:
        fig, ax = plt.subplots(ncols=2, figsize=(10, 10), sharex=True, sharey=True)
        ax[0].imshow(image_cropped)
        ax[0].set_title('Original Image')
        ax[0].axis('off')
        ax[1].imshow(color_labels)
        ax[1].set_title('Final Features')

        # for id, row in df.iterrows():
        #     ax[0].annotate(f"A={int(row['area'])}\nW={int(row['width'])}\nL={int(row['length'])}", xy=(
        #         row['x'], row['y']), fontsize=8, color="white")

        for id, row in df.iterrows():
            ax[1].annotate(f"{int(row['label'])}", xy=(
                row['x'], row['y']), fontsize=8, color="white")

        plt.show()

    df['Image'] = file.name

    # print(df)
    df.to_csv(file.file_path + "_ROIs.txt", sep="\t", index=False)
    io.imsave(file.file_path + '_ROIs.png', util.img_as_ubyte(color_labels))
