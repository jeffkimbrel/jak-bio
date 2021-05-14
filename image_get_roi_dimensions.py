#!/usr/bin/env python

import sys
import os
import pandas as pd
import numpy as np
import argparse
import json
import datetime

from jakomics import colors, utilities

from skimage import exposure
import matplotlib.pyplot as plt
from skimage.filters import gaussian, threshold_local
from skimage.transform import rescale

import napari
import skimage.data
import skimage.filters
from napari.layers import Image
from napari.types import ImageData

from magicgui import magicgui

from skimage import (
    color, feature, io, measure, morphology, segmentation, util
)

import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Return ROI Information for a list or directory of images',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with fasta files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual fasta files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('-o',
                    '--out',
                    help="File to write results to",
                    default="ROIs.txt",
                    required=False)

parser.add_argument('-j',
                    '--json',
                    help="Use parameters from a file",
                    default=None,
                    required=False)


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
                    default=2,
                    type=int)

parser.add_argument('-bt', '--background_treshold',
                    help="Background threshold",
                    required=False,
                    default=0.8,
                    type=float)

parser.add_argument('--low',
                    help="Low",
                    required=False,
                    default=0.2,
                    type=float)

parser.add_argument('--high',
                    help="High",
                    required=False,
                    default=99.8,
                    type=float)

parser.add_argument('-at', '--area_treshold',
                    help="Area threshold (in pixels)",
                    required=False,
                    default=15,
                    type=int)

parser.add_argument('-ba', '--blur_amount',
                    help="Blur amount",
                    required=False,
                    default=1,
                    type=int)

parser.add_argument('-bf', '--blur_factor',
                    help="Blur factor",
                    required=False,
                    default=3,
                    type=int)

parser.add_argument('-bs', '--block_size',
                    help="Block size",
                    required=False,
                    default=21,
                    type=int)

parser.add_argument('--interactive',
                    '-i',
                    action='store_true',
                    help='Interactively pick parameters')

args = parser.parse_args()
args.crop_lrtb = args.crop_lrtb.split(",")

if args.json != None:
    with open(args.json) as json_file:
        params = json.load(json_file)

    print(f'{colors.bcolors.YELLOW}Overwriting parameters with {args.json}: {colors.bcolors.END}')
    print(
        f"{colors.bcolors.YELLOW}crop_lrtb: {args.crop_lrtb} -> {params['crop_lrtb']}{colors.bcolors.END}")
    print(
        f"{colors.bcolors.YELLOW}scale: {args.scale} -> {params['scale']}{colors.bcolors.END}")
    print(
        f"{colors.bcolors.YELLOW}min_distance: {args.min_distance} -> {params['min_distance']}{colors.bcolors.END}")
    print(
        f"{colors.bcolors.YELLOW}background_treshold: {args.background_treshold} -> {params['background_treshold']}{colors.bcolors.END}")
    print(
        f"{colors.bcolors.YELLOW}rescale low: {args.low} -> {params['low']}{colors.bcolors.END}")
    print(
        f"{colors.bcolors.YELLOW}rescale high: {args.high} -> {params['high']}{colors.bcolors.END}")
    print(
        f"{colors.bcolors.YELLOW}area_treshold: {args.area_treshold} -> {params['area_treshold']}{colors.bcolors.END}")
    print(
        f"{colors.bcolors.YELLOW}blur_amount: {args.blur_amount} -> {params['blur_amount']}{colors.bcolors.END}")
    print(
        f"{colors.bcolors.YELLOW}blur_factor: {args.blur_factor} -> {params['blur_factor']}{colors.bcolors.END}")
    print(
        f"{colors.bcolors.YELLOW}block_size: {args.block_size} -> {params['block_size']}{colors.bcolors.END}")

    args.crop_lrtb = params['crop_lrtb']
    args.scale = params['scale']
    args.min_distance = params['min_distance']
    args.background_treshold = params['background_treshold']
    args.low = params['low']
    args.high = params['high']
    args.area_treshold = params['area_treshold']
    args.blur_amount = params['blur_amount']
    args.blur_factor = params['blur_factor']
    args.block_size = params['block_size']

# FUNCTIONS ###################################################################


def blur(image, factor=4, sigma=1, mode='nearest'):
    g = rescale(image, factor, anti_aliasing=False)
    g = gaussian(g, sigma=sigma, mode=mode)
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


def segment(i, block_size, distance):
    img_threshold = threshold_local(i,
                                    block_size,
                                    method='mean',
                                    mode='reflect',
                                    cval=0)

    img_binary = img_gaussian > img_threshold

    local_maxi = feature.peak_local_max(i,
                                        indices=False,
                                        min_distance=distance)

    markers = measure.label(local_maxi)
    img_watershed = segmentation.watershed(-i,
                                           markers,
                                           mask=img_binary,
                                           watershed_line=True)

    return img_watershed


# MAIN ########################################################################

if __name__ == "__main__":

    jak_utils.header()

    cols = ['label', 'area', 'perimeter', 'y', 'x', 'length', 'width',
            'scaled_length', 'scaled_width', 'scaled_perimeter', 'scaled_area', 'Image']

    results = pd.DataFrame(columns=cols)

    images = utilities.get_files(args.files, args.in_dir, ['tif', 'tiff', 'jpg'])

    if args.interactive and len(images) > 1:
        print(f'{colors.bcolors.RED}Please pick only one image when using interactive mode{colors.bcolors.END}', file=sys.stderr)
        sys.exit()

    for file in images:

        file.view()

        img_orig = io.imread(file.file_path, as_gray=True)
        img_cropped = img_orig[int(args.crop_lrtb[2]): img_orig.shape[0]-int(args.crop_lrtb[3]),
                               int(args.crop_lrtb[0]): img_orig.shape[1]-int(args.crop_lrtb[1])]  # y axis, x-axis

    #####################################################################################

        if args.interactive:
            with napari.gui_qt():
                viewer = napari.Viewer()
                viewer.add_image(img_cropped, name="Cropped Original", blending='additive')

                @magicgui(
                    auto_call=True,
                    background={"widget_type": "FloatSlider", "max": 1, "min": 0},
                    low={"widget_type": "FloatSlider", "max": 10, "min": 0.01},
                    high={"widget_type": "FloatSlider", "max": 100.0, "min": 90.0, },
                )
                def fix_contrast_gui(layer: Image,
                                     low: float = args.low,
                                     high: float = args.high,
                                     background: float = args.background_treshold
                                     ) -> ImageData:

                    if layer:
                        img_background = normalize_range(layer.data)
                        args.background_treshold = background
                        args.low = low
                        args.high = high

                        return(fix_contrast(img_background,
                                            background=args.background_treshold,
                                            low=args.low,
                                            high=args.high))

                viewer.window.add_dock_widget(fix_contrast_gui.native,
                                              area='bottom', name='Background')
                viewer.layers.events.changed.connect(fix_contrast_gui.reset_choices)

        img_background = normalize_range(img_cropped)
        img_background = fix_contrast(img_background,
                                      background=args.background_treshold,
                                      low=args.low,
                                      high=args.high)

    #####################################################################################

        if args.interactive:
            with napari.gui_qt():
                viewer = napari.Viewer()
                viewer.add_image(img_background,
                                 name="Background Adjusted",
                                 blending='additive')
                viewer.add_image(img_cropped,
                                 name="Cropped Original",
                                 blending='additive')

                @ magicgui(
                    auto_call=True,
                    sigma={"widget_type": "SpinBox", "max": 5, "min": 0},
                    factor={"widget_type": "SpinBox", "max": 9, "min": 1},

                )
                def blur_gui(layer: Image,
                             sigma: float = args.blur_amount,
                             factor: float = args.blur_factor) -> ImageData:

                    if layer:
                        args.blur_amount = sigma
                        args.blur_factor = factor
                        return blur(img_background, sigma=args.blur_amount, factor=args.blur_factor)

                viewer.window.add_dock_widget(blur_gui.native, area='bottom', name='Gaussian')
                viewer.layers.events.changed.connect(blur_gui.reset_choices)

        img_gaussian = blur(img_background,
                            sigma=args.blur_amount,
                            factor=args.blur_factor)

    #####################################################################################

        if args.interactive:
            with napari.gui_qt():
                viewer = napari.Viewer()
                viewer.add_image(img_gaussian, name="Gaussian Adjusted", blending='additive')
                viewer.add_image(img_cropped, name="Cropped Original", blending='additive')

                @ magicgui(
                    auto_call=True,
                    block_size={"widget_type": "SpinBox", "max": 35, "min": 1, "step": 2},
                    distance={"widget_type": "SpinBox", "max": 15, "min": 1},

                )
                def segment_gui(layer: Image,
                                block_size: float = args.block_size,
                                distance: float = args.min_distance) -> ImageData:
                    if layer:

                        block_size = int(block_size)
                        if block_size % 2 == 0:
                            block_size += 1

                        distance = int(distance)

                        args.block_size = block_size
                        args.min_distance = distance

                        return segment(img_gaussian, args.block_size, args.min_distance)

                viewer.window.add_dock_widget(
                    segment_gui.native, area='bottom', name='Segmentation')
                viewer.layers.events.changed.connect(segment_gui.reset_choices)

        img_watershed = segment(img_gaussian, args.block_size, args.min_distance)

        # get region properties
        def print_region_props(rp, ax):
            for feature in rp:
                if feature.area > args.area_treshold:
                    ax.text(feature.centroid[1] - (feature.major_axis_length / 2),
                            feature.centroid[0] + (feature.major_axis_length / 2),
                            feature.label,
                            fontsize=2,
                            color="white")

            return(ax)

        props = measure.regionprops(img_watershed)

        entries = []
        for p in props:
            entry = [p['label'], p['area'], p['perimeter'], *p['centroid'],
                     p['major_axis_length'], p['minor_axis_length']]
            entries.append(entry)

        df = pd.DataFrame(entries, columns=['label', 'area',
                                            'perimeter', 'y', 'x', 'length', 'width'])

        df['scaled_length'] = df['length'] * args.scale
        df['scaled_width'] = df['width'] * args.scale
        df['scaled_perimeter'] = df['perimeter'] * args.scale
        df['scaled_area'] = df['area'] * args.scale * args.scale
        df['Image'] = file.name

        df = df[df['area'] >= args.area_treshold]

        results = pd.concat([results, df])

        # plotting
        color_labels = color.label2rgb(img_watershed, img_gaussian, alpha=0.5, bg_label=0)

        fig, ax = plt.subplots()
        ax.imshow(color_labels, cmap=plt.cm.gray)
        ax = print_region_props(props, ax)
        ax.axis('image')
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout()

        plt.savefig(file.file_path + '_ROI_labels.png', dpi=600)

        io.imsave(file.file_path + '_ROIs.png', util.img_as_ubyte(color_labels))

    # write to file with comments
    if os.path.exists(args.out):
        os.remove(args.out)

    f = open(args.out, 'a')
    for c in jak_utils.header(r=True):
        print(f'# {c}', file=f)

    for arg in vars(args):
        print(f'# ARG {arg} = {getattr(args, arg)}', file=f)

    results.to_csv(f, sep="\t", index=False)

    for_json = ['crop_lrtb',
                'scale',
                'min_distance',
                'background_treshold',
                'low',
                'high',
                'area_treshold',
                'blur_amount',
                'blur_factor',
                'block_size']

    j = {}
    for arg in vars(args):
        #print(arg, getattr(args, arg))
        if arg in for_json:
            j[arg] = getattr(args, arg)

    now = datetime.datetime.now()
    if args.interactive:
        j['training_file'] = file.file_path
        json_out = os.path.abspath(file.file_path) + "." + now.strftime("%Y%m%d_%H%M%S") + ".json"
    else:
        json_out = now.strftime("%Y%m%d_%H%M%S") + ".json"

    j['timestamp'] = now.strftime("%Y-%m-%d %H:%M:%S")

    with open(json_out, 'w') as outfile:
        json.dump(j, outfile, indent=2)

    print(f'Parameter file written to {colors.bcolors.YELLOW}{json_out}{colors.bcolors.END}')
