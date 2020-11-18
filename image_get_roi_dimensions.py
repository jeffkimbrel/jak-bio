#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
import argparse
from jakomics import colors, utilities

from skimage import exposure
import matplotlib.pyplot as plt
from skimage.filters import gaussian, threshold_local
from skimage.transform import rescale
from matplotlib.widgets import Slider
from magicgui import magicgui
from magicgui._qt.widgets import QDoubleSlider
import napari
from napari.layers import Image

from skimage import (
    color, feature, filters, io, measure, morphology, segmentation, util
)

import jak_utils
jak_utils.header()

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='X',
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

parser.add_argument('-at', '--area_treshold',
                    help="Area threshold (in pixels)",
                    required=False,
                    default=15,
                    type=int)

parser.add_argument('-ba', '--blur_amount',
                    help="Blur amount",
                    required=False,
                    default=2,
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

# MAIN ########################################################################


images = utilities.get_files(args.files, args.in_dir, ['tif', 'tiff', 'jpg'])

if args.interactive and len(images) > 1:
    print(f'{colors.bcolors.RED}Please pick only one image when using interactive mode{colors.bcolors.END}', file=sys.stderr)
    sys.exit()

for file in images:

    file.view()

    img_orig = io.imread(file.file_path, as_gray=True)
    img_cropped = img_orig[int(args.crop_lrtb[2]):img_orig.shape[0]-int(args.crop_lrtb[3]),
                           int(args.crop_lrtb[0]):img_orig.shape[1]-int(args.crop_lrtb[1])]  # y axis, x-axis


#####################################################################################

    if args.interactive:
        with napari.gui_qt():
            # create a viewer and add some images
            viewer = napari.Viewer()
            viewer.add_image(img_cropped, name="Cropped Original", blending='additive')

            @magicgui(
                auto_call=True,
                background={"widget_type": QDoubleSlider, "maximum": 1, "minimum": 0},
                low={"maximum": 0.1, "minimum": 0.01},
                high={"minimum": 99.5, "maximum": 100},
            )
            def fix_contrast_gui(layer: Image, background: float = args.background_treshold, low: float = 0.2, high: float = 99.8) -> Image:
                if layer:
                    img_background = normalize_range(img_cropped)
                    args.background_treshold = background
                    return fix_contrast(img_background, background=background, low=low, high=high)

            gui = fix_contrast_gui.Gui()
            viewer.window.add_dock_widget(gui)
            viewer.layers.events.changed.connect(lambda x: gui.refresh_choices("layer"))

    img_background = normalize_range(img_cropped)
    img_background = fix_contrast(img_background, background=args.background_treshold)


#####################################################################################

    if args.interactive:
        with napari.gui_qt():
            viewer = napari.Viewer()
            viewer.add_image(img_background, name="Background Adjusted", blending='additive')
            viewer.add_image(img_cropped, name="Cropped Original", blending='additive')

            @magicgui(
                auto_call=True,
                sigma={"maximum": 5, "minimum": 0},
                factor={"maximum": 9, "minimum": 1},

            )
            def blur_gui(layer: Image, sigma: float = args.blur_amount, factor: float = args.blur_factor) -> Image:
                if layer:
                    args.blur_amount = sigma
                    args.blur_factor = factor
                    return blur(img_background, sigma=sigma, factor=factor)

            gui = blur_gui.Gui()
            viewer.window.add_dock_widget(gui)
            viewer.layers.events.changed.connect(lambda x: gui.refresh_choices("layer"))

    img_gaussian = blur(img_background, sigma=args.blur_amount, factor=args.blur_factor)


#####################################################################################

    if args.interactive:
        with napari.gui_qt():
            viewer = napari.Viewer()
            viewer.add_image(img_gaussian, name="Gaussian Adjusted", blending='additive')
            viewer.add_image(img_cropped, name="Cropped Original", blending='additive')

            @magicgui(
                auto_call=True,
                block_size={"maximum": 35, "minimum": 1},
                distance={"maximum": 15, "minimum": 1},

            )
            def segment_gui(layer: Image, block_size: float = args.block_size, distance: float = args.min_distance) -> Image:
                if layer:

                    block_size = int(block_size)
                    if block_size % 2 == 0:
                        block_size += 1

                    distance = int(distance)

                    args.block_size = block_size
                    args.min_distance = distance

                    # feature detection
                    img_threshold = threshold_local(img_gaussian,
                                                    block_size,
                                                    method='mean',
                                                    mode='reflect',
                                                    cval=0)

                    img_binary = img_gaussian > img_threshold

                    local_maxi = feature.peak_local_max(img_gaussian,
                                                        indices=False,
                                                        min_distance=distance)

                    markers = measure.label(local_maxi)
                    img_watershed = segmentation.watershed(-img_gaussian,
                                                           markers,
                                                           mask=img_binary,
                                                           watershed_line=True)

                    return img_watershed

            gui = segment_gui.Gui()
            viewer.window.add_dock_widget(gui)
            viewer.layers.events.changed.connect(lambda x: gui.refresh_choices("layer"))

    img_threshold = threshold_local(img_gaussian,
                                    args.block_size,
                                    method='mean',
                                    mode='reflect',
                                    cval=0)

    img_binary = img_gaussian > img_threshold

    local_maxi = feature.peak_local_max(img_gaussian,
                                        indices=False,
                                        min_distance=args.min_distance)

    markers = measure.label(local_maxi)
    img_watershed = segmentation.watershed(-img_gaussian,
                                           markers,
                                           mask=img_binary,
                                           watershed_line=True)

    # get region properties

    props = measure.regionprops(img_watershed)

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

    color_labels = color.label2rgb(img_watershed, img_gaussian, alpha=0.5, bg_label=0)

    # if args.visual:
    #     fig, ax = plt.subplots(ncols=2, figsize=(10, 10), sharex=True, sharey=True)
    #     ax[0].imshow(img_cropped)
    #     ax[0].set_title('Original Image')
    #     ax[0].axis('off')
    #     ax[1].imshow(color_labels)
    #     ax[1].set_title('Final Features')
    #
    #     # for id, row in df.iterrows():
    #     #     ax[0].annotate(f"A={int(row['area'])}\nW={int(row['width'])}\nL={int(row['length'])}", xy=(
    #     #         row['x'], row['y']), fontsize=8, color="white")
    #
    #     for id, row in df.iterrows():
    #         ax[1].annotate(f"{int(row['label'])}", xy=(
    #             row['x'], row['y']), fontsize=8, color="white")
    #
    #     plt.show()

    df['Image'] = file.name

    # print(df)
    df.to_csv(file.file_path + "_ROIs.txt", sep="\t", index=False)
    io.imsave(file.file_path + '_ROIs.png', util.img_as_ubyte(color_labels))

for arg in vars(args):
    print(arg, getattr(args, arg))
