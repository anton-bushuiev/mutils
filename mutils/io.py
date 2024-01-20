import os
from pathlib import Path

import pylab
import pickle
import dill
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from mutils.definitions import *


def save_list_in_chunks(lst, chunk_size, out_dir_path):
    os.makedirs(out_dir_path, exist_ok=True)
    for chunk_id, i in enumerate(range(0, len(lst), chunk_size)):
        chunk = lst[i:i+chunk_size]
        with open(os.path.join(out_dir_path, str(chunk_id)), 'w') as file:
            file.write('\n'.join(chunk) + '\n')


def read_pickle(pth):
    with open(pth, 'rb') as f:
        return pickle.load(f)


def write_pickle(obj, pth):
    with open(pth, 'wb') as f:
        pickle.dump(obj, f)


def read_dill(pth):
    with open(pth, 'rb') as f:
        return dill.load(f)


def color_generator(n_colors, cmap='hsv'):
    if cmap == 'plotly':  # https://stackoverflow.com/questions/41761654/plotly-where-can-i-find-the-default-color-palette-used-in-plotly-package
        return iter((n_colors // 10 + 1) * [
            (0.38823529411764707, 0.43137254901960786, 0.9803921568627451, 1.0),
            (0.0, 0.8, 0.5882352941176471, 1.0),
            (0.8392156862745098, 0.15294117647058825, 0.1568627450980392, 1.0),
            (0.09019607843137255, 0.7450980392156863, 0.8117647058823529, 1.0),
            (0.7372549019607844, 0.7411764705882353, 0.13333333333333333, 1.0),
            (0.5803921568627451, 0.403921568627451, 0.7411764705882353, 1.0),
            (0.5490196078431373, 0.33725490196078434, 0.29411764705882354, 1.0),
            (0.8901960784313725, 0.4666666666666667, 0.7607843137254902, 1.0),
            (0.4980392156862745, 0.4980392156862745, 0.4980392156862745, 1.0),
            (0.9372549019607843, 0.3333333333333333, 0.23137254901960785, 1.0),
            (0.0859375, 0.46484375, 0.3671875, 1.0)
        ][:n_colors])
    return (pylab.get_cmap(cmap)(1. * i / n_colors) for i in range(n_colors))


def rgb_to_hex(r, g, b):
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)


def get_plotly_palette(reversed_order=False, as_hex=False):
    palette = list(color_generator(11, cmap='plotly'))
    if as_hex:
        palette = [rgb_to_hex(int(256 * c[0]), int(256 * c[1]), int(256 * c[2])) for c in palette]
    if reversed_order:
        return list(reversed(palette))
    return palette


def init_plotting(figsize=(6, 2), font_scale=1.3):
    # Set default figure size
    plt.show()  # Does not work without this line for some reason
    sns.set(rc={'figure.figsize': figsize})
    # Set default style and  font scale
    sns.set_style('whitegrid')
    sns.set_context('paper', font_scale=font_scale)
    sns.set_palette(get_plotly_palette())
    # mpl.rcParams['font.family'] = 'Verdana'
    # mpl.rcParams['grid.color'] = 'gainsboro'


def savefig(name, path, extension='.pdf'):
    if isinstance(path, str):
        path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    path = path / name
    path = path.with_suffix(extension)

    # plt.savefig(path, bbox_inches='tight', pad_inches=0.05)
    plt.savefig(path, bbox_inches='tight')


def append_to_stem(pth: Path, s, sep='_'):
    """
    path/to/file.txt -> path/to/file{sep}{s}.txt
    """
    return pth.parents[0] / f'{pth.stem}{sep}{s}{pth.suffix}'
