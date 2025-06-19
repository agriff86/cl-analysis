#!/usr/bin/env python

__version__ = "0.1.0"

import sys
import signal

signal.signal(signal.SIGINT, lambda signal, frame: sys.exit(0))
import logging

logging.basicConfig(format="%(name)s: %(message)s")
log = logging.getLogger(sys.argv[0])
import os
import traceback
import re
import itertools
import argparse
import datetime as dt
import numpy as np
import netCDF4
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy import ndimage
import netCDF4
import datetime


#
# ... analysis routines
#


def load_ceilometer_file(
    fname, maxrange=3000, vmin=0, vmax=1, tmin=None, tmax=None, as_xarray=True
):
    """
    Load a Level 1 data file produced by a Vaisala ceilometer

    Arguments
    ---------
        *fname*: string
            The path to the file to load

        *maxrange*: float, default 3000.0
            Maximum range (m) for backscatter.  Data received from further away is excluded.

        *as_xarray*: bool, default True
            If True, return the data as an xarray.  Otherwise, data is returned as Dict of arrays

    Notes
    -----
    Tested on data from CL31
    Loads arrays into memory, don't use on massive files.
    """

    # I want to use Xarray to process lidar data, but XArray can't open
    # this type of netCDF file yet (discussion here:
    # https://github.com/pydata/xarray/issues/2368 )
    # The approach taken is to read using netCDF4, and then create an XArray dataset.
    # WARNING - loads all variables into memory, so don't use this on massive files
    nc = netCDF4.Dataset(fname)
    # The file reports units of "days since 1970-01-01 00:00:00.000"
        # BUT data appear to be in units of "seconds since 1970-01-01 00:00:00.000"
    time_units = "seconds since 1970-01-01 00:00:00.000"
    t = netCDF4.num2date(
        nc.variables["time"][:], time_units, only_use_cftime_datetimes=False
    )
    r = nc.variables["range"][:]
    range_slice = slice(0, np.argmax(r > maxrange))
    r = r[range_slice]
    rcs_910 = nc.variables["rcs_910"][:, range_slice]
    cloudbase = nc.variables["cloud_data"][:, 0]
    cloudbase[nc.variables["cloud_status"][:, 0] == False] = np.nan

    # slice to the time of interest
    if tmax is not None:
        idx_tmax = np.argmax(t >= tmax)
        # 0 means no records need to be removed
        if idx_tmax > 0:
            log.debug(f'Slicing because tmax is not None, [{t[0]},{t[-1]}], t[:{idx_tmax}]')
            rcs_910 = rcs_910[:idx_tmax, :]
            t = t[:idx_tmax]
            cloudbase = cloudbase[:idx_tmax]
    if tmin is not None and len(t) > 0:
        idx_tmin = np.argmax(t >= tmin)
        if idx_tmin > 0:
            rcs_910 = rcs_910[idx_tmin:, :]
            t = t[idx_tmin:]
            log.debug(f'Slicing because tmin is not None, [{t[0]},{t[-1]}] t[{idx_tmin}:]')
            cloudbase = cloudbase[idx_tmin:]
    if len(t) == 0:
        log.debug(f'No records in data')
    else:
        log.debug(f'Time range in data: [{t[0]},{t[-1]}]')

    if as_xarray:
        data = xr.Dataset(
            {
                "cloudbase": ("time", cloudbase),
                "rcs_910": (("time", "range"), rcs_910),
                "range": r,
                "time": t,
            }
        )
    else:
        data = {}
        data["cloudbase"] = cloudbase
        data["rcs_910"] = rcs_910
        data["time"] = t
        data["range"] = r

    return data


def concat_data(data_list):
    # drop zero-length items
    data_list = [itm for itm in data_list if len(itm['time'])>0]
    # TODO: handle non-xarray input
    def getkey(itm):
        return itm['time'].values[0]
    # sort according to the first time in each dataset
    data_list.sort(key=getkey)
    return xr.concat(data_list, dim="time")

def detect_pbl_gradient_method(
    data,
    min_pbl=200,
    max_pbl=2000,
    sigma_smoothing=3,
    pbl_to_cloudbase_min=100,
    take_log_of_rcs=False,
):
    """TODO

    scipy.ndimage.gaussian_filter

    rcs is indexed by rcs[time, range]
    """
    t = data["time"]
    r = data["range"]
    rcs = data["rcs_910"]
    cloudbase = data["cloudbase"]

    # handle xarray inputs by unpacking them to plain arrays
    data_is_xarray = hasattr(t, "values")
    if data_is_xarray:
        t = t.values
        r = r.values
        rcs = rcs.values
        cloudbase = cloudbase.values

    if take_log_of_rcs:
        with np.errstate(all="ignore"):
            rcs = np.log10(rcs)

    rcs_smooth = ndimage.gaussian_filter(rcs, sigma=sigma_smoothing)
    # negative gradient in vertical direction
    rcs_gradient = -np.gradient(rcs_smooth, axis=1)
    # zero out the gradient where
    #  - range is less than the minimum requested
    #  - range is greater than the maximum requested
    #  - PBL is higher than the cloud base
    rcs_gradient[:, r < min_pbl] = 0
    rcs_gradient[:, r > max_pbl] = 0

    idxmax = np.argmax(rcs_gradient, axis=1)
    pbl = r[idxmax]

    rcs_smooth = ndimage.gaussian_filter(rcs, sigma=3)
    # negative gradient in vertical direction
    rcs_gradient = -np.gradient(rcs_smooth, axis=1)
    # zero out the gradient where
    #  - range is less than the minimum requested
    #  - range is greater than the maximum requested
    #  - PBL within some threshold of the cloud base
    rcs_gradient[:, r < min_pbl] = 0
    rcs_gradient[:, r > max_pbl] = 0
    cloudbase_tiled = np.tile(cloudbase[:, np.newaxis], [1, len(r)])
    r_tiled = np.tile(r, [len(t), 1])

    # ignore warning about comparision between NaN and value
    with np.errstate(invalid="ignore"):
        rcs_gradient[r_tiled >= cloudbase_tiled - pbl_to_cloudbase_min] = 0
    # PBL is the largest gradient left, after exclusions
    idxmax = np.argmax(rcs_gradient, axis=1)
    pbl = r[idxmax].astype(float)
    rcs_gradient_at_pblh = np.array(
        [itm[idx] for itm, idx in zip(rcs_gradient, idxmax)]
    )

    if data_is_xarray:
        data["pbl_height"] = (("time",), pbl)
        data["vertical_gradient_at_pblh"] = (("time",), rcs_gradient_at_pblh)
        # include the intermediate fields
        data["rcs_gradient"] = (("time", "range"), rcs_gradient)
        data["rcs_smooth"] = (("time", "range"), rcs_smooth)
    else:
        data["pbl_height"] = pbl
        data["vertical_gradient_at_pblh"] = rcs_gradient_at_pblh
    return data


def draw_plot(data, maxrange=3000, tzoffset=0):
    """
    Draw a plot showing backscatter with diagnosed points overplotted
    """
    t = data["time"].values.copy()
    if tzoffset != 0:
        t += np.timedelta64(tzoffset, 'h')
        time_units_str = f"UTC {tzoffset:+}"
    else:
        time_units_str = "UTC"
    r = data["range"]
    rcs = data["rcs_910"]
    cloudbase = data["cloudbase"]
    pbl = data["pbl_height"]
    rcs_gradient_at_pblh = data["vertical_gradient_at_pblh"]

    # handle xarray inputs by unpacking them to plain arrays
    data_is_xarray = hasattr(t, "values")
    if data_is_xarray:
        t = t.values
        r = r.values
        rcs = rcs.values
        cloudbase = cloudbase.values
        pbl = pbl.values
        rcs_gradient_at_pblh = rcs_gradient_at_pblh.values

    with np.errstate(all="ignore"):
        logrcs = np.log10(rcs)

    # generate the plot
    cmap = mpl.cm.get_cmap("viridis")
    cmap.set_under("#222222")
    cmap.set_bad("#222222")

    # clip logcrs low values (instead of showing empty areas on the plot)

    fig, ax = plt.subplots(figsize=[10, 3.5])
    msh = ax.pcolormesh(t, r, logrcs.T, cmap=cmap)
    ax.plot(
        t,
        cloudbase,
        "s",
        linestyle="None",
        markersize=2,
        markerfacecolor="red",
        markeredgecolor="red",
        label="Cloud base (reported by instrument)",
    )
    ax.plot(
        t,
        pbl,
        ".",
        linestyle="None",
        markersize=2,
        markerfacecolor="blue",
        markeredgecolor="blue",
        label="PBL height (gradient method)",
    )
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.05), ncol=2, frameon=False)

    ax.set_ylim([0, maxrange])
    ax.set_ylabel("Range from lidar (m)")
    ax.set_xlabel(time_units_str)
    fig.colorbar(
        msh,
        label="Log(Range-corrected backscatter) (a.u.)",
        extend="both",
        fraction=0.05,
    )

    return fig


#
# Command line interface
#


def main():
    def validate_timestamp(s):
        try:
            return datetime.datetime.strptime(s, "%Y%m%d%H%M")
        except ValueError:
            msg = "not a valid timestamp: {0!r}".format(s)
            raise argparse.ArgumentTypeError(msg)

    parser = argparse.ArgumentParser(
        description="Detect PBL in Vaisala CL51 NetCDF files and generate plots. "
                    "The name of the output file determines what action is taken: if the "
                    "output ends with .csv, then boundary layer and cloudbase information "
                    "is written to a csv file.  Similarly, if the output has the special name "
                    "'stdout' the same comma-separated data is sent to standard output. "
                    "Plots are generated if the output file ends with .pdf or .png."
    )

#    parser.add_argument(
#        "-q", dest="quiet", action="store_true", help="run quietly (suppress output)"
#    )
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        help="print debugging information",
    )

    parser.add_argument(
        "-st",
        "--start_time",
        help="Start time for display, same timezone as source files, format YYYYMMDDHHMM",
        required=False,
        type=validate_timestamp,
    )
    parser.add_argument(
        "-et",
        "--end_time",
        help="End time for display, same timezone as source files, format YYYYMMDDHHMM",
        required=False,
        type=validate_timestamp,
    )
    parser.add_argument(
        "-tz",
        "--timezone_offset",
        help='Timezone Offset.  Add this many hours to the time in the source files before '
        'displaying plots (e.g. for conversion from UTC to local time).  '
        'NOTE: does not affect CSV output.',
        required=False,
        type=int,
        default=0
    )

    parser.add_argument("input", help="input files", nargs=argparse.ONE_OR_MORE)
    parser.add_argument("output", help="output file")
    args = parser.parse_args()

    print(args)

    if args.debug:
        log.setLevel("DEBUG")

    input_ = args.input
    output = args.output

    if output.lower().endswith(".csv") or output == "stdout":
        generate_plots = False
    elif output.lower().endswith(".pdf") or output.lower().endswith(".png"):
        generate_plots = True
    else:
        print("Output file needs to be .csv, .pdf, .png, or 'stdout'", file=sys.stderr)
        sys.exit(1)

    ds = concat_data([load_ceilometer_file(fn, tmin=args.start_time, tmax=args.end_time) for fn in input_])
    ds = detect_pbl_gradient_method(ds)
    if generate_plots:
        fig = draw_plot(ds, tzoffset=args.timezone_offset)
        fig.savefig(output, bbox_inches="tight")
    else:
        # if we are not generating plots, just dump data to csv
        # (file or stdout)
        if output == "stdout":
            file_descriptor = sys.stdout
        else:
            file_descriptor = open(output, "wt")
        (
            ds[["cloudbase", "pbl_height", "vertical_gradient_at_pblh"]]
            .to_pandas()
            .to_csv(file_descriptor)
        )
        sys.stdout.flush()

    sys.exit(0)

    if os.path.isdir(input_):
        for file_ in sorted([fsencode(x) for x in os.listdir(input_)]):
            if not (file_.endswith(b".dat") or file_.endswith(b".DAT")):
                continue
            input_filename = os.path.join(input_, file_)
            output_filename = os.path.join(output, os.path.splitext(file_)[0] + b".nc")
            if not args.quiet:
                print(fsdecode(input_filename))
            try:
                dd = read_input(input_filename, {"check": args.check})
                if len(dd) > 0:
                    write_output(dd, output_filename)
            except Exception as e:
                log.error(e)
                log.debug(traceback.format_exc())
    else:
        try:
            dd = read_input(input_, {"check": args.check})
            if len(dd) > 0:
                write_output(dd, output)
        except Exception as e:
            log.error(e)
            log.debug(traceback.format_exc())


if __name__ == "__main__":
    main()
