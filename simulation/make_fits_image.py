"""
make_fits_image.py

Takes a waders (psimulCCDimg) output ROOT file, generates a noise matrix,
adds it to the signal, and saves a FITS image and a boolean signal mask.

Outputs:
  <stem>_image.fits  — signal + noise in eV  (for visualisation)
  <stem>_data.npz    — three arrays for clustering validation:
                         'image'  float32 (rows, cols) eV  — signal + noise
                         'signal' float32 (rows, cols) eV  — pure signal
                         'mask'   bool    (rows, cols)     — True = signal pixel

Dependencies: pip install uproot awkward numpy astropy

Usage:
  python make_fits_image.py --input out_simulCCDimg_XXX.root [--outdir .] [--crop r0,r1,c0,c1]
"""

import argparse
import os
import numpy as np
import awkward as ak
from astropy.io import fits
import uproot

# defaults matching config_CCDSensor_PV_208a81z_10_0.json
CCD_ROWS       = 6000
CCD_COLS       = 1500
E2EV           = 3.74          # eV per electron
NOISE_SIGMA    = 0.25          # e-/pixel
NOISE_PEDESTAL = 0.0           # e-/pixel
DARKCURRENT    = 0.001         # e-/pixel/day
EXP_TIME       = 0.3333333333  # days
SEED           = 321


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input",          required=True)
    p.add_argument("--outdir",         default=".")
    p.add_argument("--ccd_rows",       type=int,   default=CCD_ROWS)
    p.add_argument("--ccd_cols",       type=int,   default=CCD_COLS)
    p.add_argument("--e2eV",           type=float, default=E2EV)
    p.add_argument("--noise_sigma",    type=float, default=NOISE_SIGMA)
    p.add_argument("--noise_pedestal", type=float, default=NOISE_PEDESTAL)
    p.add_argument("--darkcurrent",    type=float, default=DARKCURRENT)
    p.add_argument("--exp_time",       type=float, default=EXP_TIME)
    p.add_argument("--seed",           type=int,   default=SEED)
    p.add_argument("--crop",           default=None,
                   help="row_min,row_max,col_min,col_max  (0-indexed, exclusive end)")
    return p.parse_args()


def main():
    args  = parse_args()
    rng   = np.random.default_rng(seed=args.seed)
    shape = (args.ccd_rows, args.ccd_cols)

    os.makedirs(args.outdir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(args.input))[0]

    # ── 1. read waders output, flatten all nesting into 1D arrays ─────────────
    print(f"Reading {args.input}")
    with uproot.open(args.input) as f:
        tree     = f["pixelizedEvent"]
        px_raw   = tree["pixels_x"].array()
        py_raw   = tree["pixels_y"].array()
        edep_raw = tree["pixels_Edep"].array()   # keV

    # flatten completely regardless of nesting depth (event → ccd → pixel)
    px   = ak.to_numpy(ak.flatten(px_raw,   axis=None)).astype(int)
    py   = ak.to_numpy(ak.flatten(py_raw,   axis=None)).astype(int)
    edep = ak.to_numpy(ak.flatten(edep_raw, axis=None)).astype(np.float32) * 1e3  # keV → eV

    print(f"  Total signal pixels across all events: {len(px):,}")

    # ── 2. stamp signal onto CCD grid ─────────────────────────────────────────
    # pixels_x = x_pixel = col index [0, N_x=1500)
    # pixels_y = y_pixel = row index [0, N_y=6000)
    # signal_eV shape = (rows, cols) = (N_y, N_x) → index as [py, px]
    signal_eV = np.zeros(shape, dtype=np.float32)
    valid = (py >= 0) & (py < shape[0]) & (px >= 0) & (px < shape[1])
    np.add.at(signal_eV, (py[valid], px[valid]), edep[valid])

    mask = signal_eV > 0
    print(f"  Non-zero signal pixels on grid: {mask.sum():,}")

    # ── 3. generate noise matrix ──────────────────────────────────────────────
    # ElectronicNoise: Gaussian(pedestal, sigma) [e-] → eV
    electronic = rng.normal(args.noise_pedestal, args.noise_sigma, shape).astype(np.float32) * args.e2eV
    # DarkCurrent: Poisson(darkcurrent * exp_time) [e-] → eV
    dark = rng.poisson(args.darkcurrent * args.exp_time,
                       size=(shape[1], shape[0])).T.astype(np.float32) * args.e2eV
    noise_eV = electronic + dark

    # ── 4. add noise to signal ────────────────────────────────────────────────
    image_eV = (signal_eV + noise_eV).astype(np.float32)

    # ── 5. optional crop ──────────────────────────────────────────────────────
    if args.crop:
        r0, r1, c0, c1 = [int(x) for x in args.crop.split(",")]
        image_eV  = image_eV [r0:r1, c0:c1]
        signal_eV = signal_eV[r0:r1, c0:c1]
        mask      = mask     [r0:r1, c0:c1]
        print(f"  Cropped to rows {r0}:{r1}, cols {c0}:{c1} → {image_eV.shape}")

    # ── 6. save ───────────────────────────────────────────────────────────────
    npz_path  = os.path.join(args.outdir, f"{stem}_data.npz")
    fits_path = os.path.join(args.outdir, f"{stem}_image.fits")

    np.savez_compressed(npz_path, image=image_eV, signal=signal_eV, mask=mask)
    print(f"Saved {npz_path}")

    _nrows = image_eV.shape[0]
    _ncols = image_eV.shape[1]

    hdr = fits.Header()
    hdr["INFILE"]  = os.path.basename(args.input)
    hdr["BUNIT"]   = "eV"
    hdr["NEVENTS"] = len(tree) if hasattr(tree, '__len__') else -1
    hdr["E2EV"]    = args.e2eV
    hdr["SIGMA"]   = args.noise_sigma
    hdr["DARKCRNT"]= args.darkcurrent
    hdr["EXPTIME"] = args.exp_time
    hdr["SEED"]    = args.seed
    # ── panaSKImg / FitsRawData compatibility headers ─────────────────────────
    # convention keys used in moduletest JSON: NDCM, NCOL, NROW, NPBIN, NSBIN,
    #   AMPL, VCKDIRN, TINTEGR, EXPOSURE, MREAD, DATEINI, DATEEND
    hdr["NDCM"]    = (1,        "Number of skipper samples (1 = not a skipper)")
    hdr["NCOL"]    = (_ncols,   "Number of CCD columns")
    hdr["NROW"]    = (_nrows,   "Number of CCD rows")
    hdr["NPBIN"]   = (1,        "Parallel binning factor")
    hdr["NSBIN"]   = (1,        "Serial binning factor")
    hdr["AMPL"]    = ("SIM",    "Amplifier identifier")
    hdr["VCKDIRN"] = (1,        "Vertical clock direction")
    hdr["TINTEGR"] = (args.exp_time, "Integration time (days)")
    hdr["EXPOSURE"]= (args.exp_time, "Exposure time (days)")
    hdr["MREAD"]   = (0.0,      "Mean readout noise (unused for simulation)")
    hdr["DATEINI"] = ("2024-01-01T00:00:00", "Exposure start (simulation placeholder)")
    hdr["DATEEND"] = ("2024-01-01T00:00:01", "Exposure end   (simulation placeholder)")
    # panaSKImg reads data from extension=1 and headers from ext_header=0.
    # Primary HDU (ext 0): headers only, no data.
    # ImageHDU  (ext 1): image data, no custom headers.
    hdul = fits.HDUList([fits.PrimaryHDU(header=hdr), fits.ImageHDU(data=image_eV)])
    hdul.writeto(fits_path, overwrite=True)
    print(f"Saved {fits_path}")


if __name__ == "__main__":
    main()
