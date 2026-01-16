# def calibrate(file_path, json_name, keep_plot=False):
#     """
#     Take a raw FITS image, perform pedestal subtraction, fit dark current,
#     calibrate, and save the resulting calibrated FITS file.

#     Parameters
#     ----------
#     file_path : str
#         Path to the raw FITS image file.
#     json_name : str
#         Name of the JSON configuration file (must exist in pysimdamicm/json/).
#     keep_plot : bool, optional
#         If True, saves dark current fit plots as PNGs for each amplifier.

#     Returns
#     -------
#     out_fits : str
#         Path to the saved calibrated FITS file.
#     """
#     import os
#     import ROOT
#     import pysimdamicm as ccd
#     from pysimdamicm import io, processes

#     ROOT.gROOT.SetBatch(True)

#     if not os.path.isfile(file_path):
#         raise FileNotFoundError(f"Input FITS not found: {file_path}")

#     cfg_file = os.path.join(ccd.__path__[0], "json", json_name)
#     if not os.path.isfile(cfg_file):
#         raise FileNotFoundError(f"JSON config not found in pysimdamicm/json/: {json_name}")

#     print("[1/8] Loading configuration...")
#     cfg = ccd.utils.config.Config(cfg_file, simulations=False)

#     print("[2/8] Building RawData object...")
#     rdata = io.rawdata.BuilderRawData(file_path, cfg.configuration["input"])

#     print("[3/8] Preparing data...")
#     rdata.prepare_data()

#     # --- Compress skipper data ---
#     print("[4/8] Compressing skipper...")
#     comp = ccd.processes.skipper_analysis.CompressSkipperProcess()
#     comp.func_to_compress = ['mean']
#     comp.id_skip_start = 2
#     comp.execute_process(rdata)

#     # --- Pedestal subtraction ---
#     print("[5/8] Running pedestal subtraction...")
#     ped = ccd.processes.skipper_analysis.PedestalSubtractionProcess()
#     ped.in_overscan = True
#     ped.axis = "row"
#     ped.n_sigma_to_mask = 10
#     ped.use_mad = True
#     ped.method = "gauss_fit"

#     for amp in rdata.amplifier.keys():
#         print(f"    -> Pedestal subtraction for {amp}")
#         rdata.set_amplifier(amp)
#         ped.execute_process(rdata)
#     print("    Pedestal subtraction complete.")

#     # --- Fit dark current ---
#     print("[6/8] Fitting dark current...")
#     dcfit = ccd.processes.skipper_analysis.FitDarkCurrentProcess()
#     rdata.n_run, rdata.n_image = 0, 1

#     for amp in rdata.amplifier.keys():
#         rdata.set_amplifier(amp)
#         dcfit.binning_size = 0.10
#         dcfit.x_min = -1.0
#         dcfit.x_max = 6.0
#         dcfit.sigma_max = 2.0
#         dcfit.n_peaks = 2
#         dcfit.fit_options = "QS"
#         dcfit.do_calibration = True
#         dcfit.calibration = 1.0
#         dcfit.gain_min = 0.8
#         dcfit.gain_max = 1.8

#         dcfit.execute_process(rdata)
#         print(f"    -> Dark current fit complete for {amp}")

#         if keep_plot:
#             c = ROOT.TCanvas(f"c_{amp}", f"Dark Current Fit {amp}", 800, 600)
#             h = dcfit.pcd_hist
#             f = dcfit.fitfunc
#             h.SetTitle(f"Dark Current Fit {amp};Q [e⁻];Counts")
#             h.Draw("HIST")
#             f.SetLineColor(ROOT.kRed)
#             f.Draw("SAME")
#             out_png = f"dark_current_fit_{amp}.png"
#             c.SaveAs(out_png)
#             c.Close()
#             print(f"       Saved fit plot: {out_png}")

#     # --- Calibration ---
#     print("[7/8] Applying calibration...")
#     cal = processes.skipper_analysis.CalibrationProcess()
#     cal.from_dc_fit = True
#     cal.image = "mean_compressed_pedestal_subtracted"
#     cal.execute_process(rdata)
#     print("    Calibration complete.")

#     # --- Save calibrated FITS ---
#     print("[8/8] Saving calibrated FITS...")
#     out_dir = getattr(rdata, "output", os.path.dirname(file_path))
#     out_fits = os.path.join(out_dir, "image_calibrated_e.fits")

#     rdata.SaveAsFits(
#         out_fits,
#         image_names=["image_mean_compressed_pedestal_subtracted"],
#         naming=["CAL_E"]
#     )

#     print(f"Saved calibrated FITS: {out_fits}")
#     return out_fits


def calibrate(file_path, cfg_path, keep_plot=False,
              save=True, return_array=False, out_path=None):
    """
    Take a raw FITS image, perform pedestal subtraction, fit dark current,
    calibrate, and optionally save and/or return the calibrated image.

    Parameters
    ----------
    file_path : str
        Path to the raw FITS image file.
    json_name : str
        Name of the JSON configuration file (must exist in pysimdamicm/json/).
    keep_plot : bool, optional
        If True, saves dark current fit plots as PNGs for each amplifier.
    save : bool, optional
        If True, writes the calibrated image to a FITS file on disk. (default True)
    return_array : bool, optional
        If True, returns the calibrated image as a NumPy array. (default False)
    out_path : str or None, optional
        If provided (and save=True), write the output FITS to this exact path.

    Returns
    -------
    If save and not return_array:
        out_fits : str
    If return_array and not save:
        calibrated_array : np.ndarray
    If both save and return_array:
        (out_fits, calibrated_array) : (str, np.ndarray)
    """
    import os
    import ROOT
    import pysimdamicm as ccd
    from pysimdamicm import io, processes

    # only for optional array return
    import numpy as np
    from astropy.io import fits
    import tempfile

    ROOT.gROOT.SetBatch(True)

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Input FITS not found: {file_path}")

    cfg_file = cfg_path
    if not os.path.isfile(cfg_file):
        raise FileNotFoundError(f"JSON config not found in {cfg_path}")

    print("[1/8] Loading configuration...")
    cfg = ccd.utils.config.Config(cfg_file, simulations=False)

    print("[2/8] Building RawData object...")
    rdata = io.rawdata.BuilderRawData(file_path, cfg.configuration["input"])

    print("[3/8] Preparing data...")
    rdata.prepare_data()

    # --- Compress skipper data ---
    print("[4/8] Compressing skipper...")
    comp = ccd.processes.skipper_analysis.CompressSkipperProcess()
    comp.func_to_compress = ['mean']
    comp.id_skip_start = 2
    comp.execute_process(rdata)

    # --- Pedestal subtraction ---
    print("[5/8] Running pedestal subtraction...")
    ped = ccd.processes.skipper_analysis.PedestalSubtractionProcess()
    ped.in_overscan = True
    ped.axis = "row"
    ped.n_sigma_to_mask = 10
    ped.use_mad = True
    ped.method = "gauss_fit"

    for amp in rdata.amplifier.keys():
        print(f"    -> Pedestal subtraction for {amp}")
        rdata.set_amplifier(amp)
        ped.execute_process(rdata)
    print("    Pedestal subtraction complete.")

    # --- Fit dark current ---
    print("[6/8] Fitting dark current...")
    dcfit = ccd.processes.skipper_analysis.FitDarkCurrentProcess()
    rdata.n_run, rdata.n_image = 0, 1

    for amp in rdata.amplifier.keys():
        rdata.set_amplifier(amp)
        dcfit.binning_size = 0.10
        dcfit.x_min = -1.0
        dcfit.x_max = 6.0
        dcfit.sigma_max = 2.0
        dcfit.n_peaks = 2
        dcfit.fit_options = "QS"
        dcfit.do_calibration = True
        dcfit.calibration = 1.0
        dcfit.gain_min = 0.8
        dcfit.gain_max = 1.8

        dcfit.execute_process(rdata)
        print(f"    -> Dark current fit complete for {amp}")

        if keep_plot:
            c = ROOT.TCanvas(f"c_{amp}", f"Dark Current Fit {amp}", 800, 600)
            h = dcfit.pcd_hist
            f = dcfit.fitfunc
            h.SetTitle(f"Dark Current Fit {amp};Q [e⁻];Counts")
            h.Draw("HIST")
            f.SetLineColor(ROOT.kRed)
            f.Draw("SAME")
            out_png = f"dark_current_fit_{amp}.png"
            c.SaveAs(out_png)
            c.Close()
            print(f"       Saved fit plot: {out_png}")

    # --- Calibration ---
    print("[7/8] Applying calibration...")
    cal = processes.skipper_analysis.CalibrationProcess()
    cal.from_dc_fit = True
    cal.image = "mean_compressed_pedestal_subtracted"
    cal.execute_process(rdata)
    print("    Calibration complete.")

    # --- Save and/or return array ---
    print("[8/8] Finalizing output...")

    image_key = "image_mean_compressed_pedestal_subtracted"

    out_fits = None
    if save:
        if out_path:
            out_dir = os.path.dirname(out_path) or os.path.dirname(file_path)
            os.makedirs(out_dir, exist_ok=True)
            out_fits = out_path
        else:
            out_dir = getattr(rdata, "output", os.path.dirname(file_path)) or os.path.dirname(file_path)
            os.makedirs(out_dir, exist_ok=True)
            out_fits = os.path.join(out_dir, "image_calibrated_e.fits")

        rdata.SaveAsFits(
            out_fits,
            image_names=[image_key],
            naming=["CAL_E"]
        )
        print(f"Saved calibrated FITS: {out_fits}")

    calibrated_array = None
    if return_array:
        read_path = out_fits
        tmp_path = None

        if not save:
            tmp_dir = getattr(rdata, "output", os.path.dirname(file_path)) or os.path.dirname(file_path)
            os.makedirs(tmp_dir, exist_ok=True)
            fd, tmp_path = tempfile.mkstemp(prefix="calib_tmp_", suffix=".fits", dir=tmp_dir)
            os.close(fd)
            rdata.SaveAsFits(
                tmp_path,
                image_names=[image_key],
                naming=["CAL_E"]
            )
            read_path = tmp_path

        with fits.open(read_path) as hdul:
            calibrated_array = np.array(hdul[0].data, copy=True)

        if tmp_path is not None:
            try:
                os.remove(tmp_path)
            except OSError:
                pass

    # Return according to selected behavior
    if save and not return_array:
        return out_fits
    if return_array and not save:
        return calibrated_array
    if save and return_array:
        return out_fits, calibrated_array





