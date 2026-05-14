from astropy.io import fits

# -------------------------------
# INPUT / OUTPUT
# -------------------------------
infile  = "/home/dsalma/damicm-ml/science_run_data/compose_103_ccdA_final.fits"
outfile = "/home/dsalma/damicm-ml/science_run_data/cut_img/compose_103_ccdA_cut_6300x8000.fits"

# -------------------------------
# CROP REGION (x = cols, y = rows)
# -------------------------------
x_min = 0
x_max = 6300
y_min = 0
y_max = 8000


n_cols = x_max - x_min
n_rows = y_max - y_min

# -------------------------------
# CROP
# -------------------------------
with fits.open(infile) as hdul:
    new_hdus = []

    for i, hdu in enumerate(hdul):
        header = hdu.header.copy()
        data   = hdu.data

        if data is not None:
            # FITS images are indexed [row, col] = [y, x]
            data_cut = data[y_min:y_max, x_min:x_max]

            # Update geometry keywords if present
            for key, val in [("NCOL", n_cols), ("NROW", n_rows)]:
                if key in header:
                    header[key] = val

            # Shift WCS reference pixel if present
            if "CRPIX1" in header:
                header["CRPIX1"] -= x_min
            if "CRPIX2" in header:
                header["CRPIX2"] -= y_min

            if i == 0:
                new_hdus.append(fits.PrimaryHDU(data=data_cut, header=header))
            else:
                new_hdus.append(
                    fits.CompImageHDU(
                        data=data_cut,
                        header=header,
                        compression_type="RICE_1"
                    )
                )
        else:
            # No image data → keep HDU structure
            if i == 0:
                # Primary header still needs updated geometry
                for key, val in [("NCOL", n_cols), ("NROW", n_rows)]:
                    if key in header:
                        header[key] = val
                new_hdus.append(fits.PrimaryHDU(header=header))
            else:
                new_hdus.append(hdu.copy())

    fits.HDUList(new_hdus).writeto(outfile, overwrite=True)
