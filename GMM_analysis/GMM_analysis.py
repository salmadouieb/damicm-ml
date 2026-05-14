#--------------------------------
# Clustering with GMM
# Full analysis script with root file output
# Imports
from astropy.io import fits
from sklearn.mixture import GaussianMixture
import pandas as pd
import numpy as np
from scipy.ndimage import label as cc_label
from sklearn.preprocessing import StandardScaler

# Helpers
def blobs_from_gmm(X_raw, gmm_labels, frame_size, cluster_id):
    """
    X_raw      : (N,3) array [y, x, E]
    gmm_labels : (N,) GMM labels (same order as X_raw)
    frame_size : (H, W)
    cluster_id : GMM component to extract

    Returns:
      blobs   : dict {blob_id: array [[y,x,E], ...]}
      labeled : (H,W) blob-ID image
    """
    H, W = frame_size

    sel = (gmm_labels == cluster_id)
    if not np.any(sel):
        labeled = np.zeros((H, W), dtype=int)
        return {}, labeled

    y = X_raw[sel, 0].astype(int)
    x = X_raw[sel, 1].astype(int)
    E = X_raw[sel, 2].astype(float)

    img = np.zeros((H, W), dtype=np.uint8)
    img[y, x] = 1

    # HARD-CODED 8-connectivity
    labeled = cc_label(img, structure=np.ones((3, 3), dtype=int))[0]

    blob_ids = labeled[y, x]
    blobs = {}

    for bid in np.unique(blob_ids):
        if bid == 0:
            continue
        m = (blob_ids == bid)
        blobs[bid] = np.column_stack([y[m], x[m], E[m]])  # [y,x,E]

    return blobs, labeled

def weighted_std_xy(points):
    """
    points: array of shape (N,3) with columns [y, x, E]
    Returns: wSTD_y, wSTD_x
    """
    y = points[:, 0]
    x = points[:, 1]
    E = points[:, 2]

    Etot = E.sum()
    if Etot == 0:
        return np.nan, np.nan

    y_bar = np.sum(E * y) / Etot
    x_bar = np.sum(E * x) / Etot

    wSTD_y = np.sqrt(np.sum(E * (y - y_bar)**2) / Etot)
    wSTD_x = np.sqrt(np.sum(E * (x - x_bar)**2) / Etot)

    return wSTD_y, wSTD_x

def make_samples(masked_img):
    """
    Build N x 3 samples from a thresholded CCD image (in keV).
    Uses raw E (no log) as the third feature.

    Returns:
      X_raw : (N,3) with columns [y, x, E]
      idx_yx: (N,2) integer pixel coords [y, x]
    """
    valid = np.isfinite(masked_img) & (masked_img > 0)
    ys, xs = np.nonzero(valid)
    E = masked_img[ys, xs].astype(np.float64)

    X_raw = np.column_stack([ys, xs, E])
    idx_yx = np.column_stack([ys, xs])
    return X_raw, idx_yx

REG_COVAR = 1e-6
RANDOM_STATE = 209

def fit_gmm(X_raw, reg_covar=REG_COVAR, random_state=RANDOM_STATE):
    """
    Standardize features and fit a 2-component full-covariance GMM.
    Returns labels in the original sample order.
    """
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X_raw)

    gmm = GaussianMixture(
        n_components=2,
        covariance_type="full",
        reg_covar=reg_covar,
        random_state=random_state
    )
    gmm.fit(X_std)
    return gmm.predict(X_std)


# input file, this will eventually be a flag
fits_path = "/home/dsalma/damicm-ml/simulation/test_2/out_simulCCDimg_logb9c4085_CCDSensor_PV_55a26z1_5000_s0_image.fits"
with fits.open(fits_path) as hdul:
    cal_img = hdul[1].data

# convert from eV (simulation output) to keV
cal_img = cal_img.T / 1000

X_raw_cut, idx_yx = make_samples(cal_img)
labels = fit_gmm(X_raw_cut)

blobs, labeled = blobs_from_gmm(
    X_raw=X_raw_cut,
    gmm_labels=labels,
    frame_size=cal_img.shape,
    cluster_id=1
)

df_blobs = pd.DataFrame(
    {
        "blob_id": list(blobs.keys()),
        "points": list(blobs.values())
    }
)

df_blobs["cluster_size"] = df_blobs["points"].apply(len)

df_blobs["energy"] = df_blobs["points"].apply(
    lambda pts: np.sum(np.asarray(pts)[:, 2])
)

df_blobs[["wSTD_y", "wSTD_x"]] = (
    df_blobs["points"]
    .apply(lambda pts: pd.Series(weighted_std_xy(pts)))
)

df_blobs["wSTD_xy"] = np.sqrt(
    0.5 * (df_blobs["wSTD_x"]**2 + df_blobs["wSTD_y"]**2)
)