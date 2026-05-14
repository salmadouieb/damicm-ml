#--------------------------------
# Clustering with GMM
# Full analysis script with root file output
# Imports
import time
from astropy.io import fits
from sklearn.mixture import GaussianMixture
import pandas as pd
import numpy as np
from scipy.ndimage import label as cc_label
from sklearn.preprocessing import StandardScaler
import uproot
import awkward as ak

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


# ---------------------------------------------------------------------------
# Cluster class — mirrors data_formats.py Cluster exactly frm pysimdamicm
# ---------------------------------------------------------------------------

class Cluster:

    def __init__(self, cluster_id, pixels_x, pixels_y, pixels_E):
        self.cluster_id = int(cluster_id)

        # 1-indexed, matching: self.pixels_x = x[mask].astype(int) + 1
        self.pixels_x = pixels_x.astype(int) + 1
        self.pixels_y = pixels_y.astype(int) + 1
        self.pixels_E = pixels_E.astype(float)

        # STD_XY computed in __init__ before get_properties_on_axes
        self.STD_XY = (
            ((self.pixels_x - self.pixels_x.mean())**2 +
             (self.pixels_y - self.pixels_y.mean())**2).sum()
            / len(self.pixels_x)
        ) ** 0.5

        self.Energy = self.pixels_E.sum()
        self.Npix   = int(len(self.pixels_x))

    def get_cluster_properties(self, is_simulation=False, __DEBUG__=False, get_fitted_STD=False):
        self.is_simulation = is_simulation
        self.get_properties_on_axes('x')
        self.get_properties_on_axes('y')
        ###
        setattr(self,"wSTD_XY",np.sqrt((getattr(self,"wSTD_X")**2 + getattr(self,"wSTD_Y")**2)/2.))
        if get_fitted_STD:
            self.get_fitted_STD(__DEBUG__)
        self.get_charge_properties()

    def get_properties_on_axes(self, axis, stats_list=['mean', 'min', 'max']):
        values = getattr(self, 'pixels_{}'.format(axis))

        for stat_name in stats_list:
            statistic_isnt = getattr(np, stat_name)
            setattr(self, stat_name + axis.upper(), float(statistic_isnt(values)))

        setattr(self, 'RMS' + axis.upper(), float(np.mean(values**2.0)))
        setattr(self, 'D'   + axis.upper(), float(values.max() - values.min()) + 1)
        setattr(self, 'Pos' + axis.upper(), float(np.average(values, weights=self.pixels_E)))
        setattr(self, 'STD_' + axis.upper(), float(values.std()))

        w = self.pixels_E / self.pixels_E.sum()
        w = w / w.sum()
        wmean = (values * w).sum()
        setattr(self, 'dwSTD_' + axis.upper(),
                float(((w * (values - wmean)**2.0).sum() / float(len(w)))**0.5))

        wvalues = getattr(self, 'Pos' + axis.upper())
        setattr(self, 'wSTD_' + axis.upper(),
                float(np.sqrt((self.pixels_E * (values - wvalues)**2).sum()
                              / self.pixels_E.sum())))

    def get_charge_properties(self):
        ind = np.where(self.pixels_E == self.pixels_E.max())[0][0]
        self.Qmax  = float(self.pixels_E[ind])
        self.QmaxX = float(self.pixels_x[ind])
        self.QmaxY = float(self.pixels_y[ind])

    def get_fitted_STD(self, __DEBUG__=False):
        import ROOT
        _isbatch = ROOT.gROOT.IsBatch()
        ROOT.gROOT.SetBatch(1)

        px = ROOT.TH1D("px","px", int(self.pixels_x.max()-self.pixels_x.min())+3, self.pixels_x.min()-1.5, self.pixels_x.max()+1.5)
        py = ROOT.TH1D("py","py", int(self.pixels_y.max()-self.pixels_y.min())+3, self.pixels_y.min()-1.5, self.pixels_y.max()+1.5)

        for x,y,e in zip(self.pixels_x,self.pixels_y,self.pixels_E):
            _ = px.Fill(x,e/self.pixels_E.max())
            _ = py.Fill(y,e/self.pixels_E.max())

        setattr(self,"projectionX", np.zeros(px.GetNbinsX()).astype(float))
        setattr(self,"projectionX_cols", np.linspace(self.pixels_x.min()-1.5, self.pixels_x.max()+1.5,int(self.pixels_x.max()-self.pixels_x.min())+3).astype(float))
        for i in range(px.GetNbinsX()):
            self.projectionX[i] = px.GetBinContent(i+1)
        setattr(self,"projectionY", np.zeros(py.GetNbinsX()).astype(float))
        setattr(self,"projectionY_rows", np.linspace(self.pixels_y.min()-1.5, self.pixels_y.max()+1.5,int(self.pixels_y.max()-self.pixels_y.min())+3).astype(float))
        for i in range(py.GetNbinsX()):
            self.projectionY[i] = py.GetBinContent(i+1)

        fit_opt = "Q L"
        # fitting x projection to gauss
        gausx = ROOT.TF1("gaus_x","gaus")
        gausx.SetLineColor(2)
        gausx.SetParameter(1,self.PosX)
        gausx.SetParLimits(1,self.minX,self.maxX)
        px.SetTitle(";pixels_x;counts")
        px.Fit(gausx,fit_opt)
        setattr(self,"fSTD_X", float(gausx.GetParameter(2)))
        setattr(self,"fPosX",  float(gausx.GetParameter(1)))

        # fitting y projection to gauss
        gausy = ROOT.TF1("gaus_y","gaus")
        gausy.SetLineColor(2)
        gausy.SetParameter(1,self.PosY)
        gausy.SetParLimits(1,self.minY,self.maxY)
        py.SetTitle(";pixels_y;counts")
        py.Fit(gausy,fit_opt)
        setattr(self,"fSTD_Y",float(gausy.GetParameter(2)))
        setattr(self,"fPosY", float(gausy.GetParameter(1)))

        ROOT.gROOT.SetBatch(_isbatch)


# ---------------------------------------------------------------------------
# Load FITS image and run GMM
# ---------------------------------------------------------------------------

t_start = time.time()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Path to input FITS file")
args = parser.parse_args()
fits_path = args.input

with fits.open(fits_path) as hdul:
    cal_img = hdul[1].data
print("Image successfully loaded from:", fits_path)

# convert from eV (simulation output) to keV
cal_img = cal_img.T / 1000

print("Begin clustering process...")

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
        "points":  list(blobs.values())
    }
)

df_blobs["cluster_size"] = df_blobs["points"].apply(len)
df_blobs = df_blobs[df_blobs["cluster_size"] > 1].reset_index(drop=True)

print("Clustering complete. Found {} clusters!".format(len(df_blobs)))


# ---------------------------------------------------------------------------
# Instantiate Cluster objects and compute features
# ---------------------------------------------------------------------------

print("Computing features...")

cluster_list = []
for _, row in df_blobs.iterrows():
    pts = row["points"]
    cls = Cluster(
        cluster_id = row["blob_id"],
        pixels_x   = pts[:, 1].astype(float),
        pixels_y   = pts[:, 0].astype(float),
        pixels_E   = pts[:, 2].astype(float),
    )
    cls.get_cluster_properties(get_fitted_STD=True)
    cluster_list.append(cls)

# Compute closest cluster -- same logic from reconstruction.py get_closest_cluster.
# Note: original only runs this for clusters with has_seed. We have no has_seed
# so we treat all clusters as seeded (mask = all True).
PosX = np.array([c.PosX for c in cluster_list])
PosY = np.array([c.PosY for c in cluster_list])

if len(cluster_list) <= 1:
    for cls in cluster_list:
        cls.closest_cls_dist   = float(-1.0)
        cls.closest_cls_dist_X = float(-1.0)
        cls.closest_cls_dist_Y = float(-1.0)
        cls.closest_cls_id     = float(-1.0)
else:
    for idx, cls in enumerate(cluster_list):
        distances = np.sqrt((PosX - cls.PosX)**2.0 + (PosY - cls.PosY)**2.0)
        # exclude self using same logic as original
        mask = np.logical_and(PosX != cls.PosX, PosY != cls.PosY)
        distances = np.where(mask, distances, np.inf)
        cls.closest_cls_dist   = float(distances[np.argmin(distances)])
        cls.closest_cls_dist_X = float(abs(PosX[np.argmin(distances)] - cls.PosX))
        cls.closest_cls_dist_Y = float(abs(PosY[np.argmin(distances)] - cls.PosY))
        cls.closest_cls_id     = float(np.argmin(distances) + 1)  # 1-based index, matching original

print("Features computed.")

# ---------------------------------------------------------------------------
# Assemble results into dataframe
# ---------------------------------------------------------------------------

SCALAR_ATTRS = [
    "cluster_id", "Npix", "Energy",
    "DX", "DY", "meanX", "meanY", "minX", "minY", "maxX", "maxY",
    "RMSX", "RMSY", "STD_X", "STD_Y", "STD_XY",
    "PosX", "PosY",
    "wSTD_X", "wSTD_Y", "wSTD_XY", "dwSTD_X", "dwSTD_Y",
    "Qmax", "QmaxX", "QmaxY",
    "closest_cls_id", "closest_cls_dist", "closest_cls_dist_X", "closest_cls_dist_Y",
    "fPosX", "fPosY", "fSTD_X", "fSTD_Y",
]

for attr in SCALAR_ATTRS:
    df_blobs[attr] = [getattr(c, attr) for c in cluster_list]

df_blobs["pixels_x"]         = [c.pixels_x         for c in cluster_list]
df_blobs["pixels_y"]         = [c.pixels_y         for c in cluster_list]
df_blobs["pixels_E"]         = [c.pixels_E         for c in cluster_list]
df_blobs["projectionX"]      = [c.projectionX      for c in cluster_list]
df_blobs["projectionX_cols"] = [c.projectionX_cols for c in cluster_list]
df_blobs["projectionY"]      = [c.projectionY      for c in cluster_list]
df_blobs["projectionY_rows"] = [c.projectionY_rows for c in cluster_list]

# ---------------------------------------------------------------------------
# Write ROOT file with uproot
# ---------------------------------------------------------------------------

SCALAR_FLOAT_COLS = [
    "DX", "DY", "Energy",
    "PosX", "PosY", "Qmax", "QmaxX", "QmaxY",
    "RMSX", "RMSY",
    "STD_X", "STD_XY", "STD_Y",
    "closest_cls_dist", "closest_cls_dist_X", "closest_cls_dist_Y",
    "closest_cls_id",
    "dwSTD_X", "dwSTD_Y",
    "fPosX", "fPosY", "fSTD_X", "fSTD_Y",
    "meanX", "meanY", "minX", "minY", "maxX", "maxY",
    "wSTD_X", "wSTD_Y", "wSTD_XY",
]

SCALAR_INT_COLS = [
    "blob_id", "cluster_id", "cluster_size", "Npix",
]

JAGGED_FLOAT_COLS = [
    "pixels_E", "pixels_x", "pixels_y",
    "projectionX", "projectionX_cols", "projectionY", "projectionY_rows",
]

out_path = "gmm_clusters.root"

branch_data = {}
for col in SCALAR_FLOAT_COLS:
    branch_data[col] = df_blobs[col].to_numpy(dtype=np.float32)
for col in SCALAR_INT_COLS:
    branch_data[col] = df_blobs[col].to_numpy(dtype=np.int32)
for col in JAGGED_FLOAT_COLS:
    branch_data[col] = ak.Array(df_blobs[col].tolist())

with uproot.recreate(out_path) as root_file:
    root_file["clustersRec"] = branch_data

print(f"Wrote {len(df_blobs)} clusters to {out_path} (tree: clustersRec). Bosh!")
print(f"Done!! :) Total time: {time.time() - t_start:.1f}s")
