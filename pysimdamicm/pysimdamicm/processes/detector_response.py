""":mod:`pysimdamicm.detector_response`

A package to handle different digitilize process involved in the response 
of the DAMIC-M detector (process involved either in sensor or electronics 
response).

Note
----
Geant4 system units as default, see :obj:`pysimdamicm.io.G4utils.Units` 

.. moduleauthor:: Nuria Castello-Mor
"""
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.stats import binned_statistic_2d
import pandas as pd
from glob import glob
from astropy.io import fits
import warnings
#warnings.simplefilter('ignore', category=fits.verify.VerifyWarning)
import os
import random

from pysimdamicm import __path__ as __module_path__
from pysimdamicm.utils.units import Units
u = Units()

##### module DEBUG VARIABLES
__verbose__ = False
__debugs__ = False

from time import time
import random

###############################################################################################
#####       ABSTRACT CLASS TO DEFINE THE DIFFERENT PROCESS 
#####           PART OF THE DETECTOR RESPONSE
###############################################################################################
class DigitizeProcess(metaclass=ABCMeta):
    """Abstract Class for the digitilize process.
    
    Each process involved in the response of the DAMIC-M detector (based on CCDs) will 
    be created from this abstract class.

    Example
    -------
    
    1. Diffusion process
        >>> import pysimdamicm as dam
        >>> diff = dam.detector_response.Diffusion()
        >>> diffusion.__name__
        'Diffusion'
        >>> diffusion.__sequence_id__
        10
        >>> diffusion.__units__
        {'z_offset': 0.001, 'B': 1000.0, 'A': 1e-06, 'alpha': 1.0, 'fanofactor': 1}
        >>> diff_params = {'A':159.0, 'B': 0.0001, 'fanofactor':0.10 }
        >>> diff.set_parameters(**diff_params)
        >>>  

    2. DarkCurrent process
        >>> import pysimdamicm as dam
        >>> dc = dam.detector_response.DarkCurrent()
        >>> dc.__name__
        'DarkCurrent'
        >>> dc.__sequence_id__
        40
        >>> dc.__units__
        {'img_size': 1, 'mode': 1, 'darkcurrent': 1.1574074074074073e-05, 'exp_time': 86400.0,'min_enlarge': 1}
        >>> dc_params = {'exp_time':8/24., 'darkcurrent': 0.001 , 'mode':''}
        >>> dc.set_parameters(**dc_params)
        Running on mode: full-image
    >>> dc.dark_current_image
    array([[0., 0., 0., ..., 0., 0., 0.],
           [0., 0., 0., ..., 0., 0., 0.],
           [0., 0., 0., ..., 0., 0., 0.],
           ...,
           [0., 0., 0., ..., 0., 0., 0.],
           [0., 0., 0., ..., 0., 0., 0.],
           [0., 0., 0., ..., 0., 0., 0.]])
    >>> dc.dark_current_image.mean()
       0.0012602795833333822
 
    """
    __seed__ = 321
    __units__ = {}

    def __init__(self,seed=None):
        if seed is not None:
            DigitizeProcess.__seed__ = int(seed)
        
        ## Set the seed for the random generation process
        ## np.random.seed(DigitizeProcess.__seed)
        self.__verbose__ = False
        self.__display__ = False
        self.__DEBUG__   = False 
        self.save_image  = False
        self.save_plots  = False

        self.__units__ = {
                "__display__":1,
                "__DEBUG__"  :1,
                "__verbose__":1,
                "save_image" :1,
                "save_plots" :1,
                }


    #def __str__(self):
    #    """String shown when you use the print() function
    #    """
    #    msg = '<Digitize Process: {}>\n'.format(self.__class__.__name__)
    #    # FIXME Problem with attr.find('_') != 0 --> RecursionError WTF??
    #    for attribute in filter(lambda attr: attr. find('_') == -1,dir(self)):
    #        try:
    #            msg += ' {}: {}\n'.format(attribute,getattr(self,attribute))
    #        except AttributeError:
    #            pass
    #    return msg[:-1]

    #def __repr__(self):
    #    """String shown when you push return over an instance of this class
    #    """
    #    return self.__str__()

    def set_parameters(self, **model_parameters):
        """Runtime method to set all atributes of the digitize process listed on the configuration
        json file.

        Those process's parameters not listed on the configuration json file will set to the default
        values.

        Parameters
        ----------
        **model_parameters : dict
            Keyword arguments for all attributes realted with the digitizer process (see examples).

        """
        for keyword,val in model_parameters.items():
            if keyword in ["units","unit"]:
                continue
            if hasattr(self,keyword):
                if type(val)==dict:
                    # some parameters are more complicates than single values, or list, and units can
                    # not be added XXX Case of the Signal Patter Recognition Process XXX
                    setattr(self,keyword,val)
                else:
                    setattr(self,keyword,val*self.__units__[keyword])
            else:
                raise(AttributeError('"{}" invalid attribute for class {}'.format(keyword,self.__class__.__name__)))

    @property
    def mode(self):
        """Property **mode** to define how the intrinsic detector should be simulated. 
        It is mainly related to the size of the noise image, and if this image should be created
        for each pair of values (event,ccd) or not. There are two modes:
            
            **full-image**  
                
                at the beginning of the process one image noise (DarkCurrent and/or ElectronicNoise) will be
                created and used for all the pairs (event,ccd). The shape of this image corresponds to the 
                shape of the sensitive detector (defined at :obj:`Units.ccd_shape`)
    
            **cropped-image**
                
                the default mode. For each event an image will be re-created with a shape smaller than the 
                sensitive detector. This mode is faster than the former, the filtering before ClusterFinder has less 
                pixels to evaluate.
        """
        return self.__mode
    @mode.setter
    def mode(self,val):
        """Setter for mode property 
        """
        self.recreate_image = False
        self.full_image_counter = 0
        
        if not self.__name__ in ['DarkCurrent', 'ElectronicNoise']:
            msm="The attribure 'mode' is only implemented for process related with the intrinsice detector noise"
            raise AttributeError(msm)

        if val == 'cropped-image':
            print('     Set mode to `cropped-image`.')
            self.execute_process = self.__execute_process_cropped_image__
        elif val == 'full-image':
            self.execute_process = self.__execute_process_full_image__
            print('     Set mode to `full-image`.')
        else:
            msm = "Mode not implemented. Possible modes are: \n"+\
                    "- full-image: process that will be executed once during reconstruction,"+\
                    " so all events will share its output (i.e. noise)\n"+\
                    "- cropped-image: process will be executed for each hit collection"
            raise NotImplementedError(msm)
        self.__mode = val

    @abstractmethod
    def execute_process(self):
        """Metaclass abstract method to be implemented for each specific
        process.
        """
        raise NotImplementedError()

    def __del__(self,hits=None):
        """
        """
        pass

###############################################################################################
#####       ELECTRON DIFFUSION MODEL
###############################################################################################
class Diffusion(DigitizeProcess):
    """Class to generate the creation of electron-hole pairs and its diffusion, 
    for silicon detector.
    
    **Charge e-h pair Creation**

    A delta-function energy lost at :math:`p(x,y,z)`, within the silicon bulk of the CCD, will
    generate :math:`N_{Q}` e-h pairs, which is given by
    
    
    .. math::
       
       N_{Q} = integer \left[ \mathcal{G}(\mu,\sigma) \\right]
    

    where 
    
    .. math::
       :nowrap:

       \\begin{equation}
         \mu = \delta E/ \epsilon_{m}
       \end{equation}

        
    .. math::
       :nowrap:

       \\begin{equation}
           \sigma = \sqrt{(F * \mu)}
       \end{equation}

    where :math:`F` is the fano factor and :math:`\epsilon_{m}` is the average 
    energy loss required to produce a charge pair (e-h pair).
    
    **Charge diffusion**
   
    The transverse diffusion of delta-function charge deposition in
    the bulk of the CCD is described by a Gaussian with
    :math:`\sigma_{xy}=\sqrt{2 D_{p} t}` where :math:`D{p}` is the
    diffusion coefficient of the holes and :math:`t` is the time
    elapssed since creation. By using the linear relationship between drift velocity of charge 
    carriers and electric field in the semiconductor, :math:`\sigma_{xy}` 
    measured by the detector is
    just a function of the starting z-position in the bulk

    .. math::

       \sigma_{xy}^{2}(z,E=0)  = - A \ln{( 1-B(z- z_{offset}) )}

    Where :math:`A` and :math:`B` are coefficients based on the
    physical parameters of the silicon detector (permitivity of
    silicon, donor charge density, operating temperature, bias across
    the CCD, and the thickness of the active region).


    It is obvious to think that the transversal diffusion has an energy dependence: a higher 
    energy deposition, stronger the created local electric field. The nal diusion model is therefore

    .. math::
       :nowrap:

       \\begin{equation}
         \sigma_{xy}(z,E) = \sigma_{xy}(z,E=0) \left( \\alpha_{0} + \\frac{\\alpha E}{ \sigma_{xy}(z=z_{offset}, E=0)} \\right)
       \end{equation}
    
    where :math:`\sigma_{xy}(z,E=0)` is shown in the previous equation.
    
    **See more details** at `DAMIC-100 notes <https://www.overleaf.com/3446968864nchrgrjbmccq>`_
    
    
    **List of Atributes:**

    Attributes
    ----------
        
        A 
            amplitude of the :math:`\sigma_{xy}-z` relationship (depends on the silicon bulk parameters: silicon permeability, donor charge density, ...)
            
        B
            shape of the :math:`\sigma_{xy}-z` relationship
        
        fanofactor
            
            value related with the conversion of absorbed energy into signal quanta, related with the fluctuation of the signal for a given absorbed energy

        z_offset

            thickness of the active region (fully depleted region)

    only_pair_creation
        
        variable to activate/deactive diffusion. 0: diffusion on, 1: diffusion off

.. warning:: 

    * y-axis as the columns (horitzontal register), 
    * x-axis as the rows (serial register), 
    * z-axis as the direction of the electric field in the dopped region of the silicon detector.
    * Geant4 units: mm, eV, second (see :obj:`pysimdamicm.io.G4utils.Units`)
    * 3.77eV as the average energy lost required to produce a charge pair in the bulk of the CCD, which can be changed from :obj:`pysimdamicm.io.G4utils.Units.e2eV`

    """
    __sequence_id__ = 10
    __name__ = 'Diffusion'

    __E_in_eV__ = True

    def __init__(self):
        super().__init__()

        self.model = "gauss"
        self.only_pair_creation = 0


        #### Parameters of the model
        self.A = (float(216.2)*u.um*u.um)
        self.B = (float(0.000886)*1/u.um)
        self.z_offset = float(0.0)*u.um
        self.fanofactor = float(0.129)
        self.alpha = float(0.00597*1e-3)*u.pixel/u.eVee
        self.alpha0 = 1.0 
        self.z_max = float(669.0)*u.um
        self.apply_pcc = 0
        self.pcc_debug =0
        self.repulsion = 0                      # 0 = off, 1 = on
        self.sigma_map_file = None              # path to sigma map file (.xlsx or .csv)
        self.sigma_interpolator = None          # store interpolation object
        self.sigma_map_loaded = False           # flag when map successfully loaded
        
        
        ###PCC function
        from scipy.interpolate import interp1d
        thick= u.ccd_thickness/u.um   # CCD thickness from json file in microns
        #print(thick)
        #z_um_data = np.array([
        #669.0, 668.5, 668.0, 667.5,        # dead layer
        #667.0, 666.5, 666.0, 665.5, 665.0, 664.5, 664.0,
        #663.5, 663.0, 662.5, 662.0, 661.5
        #])
        z_um_data=np.array([thick+0.0,thick-0.5,thick-1.0,thick-1.5,thick-2.0,thick-2.5,
            thick-3.0,thick-3.5,thick-4.0,thick-4.5,thick-5.0,thick-5.5,thick-6.0,
            thick-6.5, thick-7.0, thick-7.5])
        #print(z_um_data)
        eps_vals = np.array([ # dead layer
        0.00781727, 0.027572, 0.0789126, 0.157685, 0.247863, 0.34214,
        0.437597, 0.533334, 0.629258, 0.725383, 0.821692, 0.918136,
        0.991682, 0.999804, 0.999997, 1.0
        ])
        self.eps_func = interp1d(z_um_data, eps_vals, kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
        ###PCC FUNCTION ENDS

        self.__units__={'A':u.um**2.0,'B':1/u.um,'z_offset':u.um,'fanofactor':1,'alpha':u.pixel/u.eVee,'model':1,
                'z_max':u.um,'alpha0':1, 'only_pair_creation':1,'apply_pcc': 1, 'pcc_debug':1, 'repulsion': 1,
                'sigma_map_file': 1}


        datafile = "{}/data/p100K.npz".format(__module_path__[0])
        ehdata = np.load(datafile)['data']
        self.ehenergies = ehdata[:,0].flatten()
        self.ehprobs = ehdata[:,1::]


    from scipy.interpolate import RegularGridInterpolator
    ##Loading sigma_xy map if repuslion is on in json ##
    
    def load_sigma_map_lazy(self):
        from scipy.interpolate import RegularGridInterpolator

        """Load sigma map on first use (called inside execute_process)."""
        if getattr(self, 'repulsion', 0) != 1:
            return

        if self.sigma_map_loaded:
            return

        try:
            datafile_1 = self.sigma_map_file
            if not os.path.isabs(datafile_1):
                datafile_1 = os.path.join(__module_path__[0], 'data', datafile_1)

            if not os.path.exists(datafile_1):
                print(f"[WARNING] Sigma map file not found: {datafile_1}")
                self.repulsion = 0
                return

            if pd is None:
                raise RuntimeError("pandas required to read Excel/CSV sigma_map_file")

            print(f"[Diffusion] Lazy-loading sigma map from {datafile_1}")
            df = pd.read_excel(datafile_1, index_col=0)
            sigma_map = df.values.astype(float)
            energy_vals = df.index.to_numpy(dtype=float)
            depth_vals  = df.columns.to_numpy(dtype=float)

            self.sigma_interpolator = RegularGridInterpolator(
                (energy_vals, depth_vals),
                sigma_map,
                bounds_error=False, fill_value=None
            )
            self.sigma_map_loaded = True

            print(f"[Diffusion] Loaded σ_map: shape={sigma_map.shape}, ")
                 # f"E-range={energy_vals.min():.3f}–{energy_vals.max():.3f} keV, "
                  #f"z-range={depth_vals.min():.6f}–{depth_vals.max():.6f} m")

        except Exception as e:
            print(f"[Diffusion]  Could not load sigma map: {e}")
            self.sigma_map_loaded = False
            self.repulsion = 0
        
    def execute_process(self,hits):

        if self.model == "gauss":
            self.execute_process_gauss(hits)
        else:
            self.execute_process_karthik(hits)

    def execute_process_gauss(self,hits):
        """:obj:`Diffusion.execute_process` is a method to model the charge pair creation 
        and charge diffusion for each energy loss given by the numpy.ndarray `hits.Edep`.

        While the parameters of the model :math:`A`, :math:`B` and :math:`z_{offset}` were obtained by
        assuming units of **um**, the positions for hits are in **mm** (default units of geant4).
        **Internally the positions will be converted to um, but stored in mm.**

        Parameters
        ----------
        hits : :obj:`pysimdamicm.utils.data_foramts.G4HitCollection`
            A collection of hits (subset of hits from the geant4
            simulations) with at least the following attributes: time, and
            `__column_axis__`.
    
        Returns
        -------
        hits : :obj:`pysimdamicm.utils.data_foramts.G4HitCollection`
            with actualized values of x,y,z, Edep, time and pdg

            Units (as the inputs ones):
                * Edep in eV
                * position in mm
                * time in seconds
                * pdg adimensional

        Note
        ----
        **hits** must have **x**, **y**, **z** and **Edep** as attributes, if not an
        **AttributeError** will be raised.

        """
        self.load_sigma_map_lazy()
        ### COLLECTION OF IONIZED ELECTRONS PER ENERGY DEPOSITION
        ###########################################################################################
        ### Number of electron-hole pairs given an energy loss of Edep eV is given by a gaussian
        ###   distribution with mu and sigma as follows:
        mu = hits.Edep / u.e2eV
        sigma = (self.fanofactor * mu)**0.5
        prob_n_obs_e = np.random.normal(mu,sigma)
        ### For one particle (a set of steps), casting to int the conversion is always getting less 
        ##  pairs than the one should we have for the sum of all Edeps (per particle)
        ##  Using round, we are now assuming an step function where 50 per cent of the times we are
        ##  using an eh energy conversion lower than 3.77 and the rest, an energy conversion higher 
        ##  Which gives us an average of 3.77 which is what we want.
        n_obs_e = np.round( prob_n_obs_e ).astype(int)

        if self.apply_pcc == 1:
             z_um = hits.z *(u.mm/u.um)  # Convert mm to microns
             #print("[PCC] Applied to hits at z (µm):", z_um[:10])
             eps = self.eps_func(z_um)
             # Identify hits in the PCC transition region
             low_z_mask = (z_um >= 663.0) & (z_um <= 669.0)

             # If debug is enabled, and first time this region appears, print BEFORE scaling
             if getattr(self, 'pcc_debug', 0) == 1:
                if not hasattr(Diffusion, 'pcc_debug_printed') and np.any(low_z_mask):
                       print(f"[PCC DEBUG] Found {np.sum(low_z_mask)} hits with z ∈ [663, 669] µm")
                       print("z_um[:10]:", z_um[low_z_mask][:10])
                       print("eps[:10]:", eps[low_z_mask][:10])
                       print("n_obs_e (before):", n_obs_e[low_z_mask][:10])

             n_obs_e = np.round(n_obs_e * eps).astype(int)
             # If debug is enabled, print AFTER scaling (once)
             if getattr(self, 'pcc_debug', 0) == 1:
                if not hasattr(Diffusion, 'pcc_debug_printed') and np.any(low_z_mask):
                     print("n_obs_e (after):", n_obs_e[low_z_mask][:10])
                     Diffusion.pcc_debug_printed = True


        ##########DIFFUSION PROCESS ###################################################################
        ## Only if the process is active in the json file, i.e., if "only_pair_creation is equal" = 0

        if(self.only_pair_creation==0):

            ### GAUSSIAN DISTRIBUTION OF THE E-H PAIR AROUND THE NOMINAL POSITION
            ###########################################################################################
            ### variance proportional to z, allow a value of z, where to start for fully diffusion
            ### for values below z_offset no diffusion applied (sigma=0)
            var_xy=np.zeros_like( hits.z )
            mask_zoffset=hits.z>self.z_offset
            

            var_xy[np.array(mask_zoffset)]=-1.0*self.A*np.log(1.0-self.B*(hits.z[mask_zoffset]-self.z_offset))
            # check for negative values of var_xy
            var_xy[var_xy<=0.0]=0

            # Note that for depth EQUAL to 1000 um the log is -inf, and nan for those depth > 1000.0u
            # however the depth of the CCD, i.e. the silicon bulk will never be as much as that
            if (var_xy==np.inf).any() or (var_xy==np.nan).any():
                print(" ---- ")
                print("\033[33;m DiffuseElectrons :: Informative WARNING\033[1;m")
                print(" This success (eventid+ccdid) has some energy depositons at the maximum thickness of the CCD")
                print("\t The variance for this success has been set to the `maximum variance`, pixel size")
                print("\t The problem is comming from the simulations, where an event at the edge has no sense!")
                print(" ---- ")
                var_xy[np.logical_or(var_xy == np.inf, var_xy == np.nan)] =  cfg.xpixelsize

            ### place all e-h pair created at the nominal position
            nominal_posx_e=np.repeat(hits.x,n_obs_e)
            nominal_posy_e=np.repeat(hits.y,n_obs_e)
            varxy_e = np.repeat(var_xy,n_obs_e)    
            
            sigma_xy_e = np.sqrt(varxy_e)


            # --- Apply precomputed sigma map if repulsion model is on ---
            if getattr(self, 'repulsion', 0) == 1 and getattr(self, 'sigma_map_loaded', False):
                # Energy per electron cloud in keV
                nominal_E = np.repeat(hits.Edep.values, n_obs_e) / 1e3  # eV → keV

                # Depth per electron cloud in microns (from mm), then to meters for the map
                z_um = np.repeat(hits.z.values, n_obs_e) * (u.mm / u.um)  # mm → µm
                z_m  = z_um * 1e-6                                        # µm → m

                try:
                    sigma_from_map = self.sigma_interpolator(
                        np.column_stack([nominal_E, z_m])
                    )  # returns σ in µm if your map stores µm; adjust next line accordingly

                    # Convert map σ to millimeters for the rest of the pipeline.
                    # If your sigma_map values are µm: multiply by (u.um / u.mm) = 1e-3
                    sigma_xy_e = np.asarray(sigma_from_map, dtype=float) * (u.um / u.mm)  # µm → mm
                    #if not hasattr(self, '_repulsion_debug_detailed'):
                     # n_show = min(10, len(sigma_xy_e))
                      #print(f"[Diffusion] Showing first {n_show} (E, z, σ) samples:")
                      #for i in range(n_show):
                       #    print(f"  E = {nominal_E[i]:8.4f} keV,  z = {z_um[i]:8.2f} µm,  σ_xy = {sigma_xy_e[i]*1e3:8.4f} µm")
                    # Debug (prints once to avoid spam)

                    ##DEBUG###

                    #if not hasattr(self, '_repulsion_debug_printed'):
                     #   print(f"[Diffusion] Using σ_map for {len(sigma_xy_e)} e-h pairs "
                      #        f"(E={nominal_E.min():.3f}–{nominal_E.max():.3f} keV, "
                       #       f"z={z_um.min():.1f}–{z_um.max():.1f} µm)")
                       # print(f"[Diffusion] Sample σ_xy (mm): {sigma_xy_e[:5]}")
                       # self._repulsion_debug_printed = True
                except Exception as e:
                    print(f"[Diffusion] Repulsion map interpolation failed: {e}")

            ### add energy dependence
            else:
                #print("no repulsion ")
                ### repulsion ==0 so, execute old energy dependent model ###
                nominal_E = np.repeat(hits.Edep.values,n_obs_e)
                sigma_xy_e_zmax = np.sqrt( -1.0*self.A*np.log(1.0-self.B*(self.z_max-self.z_offset)) )
                beta = self.alpha*(u.pix2mm*(u.um/u.mm)) / sigma_xy_e_zmax
                sigma_xy_e = sigma_xy_e * (self.alpha0 + beta*nominal_E )

            ### add also the sigma_xy as an attribure of the hits object
            hits.sigma_xy = sigma_xy_e

            ### shift each e-h pair created from its nominal position
            hits.x = np.random.normal(nominal_posx_e,sigma_xy_e)
            hits.y = np.random.normal(nominal_posy_e,sigma_xy_e)

        ########NO DIFFUSION PROCESS, if "only_pair_creation" = 1 ###################

        else:
            ### place all e-h pair created at the nominal position
            nominal_posx_e=np.repeat(hits.x,n_obs_e)
            nominal_posy_e=np.repeat(hits.y,n_obs_e)
            hits.x = nominal_posx_e
            hits.y = nominal_posy_e
    
        

        ## the depth for each e-h pair created given an energy lost, is equal 
        ##  to the nominal z of the Edep
        hits.z = np.repeat(hits.z,n_obs_e)

        ### each e-h pair created corresponds to one carried charge, i.e. 1
        hits.Edep = np.ones_like(hits.x)

        ### for each e-h pair, add (if needed) its attributes: PDG and time
        if hasattr(hits,'pdg'):
            pdg_electrons = np.repeat(hits.pdg.values,n_obs_e)
            hits.pdg = pdg_electrons
        if hasattr(hits,'time'):
            time_electrons = np.repeat(hits.time.values,n_obs_e)
            hits.time = time_electrons
        
        ### include the number of electron-hole pairs for each energy deposition
        ### this will be used during pixelization
        hits.N_carried_charges = hits.Edep
        ### go to internal units of Edep, the rest are already in the internal units
        hits.Edep = hits.Edep*u.e2eV
        
        if hits.Edep.size == 0:
            hits.killed = True

    def execute_process_karthik(self,hits):
        """:obj:`Diffusion.execute_process` is a method to model the charge pair creation 
        and charge diffusion for each energy loss given by the numpy.ndarray `hits.Edep`.

        While the parameters of the model :math:`A`, :math:`B` and :math:`z_{offset}` were obtained by
        assuming units of **um**, the positions for hits are in **mm** (default units of geant4).
        **Internally the positions will be converted to um, but stored in mm.**

        Parameters
        ----------
        hits : :obj:`pysimdamicm.utils.data_foramts.G4HitCollection`
            A collection of hits (subset of hits from the geant4
            simulations) with at least the following attributes: time, and
            `__column_axis__`.
    
        Returns
        -------
        hits : :obj:`pysimdamicm.utils.data_foramts.G4HitCollection`
            with actualized values of x,y,z, Edep, time and pdg

            Units (as the inputs ones):
                * Edep in eV
                * position in mm
                * time in seconds
                * pdg adimensional

        Note
        ----
        **hits** must have **x**, **y**, **z** and **Edep** as attributes, if not an
        **AttributeError** will be raised.

        """

        ########################################################################################3 CHARGE GENERARION: KARTHICK MODEL
        self.load_sigma_map_lazy()

        ##########################
        ##### check if electron or nuclear recoil HERE by looking into the PDF
        #####  

        sumEnergy = hits.Edep.sum()
        if sumEnergy > 50.:
            #The number of electrons worth considering varies. Higher energy deposits have 
            # a very wide spread of possible numbers of electrons produced. If n is round(E/e2eV), 
            # we want to calculate what values of m give us p_m(E) / p_n(E) ~= 0.05, and make 
            # our range of test values span this interval, as this would give us a decent sample 
            # of electron numbers that have a reasonable chance of being produced. For very high 
            # energies this could result in hundreds of values being tested, but this is to be expected
            n_eh_avg = np.round(sumEnergy/u.e2eV)
            prob_const = np.sqrt(((n_eh_avg * u.e2eV - sumEnergy)/(u.e2eV * np.sqrt(n_eh_avg * 
                self.fanofactor))) ** 2 - 2 * np.log(0.05))
            upper_bound = np.round((2 * sumEnergy + self.fanofactor * prob_const ** 2 * u.e2eV + prob_const * 
                np.sqrt((self.fanofactor * prob_const * u.e2eV)** 2 + 4 * sumEnergy * self.fanofactor * u.e2eV)) / (2 * u.e2eV))
            lower_bound = np.round((2 * sumEnergy + self.fanofactor * prob_const ** 2 * u.e2eV - prob_const * 
                np.sqrt((self.fanofactor * prob_const * u.e2eV) ** 2 + 4 * sumEnergy * self.fanofactor * u.e2eV)) / (2 * u.e2eV))
            n_eh_range =  np.arange(lower_bound, upper_bound, dtype = int)

            n_eh_probs = (1 / np.sqrt(2 * np.pi * n_eh_range * self.fanofactor)) * np.exp((-1 / (2 * n_eh_range * self.fanofactor)) * 
                    ((n_eh_range * u.e2eV - sumEnergy) / u.e2eV) ** 2)
            n_eh_probs = n_eh_probs / n_eh_probs.sum()

        else:
            sumEnergy = np.round(sumEnergy / 0.05) * 0.05
            if sumEnergy < 1.1:
                sumEnergy = 1.1
            hit_energy = np.searchsorted(self.ehenergies, sumEnergy)
            n_eh_probs = self.ehprobs[hit_energy]
            n_eh_range = np.arange(1, 21)
        
        if np.sum(n_eh_probs) > 0:
            sum_obs_e = np.random.choice(n_eh_range, p = n_eh_probs)
            hits_Edep_np = hits.Edep.to_numpy().flatten()
            electron_Dist = hits_Edep_np / hits_Edep_np.sum()
            indices = np.arange(hits_Edep_np.size)
            ### Two methods for distribution of electrons. As of yet undecided which...
            # Method 1: full numbers of electrons in pixels are certain, the remainder is probabilistic
            n_obs_e = np.floor(sum_obs_e * electron_Dist)
            electron_Remainder = sum_obs_e - n_obs_e.sum(dtype = int)
            if electron_Remainder > 0:
                try:
                    electron_Dist_Remainder = (sum_obs_e * electron_Dist - n_obs_e) / (sum_obs_e * electron_Dist - n_obs_e).sum()
                    added_indices = np.histogram(np.random.choice(indices, size = electron_Remainder, p = electron_Dist_Remainder), bins = np.append(indices, indices.size))[0]
                    n_obs_e = (n_obs_e + added_indices).astype(int)
                except:
                    print("sumEnergy: {}".format(sumEnergy))
                    print("n_eh_range: {}".format(n_eh_range))
                    print("n_eh_probs: {}".format(n_eh_probs))
                    print("sum_obs_e: {}".format(sum_obs_e))
                    print("n_obs_e: {}".format(n_obs_e))
                    print("electron_Remainder: {}".format(electron_Remainder))
        else:
            sum_obs_e = 0
            n_obs_e = np.zeros(hits.Edep.to_numpy().size)


        ########################################################################################3 DIFFUSE EACH e- TO THE SURFACE
        n_obs_e = pd.Series(n_obs_e)
        if self.apply_pcc == 1:
             z_um = hits.z *(u.mm/u.um)  # Convert mm to microns
             #print("[PCC] Applied to hits at z (µm):", z_um[:10])
             eps = self.eps_func(z_um)
             # Identify hits in the PCC transition region
             low_z_mask = (z_um >= 663.0) & (z_um <= 669.0)

             # If debug is enabled, and first time this region appears, print BEFORE scaling
             if getattr(self, 'pcc_debug', 0) == 1:
                if not hasattr(Diffusion, 'pcc_debug_printed') and np.any(low_z_mask):
                       print(f"[PCC DEBUG] Found {np.sum(low_z_mask)} hits with z ∈ [663, 669] µm")
                       print("z_um[:10]:", z_um[low_z_mask][:10])
                       print("eps[:10]:", eps[low_z_mask][:10])
                       print("n_obs_e (before):", n_obs_e[low_z_mask][:10])

             n_obs_e = np.round(n_obs_e * eps).astype(int)
             # If debug is enabled, print AFTER scaling (once)
             if getattr(self, 'pcc_debug', 0) == 1:
                if not hasattr(Diffusion, 'pcc_debug_printed') and np.any(low_z_mask):
                     print("n_obs_e (after):", n_obs_e[low_z_mask][:10])
                     Diffusion.pcc_debug_printed = True
        
        ##########DIFFUSION PROCESS ###################################################################
        ## Only if the process is active in the json file, i.e., if "only_pair_creation is equal" = 0
        if(self.only_pair_creation==0):
            ### GAUSSIAN DISTRIBUTION OF THE E-H PAIR AROUND THE NOMINAL POSITION
            ###########################################################################################
            ### variance proportional to z, allow a value of z, where to start for fully diffusion
            ### for values below z_offset no diffusion applied (sigma=0)
            var_xy=np.zeros_like( hits.z )
            mask_zoffset=hits.z>self.z_offset
            
            var_xy[np.array(mask_zoffset)]=-1.0*self.A*np.log(1.0-self.B*(hits.z[mask_zoffset]-self.z_offset))
            # check for negative values of var_xy
            var_xy[var_xy<=0.0]=0

            # Note that for depth EQUAL to 1000 um the log is -inf, and nan for those depth > 1000.0u
            # however the depth of the CCD, i.e. the silicon bulk will never be as much as that
            if (var_xy==np.inf).any() or (var_xy==np.nan).any():
                print(" ---- ")
                print("\033[33;m DiffuseElectrons :: Informative WARNING\033[1;m")
                print(" This success (eventid+ccdid) has some energy depositons at the maximum thickness of the CCD")
                print("\t The variance for this success has been set to the `maximum variance`, pixel size")
                print("\t The problem is comming from the simulations, where an event at the edge has no sense!")
                print(" ---- ")
                var_xy[np.logical_or(var_xy == np.inf, var_xy == np.nan)] =  cfg.xpixelsize

            ### place all e-h pair created at the nominal position
            nominal_posx_e=np.repeat(hits.x,n_obs_e)
            nominal_posy_e=np.repeat(hits.y,n_obs_e)
            varxy_e = np.repeat(var_xy,n_obs_e) 

            sigma_xy_e = np.sqrt(varxy_e)

            # --- Apply precomputed sigma map if repulsion model is on ---
            if getattr(self, 'repulsion', 0) == 1 and getattr(self, 'sigma_map_loaded', False):
                # Energy per electron cloud in keV
                nominal_E = np.repeat(hits.Edep.values, n_obs_e) / 1e3  # eV → keV

                # Depth per electron cloud in microns (from mm), then to meters for the map
                z_um = np.repeat(hits.z.values, n_obs_e) * (u.mm / u.um)  # mm → µm
                z_m  = z_um * 1e-6                                        # µm → m

                try:
                    # returns σ in µm if your map stores µm; adjust next line accordingly
                    sigma_from_map = self.sigma_interpolator(
                                        np.column_stack([nominal_E, z_m])
                                        )
                   
                    # Convert map σ to millimeters for the rest of the pipeline.
                    # If your sigma_map values are µm: multiply by (u.um / u.mm) = 1e-3
                    sigma_xy_e = np.asarray(sigma_from_map, dtype=float) * (u.um / u.mm)  # µm → mm

                except Exception as e:
                    print(f"[Diffusion] Repulsion map interpolation failed: {e}")
            else:
                ### add energy dependence
                nominal_E = np.repeat(hits.Edep.values,n_obs_e)
                sigma_xy_e_zmax = np.sqrt( -1.0*self.A*np.log(1.0-self.B*(self.z_max-self.z_offset)) )
                beta = self.alpha*(u.pix2mm*(u.um/u.mm)) / sigma_xy_e_zmax
                sigma_xy_e = sigma_xy_e * (self.alpha0 + beta*nominal_E )

            ### add alsot the sigma_xy as an attribure of the hits object
            hits.sigma_xy = sigma_xy_e


            ### shift each e-h pair created from its nominal position
            hits.x = np.random.normal(nominal_posx_e,sigma_xy_e)
            hits.y = np.random.normal(nominal_posy_e,sigma_xy_e)

        else:
            ### place all e-h pair created at the nominal position
            nominal_posx_e=np.repeat(hits.x,n_obs_e)
            nominal_posy_e=np.repeat(hits.y,n_obs_e)
            hits.x = nominal_posx_e
            hits.y = nominal_posy_e
    


        ## the depth for each e-h pair created given an energy lost, is equal 
        ##  to the nominal z of the Edep
        hits.z = np.repeat(hits.z,n_obs_e)

        ### each e-h pair created corresponds to one carried charge, i.e. 1
        hits.Edep = np.ones_like(hits.x)

        ### for each e-h pair, add (if needed) its attributes: PDG and time
        if hasattr(hits,'pdg'):
            pdg_electrons = np.repeat(hits.pdg.values,n_obs_e)
            hits.pdg = pdg_electrons
        if hasattr(hits,'time'):
            time_electrons = np.repeat(hits.time.values,n_obs_e)
            hits.time = time_electrons
        
        ### include the number of electron-hole pairs for each energy deposition
        ### this will be used during pixelization
        hits.N_carried_charges = hits.Edep
        ### go to internal units of Edep, the rest are already in the internal units
        hits.Edep = hits.Edep*u.e2eV
        
        if hits.Edep.size == 0:
            hits.killed = True


###############################################################################################
#####       CONTINUOUS READOUT MODEL
###############################################################################################
class ContinuousReadout(DigitizeProcess):
    """Class to model the continuous readout (electronics process).

    While a pixel is readed, the collected charge (from the previous pixel) 
    suffers a shift along the columns (y-axis), and part of this collected 
    charge is
    transfered to the following pixel different from the pixel where charge
    was collected.

    **List of Attributes:**

    Attributes
    ----------

    pixel_read_time

        Time to read a pixel in seconds, default 0.001

    ampli

        Number of amplifiers used to read the full CCD, default 4

    __column_axis__

        Axes where a shift on the position will be applyied due to readout
        continuous effects, default is y

    __sequence_id__

        Unique identifier for the reconstruction process, default is 20

        Used by `ProcessManager.execute_process` to sort all activated
        process, i.e. included in the full reconstruction process contained
        by `ProcessManager.__sequence__`. 
 
    
    Note
    ----
    Assuming
        * y-axis as the columns of the CCD (horitzontal register)
        * x-axis as the rows (serial register), and 
        * z-axis as the direction of the electric field in the dopped region of the silicon detector.
        * Geant4 units as default units: mm, eV, seconds

    """
    __sequence_id__ = 20
    __column_axis__ = 'y'
    __name__ = 'ContinuousReadout'

    def __init__(self):
        super().__init__()

        # parameters of the model
        # time to read one pixel in seconds
        self.pixel_read_time = float(0.001)*u.s
        self.ampli = int(4)
        self.__init__={'pixel_read_time':u.s,'ampli':1}

    def execute_process(self, hits):
        """Class method to apply a shift on the column position due to
        readout time

        Parameters
        ----------

        hits : :obj:`pysimdamicm.utils.data_foramts.G4HitCollection`
            A collection of hits (subset of hits from the geant4
            simulations) with at least the following attributes: x, y, z and Edep.

        Raises
        ------
        AttributeError
            When hits do not have __column_axis__ as un attribute
                
        """

        if ContinuousReadout.pixel_read_time > 0.0:
            time_zero = 0.0
            pos_time_shift = np.floor( (hits.time - time_zero)/self.pixel_read_time )
        else:
            return
       
        values_axes = (getattr(hits,__column_axis__))
        setattr(hits,__column_axis__, values_axes+pos_time_shift)

        if __verbose__:
            print( "VERBOSE ContinuousReadout.execute_process: \n shift on {0} is {1}".format(
                __axes_shift__,pos_time_shift) )




###############################################################################################
#####       PIXELIZE SIGNAL PROCESS
###############################################################################################
class PixelizeSignal(DigitizeProcess):
    """Class for the pixelization process to convert continuos positions 
    (in mm) to pixel units.

    This is a generalization of a histogram2d function. A histogram divides
    the space into bins, and returns the count of the number of points in
    each bin. `PixelizeSignal.execute_process` allows the computation of 
    the sum, mean, median, or other statistics of a set of values 
    (x,y,z,Edep,[time], [pdg]) within each bin.

    This process must be applied after electron diffusion and continuous
    readout processes (if applicable).

    Assume Geant4 units: mm, eV, seconds

    **List of Attributes:**

    Attributes
    ----------

    N_pixels_in_x

        Number of pixels on the x-axis (columns) within the active region
        of the CCD, default 6000 pixels

    shift_pixel_in_x 

        Number of pixels on the x-axis to skip (non-active silicon region),
        default 0 pixels

    pixel_size_in_x

        Size of the pixel on the x-axis in mm, default 0.015 mm

    N_pixels_in_y

        Number of pixels on the y-axis (rows) within the active region of 
        the CCD, default 6000 pixels

    shift_pixel_in_y

        Number of pixels on the y-axis to skip (non-active silicon region),
        default 0

    pixel_size_in_y

        Size of the pixel on the y-axis in mm, default 0.015 mm


    Note
    ----
    Assuming
        * y-axis as the columns of the CCD (horitzontal register)
        * x-axis as the rows (serial register), and 
        * z-axis as the direction of the electric field in the dopped region of the silicon detector.
    
    So the pixelization is done on the xy-plane, without altering the depth values (x-axis)

    """
    __name__ = 'PixelizeSignal'

    ### Unique identifier for the full reconstruction process.
    __sequence_id__ = 30

    __statistic_function_for_Edep__ = 'sum'

    __statistic_function_for_depth = 'mean'
    __DoDepth = True
    __statistic_function_for_pdg = None
    __DoPDG = True
    __statistic_function_for_time = 'min'
    __DoTime = True


    __Edep_in_e__ = False

    # XXX POSSIBLE value to select part of the e-h pair created XXX
    __charge_threshold = 0.0

    ### others
    #__AXES = ['x','y']
    __AXES = ['y','x']

    def __init__(self):
        super().__init__()

        # parameters of the model
        self.shift_pixel_in_x = int(0)*u.pixel

        self.shift_pixel_in_y = int(0)*u.pixel

        self.unit = "pixel"

        self.__units__={'shift_pixel_in_x':u.pixel,'shift_pixel_in_y':u.pixel}
    
    def get_ccd_shape(self):
        """Read the CCD dimensions from the attribute ccd_shape in singleton Units

            Read as  ccd_shape = (rows,cols) = (y-axis, x-axis)
        """
        self.N_pixels_in_y,self.N_pixels_in_x = u.ccd_shape
        self.pixel_size_in_x = u.ccd_pixel_size_x
        self.pixel_size_in_y = u.ccd_pixel_size_y

    def execute_process(self,hits):
        """Function to model the pixelization of the attributes x,y,Edep,time, and pdg of 
        the `hits` object.

        Parameters
        ----------
        hits : :obj:`pysimdamicm.io.data_foramts.G4HitCollection`
            A collection of hits (subset of hits from the geant4
            simulations) with at least the following attributes: x,y, and
            Edep.

            By default, the function assumes that `hits` also have the
            following attributes pdg, time and z (depth). Unset (if needed) the
            pixeliztion of these optional attributes with
            `PixelizeSignal.__DoPDG`, `PixelizeSignal.__DoTime`, and
            `PixelizeSignal.__DoDepth`, respectively.


            Input Units:
                * position in mm
                * Edep in eV
                * time in seconds
                * pdg adimensional

        Returns
        -------
        hits : :obj:`pysimdamicm.utils.data_foramts.G4HitCollection`
            The same hits object as the input one, but with the attributes
            x,y,Edep,[z,pdg,time] updated to the values of the selected
            statistic in each bin (two-dimensional bin in the xy-plane).
            
            Output units will be the same as the input once
        """
        self.get_ccd_shape()

        #### actualize CCD image size
        #u.ccd_shape=(self.N_pixels_in_x,self.N_pixels_in_y)

        #### From hits.coordinate units of mm to pixel units by using the pixel_size along each axis,
        ####    and obtain the bin_edges for the pixelize process
        coordinate_in_pixel = {}
        bins_in_pixels = {}
        for axes in self.__AXES:
            ### pixel size is in mm!
            pixel_size = getattr(self,'pixel_size_in_'+axes)
            coord_mm_axes  = getattr(hits,axes)
            ### convert units of mm to pixel units by using the size of the pixel (in mm!)
            coordinate_in_pixel[axes] = coord_mm_axes/pixel_size
            # get edges for the pixel grid
            bins_in_pixels[axes] = list(set(coordinate_in_pixel[axes].astype(int)))
            bins_in_pixels[axes].append(max(bins_in_pixels[axes])+2)

        if __verbose__:
            for axes in self.__AXES:
                print("VERBOSE PixelizeSignal.execute_process: pixels in {0} axes {1}, have the follwoing bins {2}".format(
                    axes,coordinate_in_pixel[axes],bins_in_pixels[axes]))

        #### auxiliars values for the pixelization process
        x_sorted = sorted(bins_in_pixels['x'])
        y_sorted = sorted(bins_in_pixels['y'])

        ########################################################################################################
        #
        #   PIXELIZATION PROCESS STARTS FOR ATTRIBUTES
        #       Edep_pixel                  :: total enegy lost per pixel               |__ using the same statistic function
        #       N_carried_charges_pixel      :: total number of e-h pairs per pixel      |
        #       x_pixel                     :: x-axis position in units of pixel
        #       y_pixel                     :: y-axis position in units of pixel
        #       z_pixel                     :: mean depth (x) value in units of mm
        #       time_pixel                  :: energy-weighted time average per pixel
        #       pdg_pixel                   :: the most frequent PDG value per pixel
        #       sigma_xy_pixel              :: XXX first value on the array (mean ?) XXX
        ########################################################################################################

        #### Energy
        #### create a grid to accumulate energy depositions that belongs to the same pixel coordenates
        E_stat, y_stat, x_stat, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                hits.Edep,statistic=self.__statistic_function_for_Edep__,bins=[y_sorted,x_sorted])

        #### FILTERING THE PIXELIZED SIGNAL: only pixels with Edep_pixel > __charge_threshold
        ###############################################################################################
        #### Get COORDENATES for those pixels whose total energy is higher than a min value
        y_inds, x_inds = np.where( E_stat > self.__charge_threshold )
        setattr(hits,'Edep_pixel',E_stat[y_inds, x_inds])
        
        #### N_carried_charges
        #### create a grid to accumulate energy depositions that belongs to the same pixel coordenates
        if hasattr(hits, 'N_carried_charges'):
            Neh_stat, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.N_carried_charges,statistic=self.__statistic_function_for_Edep__,bins=[y_sorted,x_sorted])
            setattr(hits,'N_carried_charges_pixel',Neh_stat[y_inds, x_inds].astype(int))
        
        #### position x and y
        ### x,y in pixel units
        setattr(hits,'x_pixel',x_stat[x_inds].astype(int))
        setattr(hits,'y_pixel',y_stat[y_inds].astype(int))
        
        #### depth (z)
        #### Associate as pixel depth (z-axis) the mean of all z values
        if self.__DoDepth and hasattr(hits,'z'):
            z_stat, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.z, statistic=self.__statistic_function_for_depth, bins=[y_sorted,x_sorted] )
            #### getting only those with charge above the min val
            setattr(hits,'z_pixel',z_stat[y_inds, x_inds])

        ##################################################################################
        ### OTHER STATISTICAL FUNCTIONS TO DEFINE THE VALUES OF PDG/TIME/... PER PIXEL
        ##################################################################################
        ### returns the first value of an array
        get_first_value_of_an_array=lambda l: l[0]
        ### returns the most frequent value of an array
        def get_most_frequent_value(l):
            """Return the most frequent value in the array
            """
            (values,counts)=np.unique(l, return_counts=True)
            ind=np.argmax(counts)
            return values[ind]
        def _has_interaction_from_neutrons(l):
            """Return True if there is any neutron interaction
            """
            if np.any(np.array(l)==2112):
                return 1
            return 0

        def _has_interaction_from_gammas(l):
            """Return True if there is any gamma interaction
            """
            if np.any(np.array(l)==22):
                return 1
            return 0

        def _has_interaction_from_silicon(l):
            """Return True if there is any silicon interaction
            """
            isneutron = False
            if ( np.any( (np.array(l) - 1000000000) >= 140280 ) and np.any( (np.array(l) - 1000000000) <= 140309 ) ) :
                isneutron = True
            if isneutron:
                mfv = get_most_frequent_value(l)
                if mfv < 1000000000:
                    return 0
                else:
                    return 1
            return 0

        def _has_interaction_from_alpha(l):
            """Return True if there is any alpha interaction
            """
            if np.any(np.array(l)== 1000020040):
                return 1
            return 0

        def _has_interaction_from_ion(l):
            """Return True if there is any ion interaction
            """
            if np.any(np.array(l)>1000000000):
                return 1
            return 0

        def _has_interaction_from_electron(l):
            """Return True if there is any ion interaction
            """
            if np.any(np.array(l)==11) or np.any(np.array(l)==-11):
                return 1
            return 0

        #### sigma_xy
        if hasattr(hits,'sigma_xy'):
            sigmaxy_stat, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.sigma_xy, statistic=self.__statistic_function_for_depth, bins=[y_sorted,x_sorted] )
            setattr(hits,'sigma_xy_pixel',sigmaxy_stat[y_inds, x_inds])
        #### PDG
        #### AsSociate as pixel PDG value the mean of all pdg values 
        if self.__DoPDG and hasattr(hits,'pdg'):
            pdg_stat, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.pdg, statistic=get_most_frequent_value, bins=[y_sorted,x_sorted] )
            setattr(hits,'pdg_pixel',pdg_stat[y_inds, x_inds])
        
            #### FLAG TO ALERT ABOUT THE INTERACCION FROM GAMMA OR NEUTRON
            neutron_flag, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.pdg, statistic=_has_interaction_from_neutrons, bins=[y_sorted,x_sorted] )
            setattr(hits,'neutron_pixel',neutron_flag[y_inds, x_inds])
            
            gamma_flag, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.pdg, statistic=_has_interaction_from_gammas, bins=[y_sorted,x_sorted] )
            setattr(hits,'gamma_pixel',gamma_flag[y_inds, x_inds])

            #### FLAG TO ALERT ABOUT THE INTERACTION FROM SILICON, ALPHA, ION
            silicon_flag, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.pdg, statistic=_has_interaction_from_silicon, bins=[y_sorted,x_sorted] )
            setattr(hits,'silicon_pixel',silicon_flag[y_inds, x_inds])
            
            alpha_flag, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.pdg, statistic=_has_interaction_from_alpha, bins=[y_sorted,x_sorted] )
            setattr(hits,'alpha_pixel',alpha_flag[y_inds, x_inds])
        
            ion_flag, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.pdg, statistic=_has_interaction_from_ion, bins=[y_sorted,x_sorted] )
            setattr(hits,'ion_pixel',ion_flag[y_inds, x_inds])
            
            elec_flag, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.pdg, statistic=_has_interaction_from_electron, bins=[y_sorted,x_sorted] )
            setattr(hits,'electron_pixel',elec_flag[y_inds, x_inds])

        #### TIME
        ####     sum( E_{ij}*t_{ij} ) / sum(E_{ij)} where ij is all values at pixel (x_i,y_j)
        ####
        if self.__DoTime and hasattr(hits,'time'):
            ### with the function binned_statistic_2d get the sum of all time (weighted by the
            ###    energy) within a pixel, by using the statistic sum
            time_Edep = hits.Edep * hits.time
            time_stat, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    time_Edep, statistic='sum', bins=[y_sorted,x_sorted])
            ### for each pixel, just divide the sum of the weighted time by the energy for the total
            ###     energy at the pixel
            setattr(hits,'time_pixel',time_stat[y_inds, x_inds]/E_stat[y_inds,x_inds])
            
            ### Add time statistics for full decay chain
            time_stat_dev, y_ledge, x_ledge, pos_bin = binned_statistic_2d(coordinate_in_pixel['y'],coordinate_in_pixel['x'],
                    hits.time, statistic='std', bins=[y_sorted,x_sorted])
            setattr(hits,'std_time_pixel',time_stat_dev[y_inds, x_inds])
            setattr(hits,'time_img_pixel',(time_stat[y_inds,x_inds]/(E_stat[y_inds,x_inds]*u.img_exp_time)).astype(int))

        #### due to diffusion, or other process the values of x, and y can reach the non-active region
        ####   accept only points within the active region of the CCD
        self.pixels_within_active_region(hits)

        #### SET THE VARIABLE TO BE CLUSTERIZED
        hits.__attr_to_clusterize__+="_pixel"

    def pixels_within_active_region(self,hits):
        
        #print(" Restrictions ... x:{},{}, y:{},{}".format(self.shift_pixel_in_x,
        #    self.N_pixels_in_x+self.shift_pixel_in_x,
        #    self.shift_pixel_in_y, self.N_pixels_in_y+self.shift_pixel_in_y))

        ### x-axis: shift_pixel_in_x <= x_pixel <= N_pixels_in_x + shift_pixel_in_x 
        mask_x=np.logical_and(hits.x_pixel>=self.shift_pixel_in_x,hits.x_pixel<=(self.N_pixels_in_x+self.shift_pixel_in_x))

        ### y-axis: shift_pixel_in_y <= y_pixel <= N_pixels_in_y + shift_pixel_in_y
        mask_y=np.logical_and(hits.y_pixel>=self.shift_pixel_in_y,hits.y_pixel<=(self.N_pixels_in_y+self.shift_pixel_in_y))

        ### both requirements are mandatory
        mask=np.logical_and(mask_x,mask_y)
        
        for attr in filter(lambda a:a.count("pixel")>0,hits.__dict__.keys()):
            setattr(hits,attr,getattr(hits,attr)[mask])
        
        if len(hits.x_pixel)==0:
            hits.killed = True


###############################################################################################
#####       DARK CURRENT MODEL
###############################################################################################
class DarkCurrent(DigitizeProcess):
    """Dark Current process.

    Most photodetectors (such as photodiodes, CCD sensors, ...) produce a signal current which is
    more or less proportional to the incident optical power. However, even in the absence of any
    light input, there is often some tiny amount of current, known as dark current. 

    The DC current is often caused by thermal excitation (generation) of carriers --not necessarily
    directly from valence to conduction band, but possible through defect states related to crystal
    defects or impurities. The rate of such thermal processes depends not only on the active area,
    but also critically on the temperature and on the band gap energy of the material, and also on
    the operation voltage (particularly near the breakdown voltage, where impact ionization can
    occur). At high voltages, tunneling through the depletion region may also contribute.
    Dark currents may also be generated by some leakage currents which are not related to thermal
    excitation.
 
    .. math::
       :nowrap:

       \\begin{equation}
         I_{DC} \propto T^{\\frac{3}{2}} exp{ (\\frac{ -E_{g} }{K_{B}T}) }
       \end{equation}


   
    So, DC is a continuous thermal generation of electron-hole pairs in the depletion region. As the
    vias increases the depletion region grows and the dark current also grows. Once the whole wafer
    is depleted the dark current will be remain constant.
 
    Dark current is assumed to follow a poisson distribution where lambda is proportional to the 
    exposure time of the CCD image. 


    **List of Attributes:**

    Attributes
    ----------

        darkcurrent
            
            number of electrons per pixel per day

        exp_time
            
            exposure time of the dataset (time that CCDs have been exposed to ionizing events), in days

    row_readtime
        
            exposure time of a row, in days


        img_size
            
            only for running mode **cropped-image**
            Number of pixels in the x-axis to create the dark current image surrounding the hit collection of the event (squared image)

        min_enlarge
            
            minimum number of pixels to add when the image size is smaller than the cluster elongation.
            In this cases, the size of dark  current image will be <cluster_size> + <min_enlarge>
        
        mode
        
            Name of the running mode to generate the dark current image. Possible values:

                * 'full-image' 
                
                    all events will share the same dark current image with the size of the CCD
                    (in this case img_size is ignored)

                * 'cropped-image'
                
                    for each hit collection, a new DC image will be created, in this case the
                    DC image should be smaller than the real CCD size, and img_size will be used
     exp_mode 

        Method to calculate the mean DC:

                * 'exposure': the DC mean is calculated only on the base of the exp_time. Each pixel has the same exposure time.
     
                * 'exposure-per-row': the DC mean is calculated row by row, taking into account the row_readtime. 

    
    
        n_rows_overscan
            number of extra pixels for the overscan region in rows (x-axis)

        n_rows_prescan
            number of extra pixels for the prescan region in rows (x-axis)

        n_cols_overscan
            number of extra pixels for the overscan region in cols (y-axis)

        n_cols_prescan
            number of extra pixels for the prescan region in cols (y-axis)


    Example
    -------

        >>> import pysimdamicm as dam
        >>> dc = dam.detector_response.DarkCurrent()
        >>> dc.__name__
        'DarkCurrent'
        >>> dc.__sequence_id__
        40
        >>> dc.__units__
        {'img_size': 1, 'mode': 1, 'darkcurrent': 1.1574074074074073e-05, 'exp_time': 86400.0,'min_enlarge': 1}
        >>> dc_params = {'exp_time':8/24., 'darkcurrent': 0.001 , 'mode':''}
        >>> dc.set_parameters(**dc_params)
        Running on mode: full-image
    >>> dc.dark_current_image
    array([[0., 0., 0., ..., 0., 0., 0.],
           [0., 0., 0., ..., 0., 0., 0.],
           [0., 0., 0., ..., 0., 0., 0.],
           ...,
           [0., 0., 0., ..., 0., 0., 0.],
           [0., 0., 0., ..., 0., 0., 0.],
           [0., 0., 0., ..., 0., 0., 0.]])
    >>> dc.dark_current_image.mean()
       0.0012602795833333822


    """
    __name__="DarkCurrent"

    ### Unique identifier for the full reconstruction process.
    __sequence_id__ = 40

    def __init__(self):
        super().__init__()

        ### EXECUTION MODE: full-image or cropped-image
        self.mode = 'cropped-image'
        self.exp_mode = 'exposure' #or exposure-per-row    
        self.rng = np.random.default_rng(seed=self.__seed__)

        # lambda of the poisson distribution in units of electrons/pixel/hour
        self.darkcurrent=float(0.001)*u.e/u.pixel/u.day
        self.exp_time=float(8.0/24.)*u.day
        self.row_readtime=float(86400/86400)*u.day

        ### DEFINE parameters for R.O.I.
        self.img_size=int(300)
        self.min_enlarge=int(30)

       
        self.__units__={'darkcurrent':u.e/u.pixel/u.day,'exp_time':u.day,'img_size':u.pixel, 'row_readtime':u.day,'min_enlarge':u.pixel,"mode":1,"exp_mode":1}
        self.unit = "e/pixel"

    def get_overscan_regions(self):
        # when full image is created (simulate noise with also extra-regions)
        self.n_rows_overscan = u.n_rows_overscan
        self.n_rows_prescan  = u.n_rows_prescan
        self.n_cols_overscan = u.n_cols_overscan
        self.n_cols_prescan  = u.n_cols_prescan


    @property
    def dark_current_image(self):
        return self.__dark_current_image
    @dark_current_image.setter
    def dark_current_image(self,val):
        lam = val[0]*val[1]
        if self.recreate_image or self.full_image_counter == 0:
            self.ccd_shape = list(u.ccd_shape)
            self.get_overscan_regions()
            self.ccd_shape[0] += self.n_rows_overscan+self.n_rows_prescan
            self.ccd_shape[1] += self.n_cols_overscan+self.n_cols_prescan
        #insert here a lambda which is row dependent. Lambda components have increasing values
            if self.exp_mode == "exposure-per-row":
                        lam = lam + val[0] * self.row_readtime * np.arange(self.ccd_shape[0])
        #create image which has reverse axis: CCD X is rows, and CCD Y is columns. The image is transposed later to have correct shape. 
            self.__dark_current_image = self.rng.poisson(lam,size =(self.ccd_shape[1],self.ccd_shape[0]))*u.e2eV
            self.__dark_current_image = np.transpose(self.__dark_current_image)
            #self.__dark_current_image = self.rng.poisson(lam,self.ccd_shape)*u.e2eV
            self.full_image_counter += 1
            self.recreate_image = False
            print(" -- Image dark current shape: ", self.__dark_current_image.shape)
            print("    Generate Dark Current Image: mean={} eV, std={} eV".format(
                round(self.__dark_current_image.mean(),3),round(self.__dark_current_image.std(),3)))
            print(" -- Image dark current: ", self.__dark_current_image)

    def execute_process(self):
        """Process implementation is done for two different running modes:
        
            **full-image**

                see :obj:`DarkCurrent.__execute_process_full_image__`

            **cropped-image**

                see :obj:`DarkCurrent.__execute_process_cropped_image__`
        """
        pass

    def __execute_process_cropped_image__(self,hits):
        """Add simulated Dark current image into the attribure **noise** of hits.

        Running on mode **cropped-image**, i.e. one smaller image will be re-generated 
        for each event.

        This run mode is not appropiate for geant4 simulations under the mode
        'full decay', in this case use **full-image** mode.

        Parameters
        ----------
            
            hits : :obj:`pysimdamicm.io.data_format.G4HitCollection`
        """
        hits.mode = self.mode
        if not hasattr(hits,'roi'):
            hits.get_roi_fixed_length(self.img_size,self.min_enlarge)

        ### generating poisson numbers to fill the DC image
        #E_pixel_dc_noise=self.rng.poisson(self.darkcurrent * self.exp_time,hits.roi['size'])*u.e2eV

        #insert here a lambda which is row dependent. Lambda components have increasing values
        min_crop = min(hits.y_pixel)-hits.roi['shift']['y'][0] #finding the bottom,left corner of cropped image
        if self.exp_mode == "exposure-per-row":
                lam = self.darkcurrent * self.exp_time + self.darkcurrent * self.row_readtime * np.arange(min_crop,min_crop + hits.roi['size'][0])
        if self.exp_mode == "exposure":
                lam = self.darkcurrent * self.exp_time
        E_pixel_dc_noise=self.rng.poisson(lam,size= (hits.roi['size'][1],hits.roi['size'][0]))*u.e2eV
        E_pixel_dc_noise = np.transpose(E_pixel_dc_noise)

        print(" -- Image dark current shape: ", E_pixel_dc_noise.shape)
        print(" -- Image dark current: ", E_pixel_dc_noise)

        if not hasattr(hits,'noise'):
            hits.noise = []
        hits.noise.append(E_pixel_dc_noise)
     
    def __execute_process_full_image__(self,hits):
        """Execute process for mode 'full-image'
        
        In this case, the dark current image will be generated once and 
        used in all processed events.

        Parameters
        ----------
            
            hits : :obj:`pysimdamicm.io.data_format.G4HitCollection`

        """
        hits.mode = self.mode
        if not hasattr(self,'dark_current_image'):
            self.dark_current_image = (self.darkcurrent,self.exp_time)


        ### define the ROI needed to account for the total energy per pixel (hits.get_total_signal)
        ##  that attribute ROI (region of interest) is necessary when the noise image is smaller than
        ##  the total CCD region, and a translation to the smallest size must be done to account for
        ## the total charge per pixel
        ## adding `-1`: arrays index starts at 0, and  6000 pixels represents the pixel 5999
        hits.roi = {'shift': {
                    'x': (min(hits.x_pixel)-1,0), 
                    'y': (min(hits.y_pixel)-1,0)
                    },
                    'size': self.ccd_shape
                    }

        ### the detector has more than one type of noise: recorded at noise attribute
        if not hasattr(hits,'noise'):
            hits.noise = []
        hits.noise.append(self.dark_current_image)
        ### only when fits file is recorded
        hits.darkcurrent=self.darkcurrent
        hits.exp_time   = self.exp_time



###############################################################################################
#####       NOISE MODEL
###############################################################################################
class ElectronicNoise(DigitizeProcess):
    """Electronic Noise process.


    Electronic noise is assumed to follow a gaussian distribution. 


    **List of Attributes:**

    Attributes
    ----------

        pedestal
            
            mean number of electrons per pixel

        sigma
            
            standar desviation for the gaussian distribution

        nSamples

            for skipper CCDs the number of samples that one pixel have been read; in this case, the sigma of the Gaussian becomes 
            
        .. math::
               
           \sigma = \\frac{\sigma}{nSamples}

        img_size
            
            only for running mode **cropped-image**
            Number of pixels in the x-axis to create the electronic noise image surrounding the hits collections of the event (squared image)

        min_enlarge
            
            minimum number of pixels to add when the image size is smaller than the cluster elongation.
            In this cases, the size of dark  current image will be <cluster_size> + <min_enlarge>
        
        mode
        
            Name of the running mode to generate the dark current image. Possible values:

                * 'full-image' 
                
                    all events will share the same dark current image with the size of the CCD
                    (in this case img_size is ignored)

                * 'cropped-image'
                
                    for each hit collection, a new DC image will be created, in this case the
                    DC image should be smaller than the real CCD size, and img_size will be used

    """

    __name__="ElectronicNoise"
    ### Unique identifier for the full reconstruction process.
    __sequence_id__ = 50

    def __init__(self):
        """
        """
        super().__init__()
        ### EXECUTION MODE: full-image or cropped-image
        self.mode = 'cropped-image'

        self.nSamples = int(1)
        self.rng = np.random.default_rng(seed=self.__seed__)

        self.sigma   =float(0.25)*u.e
        self.pedestal=float(0.0)*u.e
        
        self.img_size=int(300)
        self.min_enlarge=int(30)

        self.units = "e/pixel"
        self.__units__={'pedestal':u.e,'sigma':u.e,'img_size':u.pixel,'min_enlarge':u.pixel,'mode':1}

    def get_overscan_regions(self):
        # when full image is created (simulate noise with also extra-regions)
        self.n_rows_overscan = u.n_rows_overscan
        self.n_rows_prescan  = u.n_rows_prescan
        self.n_cols_overscan = u.n_cols_overscan
        self.n_cols_prescan  = u.n_cols_prescan


    @property
    def electronic_noise_image(self):
        return self.__electronic_noise_image
    @electronic_noise_image.setter
    def electronic_noise_image(self,val):
        pedestal,sigma = val
        if self.recreate_image or self.full_image_counter == 0:
            self.ccd_shape = list(u.ccd_shape)
            self.get_overscan_regions()
            self.ccd_shape[0] += self.n_rows_overscan+self.n_rows_prescan
            self.ccd_shape[1] += self.n_cols_overscan+self.n_cols_prescan
            self.__electronic_noise_image =  self.rng.normal(pedestal,sigma/np.sqrt(self.nSamples),self.ccd_shape)*u.e2eV
            self.full_image_counter += 1
            self.recreate_image = False
            print(" -- Image electronic noise shape: ", self.__electronic_noise_image.shape)
            print("    Using: mu_0={} e-/pix, sigma_0={} e-/pix".format(pedestal,sigma))
            print("    Generate Electronic Noise Image: \n\t values for the data sample \n\t mean={} eV,std={} eV".format(
                self.__electronic_noise_image.mean(),self.__electronic_noise_image.std()))

    def execute_process(self,hits):
        """Process implementation is done for two different running modes:
        
            **full-image**

                see :obj:`ElectronicNoise.__execute_process_full_image__`

            **cropped-image**

                see :obj:`ElectronicNoise.__execute_process_cropped_image__`
        """
        pass

    def __execute_process_cropped_image__(self,hits):
        """Add simulated Electronic noise into the attribure **noise** of hits.

        Running on mode **cropped-image**, i.e. one smaller image will be re-generated 
        for each event.

        This run mode is not appropiate for geant4 simulations under the mode
        'full decay', in this case use **full-image** mode.

        Parameters
        ----------
            
            hits : :obj:`pysimdamicm.io.data_format.G4HitCollection`
        """
        hits.mode = self.mode
        if not hasattr(hits,'roi'):
            hits.get_roi_fixed_length(self.img_size,self.min_enlarge)

        ### generating gaussian values as electronic noise
        ###     for skipper CCDs sigma --> sigma/sqrt(nSamples)
        E_pixel_readout_noise=self.rng.normal(self.pedestal,self.sigma/np.sqrt(self.nSamples),hits.roi['size'])*u.e2eV

        if not hasattr(hits,'noise'):
            hits.noise = []
        hits.noise.append(E_pixel_readout_noise)


    def __execute_process_full_image__(self,hits):
        """Execute process for mode 'full-image'
        
        In this case, the electronic noise image will be generated once and 
        used in all processed events.

        Parameters
        ----------
            
            hits : :obj:`pysimdamicm.io.data_format.G4HitCollection`

        """

        hits.mode = self.mode
        if not hasattr(self,'electronic_noise_image'):
            self.electronic_noise_image = (self.pedestal,self.sigma)

        ### define the ROI needed to account for the total energy per pixel (hits.get_total_signal)
        ##  that attribute ROI (region of interest) is necessary when the noise image is smaller than
        ##  the total CCD region, and a translation to the smallest size must be done to account for
        ## the total charge per pixel
        ## adding `-1`: arrays index starts at 0, and  6000 pixels represents the pixel 5999
        hits.roi = {'shift': {
                    'x': (min(hits.x_pixel)-1,0), 
                    'y': (min(hits.y_pixel)-1,0)
                    },
                    'size': self.ccd_shape
                    }

        ### the detector has more than one type of noise: recorded at noise attribute
        if not hasattr(hits,'noise'):
            hits.noise = []
        hits.noise.append(self.electronic_noise_image)

        ### only when fits file is recorded
        hits.pedestal=self.pedestal
        hits.sigma   = self.sigma


###############################################################################################
#####       PIXEL SATURATION
###############################################################################################
class PixelSaturation(DigitizeProcess):
    """
    Pixel saturation is when the incident light at a pixel causes one of the color channels of the
    sensor to respond at its maximum value.

    Saturation and blooming are related phenomena that occur in all CCD image sensors under conditions 
    in which either the finite charge capacity of individual photodiodes, or the maximum charge transfer 
    capacity of the CCD, is reached.

    Once saturation occurs at a charge collection site, accumulation of additional photo-generated
    charge results in overflow, or blooming, of the excess electrons into adjacent device
    structures. 


    """

    __name__="PixelSaturation"
    ### Unique identifier for the full reconstruction process.
    __sequence_id__ = 55

    def __init__(self):
        """
        """
        super().__init__()

        self.saturation = 65535*u.ADC*u.ADC2eV
        self.__units__={'saturation':u.ADC*u.ADC2eV}

    def execute_process(self,hits):
        """Set the attribute `pix_saturation` of the singleton :obj:`pysimdamicm.io.G4Units.Units`
        
        This attribute is used when the total CCD image is computed, see method :obj:`pysimdamicm.io.data_format.get_total_signal()`

        Parameters
        ----------
            
            hits : :obj:`pysimdamicm.io.data_format.G4HitCollection`

        """
        hits.pix_saturation = self.saturation
        setattr(u,"pix_saturation",self.saturation)


###############################################################################################
#####       PIXEL SATURATION
###############################################################################################
class PasteClusterProcess(DigitizeProcess):
    """

    """

    __name__="PasteClusterProcess"
    __sequence_id__ = 60

    __distributions = ['linear','compton','uniform']

    def __init__(self):
        """
        """
        super().__init__()
        
        self.blank_outdir = None

        ### string/pattern to directory with blank images
        self.blank_dir = None
        self.extension = 0
        ### model to randomly generate the new position in rows
        self.row_cls_dist_fx = 'compton'
        ### number of clusters to paste per image
        self.clusters_per_blank = None
        # number of cluster per blank that should be pasted in a uniform way,
        # these clusters are clusters hitting the CCD during exposure time
        self.clusters_per_blank_texp = None
        
        self.in_ADU      = True
        self.calibration = np.nan
        self.sigma       = None

        # ADD parameter to control randomness on the selection of blank images
        self.random_image_start = True

        self.__units__ = {'blank_outdir':1,'blank_dir':1,'row_cls_dist_fx':1,'clusters_per_blank':1,'extension':1,'in_ADU':1,
                'calibration':1,'sigma':1, 'clusters_per_blank_texp':1, 'random_image_start':1}

    def execute_process(self,hits):
        """

        Parameters
        ----------
            
            hits : :obj:`pysimdamicm.io.data_format.G4HitCollection`

        """
        if not hasattr(self,'_blanks'):
            # Asserts
            if not self.row_cls_dist_fx in PasteClusterProcess.__distributions:
                raise NotImplementedError("The distribution {} is not implemented.".format())
            if self.clusters_per_blank <1:
                raise AttributeError("The number of cluster per image must be >0.")
            #fdir = os.path.dirname(self.blank_dir)
            #if not os.path.isdir(fdir):
            #    raise OSError("Not found directory for blank images {}".format(fdir))
            if not os.path.isdir(self.blank_outdir): 
                raise OSError("Directory {} does not exist".format(self.blank_outdir))

            ### iterator object BlankImages
            self._blanks = BlankImages(self.blank_dir,self.clusters_per_blank,self.blank_outdir,self.extension,
                    self.in_ADU,self.row_cls_dist_fx,calibration=self.calibration,sigma=self.sigma,
                    Ncls_flat=self.clusters_per_blank_texp,random_image_start=self.random_image_start)

        _ = next(self._blanks)
        self._blanks.append_clusters(hits)

    def __del__(self,hits=None):
        """
        """
        # save current image before killing the object
        try:
            self._blanks.save_blank_image()
        except AttributeError:
            print("Nothing to delete, _blanks does not exists!")


########################################################################################
#
#
########################################################################################
class BlankImages(object):
    """
    """
    def __init__(self,filenames_pattern,required_clusters_per_blank,outdir,extension=0,
            in_ADU=True,row_dist_fx='compton',calibration=np.nan,sigma=np.nan,Ncls_flat=0,
            random_image_start=True):
        
        if in_ADU:
            # define calibration with e2eV and ADC2eV (both given by the user)
            self.cal = u.e2eV/u.ADC2eV
        else:
            self.cal = 1.0
        
        self.isgauss = False
        if sigma is not None:
            self.sig     = sigma
            self.isgauss = True

        if not np.isnan(calibration):
            self.isgauss = True
            self.cal = calibration
        
        self._charge_conversion_function  = lambda isgauss: (1/u.e2eV)*np.random.normal(self.cal,self.sig) if isgauss else self.cal

        # initialization of the conversion factor
        self._charge_conversion = self._charge_conversion_function(self.isgauss)
        # set distribution function for the rows
        self.row_dist_fx = row_dist_fx

        # output directory to store the blank images after appending clusters
        self.outdir = outdir

        # list of files in the given directory following the pattern name
        #   randomly shuffle the list of files to be sued as blanks 
        self._filelist = sorted(glob(filenames_pattern))

        self._len  = len(self._filelist)
        if self._len==0:
            msm = "There is no blank images at {}".format(filenames_pattern)
            raise RuntimeError("\x1b[31m {} \x1b[m".format(msm))

        # number of clusters pasted into the blank
        self._nmax_clusters = required_clusters_per_blank
        # initialization of counter and current index (pointing to an image from the list)
        self._appended_clusters = 0
        if random_image_start:
            random.shuffle(self._filelist)
            self._current_index = np.random.randint(self._len)
        else:
            self._current_index = 0
        
        # number of clusters pasted into the blank that should follow a flat distribution
        self.Ncls_flat = Ncls_flat

        # read the blank data, and keep the extension (needed to read the header before storing final image)
        self._extension = extension
        self._blank = fits.getdata(self._filelist[self._current_index % self._len],self._extension)#.astype(int)
        # monte carlo truee image, any cluster pasted on the _blank images, will be also pasted in
        # this "zero-charged" image, being the monte carlo true simulation
        # this will also be recorded in the output fits file image
        self._MCT      = np.zeros_like(self._blank).astype(float)
        self._MCT_ID   = np.zeros_like(self._blank).astype(int)
        self._MCT_EMCT = np.zeros_like(self._blank).astype(float)
        self._MCT_ERAW = np.zeros_like(self._blank).astype(float)

    def __next__(self):
        """
        """
        if self._appended_clusters == self._nmax_clusters:
            # save image before change the iter
            self.save_blank_image()
            # increase iter to the next image
            self._current_index += 1
            # initialize the clusters appended in the image
            self._appended_clusters = 0
            self._blank = fits.getdata(self._filelist[self._current_index % self._len],self._extension) #.astype(int)
            # initialize also the monte carlo true image
            self._MCT   = np.zeros_like(self._blank).astype(float)
            self._MCT_ID= np.zeros_like(self._blank).astype(int)
            self._MCT_EMCT = np.zeros_like(self._blank).astype(float)
            self._MCT_ERAW = np.zeros_like(self._blank).astype(float)
            # restart calibration constant
            self._charge_conversion = self._charge_conversion_function(self.isgauss)
        return self._blank
    
    def save_blank_image(self,overwrite=True):
        """
        """
        image_file_name = self._filelist[self._current_index % self._len]

        # get header from blank image
        header = fits.getheader(image_file_name,self._extension)
        # add some parameters
        header['NMAX'] = (self._nmax_clusters,"Maximum number of evets")
        header['NEVTS'] = (self._appended_clusters,"Number of appended events")
        header['eV2ADU'] = (self._charge_conversion*u.e2eV,"Conversion used: ADU/e-")
        # create the HDU list: using header from blank
        hdul = fits.HDUList([fits.PrimaryHDU(data=self._blank, header=header)])
        hdul.append(fits.ImageHDU(data=self._MCT,name="MonteCarloTruth"))
        # append another extension with the cluster id
        hdul.append(fits.ImageHDU(data=self._MCT_ID,name="MCT_ID"))
        # append another extension with the energy/charge from the blank image before cluster has been injected
        hdul.append(fits.ImageHDU(data=self._MCT_ERAW,name="MCT_ERAW"))
        # append another extension with the total energy of the cluster
        hdul.append(fits.ImageHDU(data=self._MCT_EMCT,name="MCT_EMCT"))

        # save image
        fitsname = "{}/psimulCCDimg_{}".format(self.outdir,
                os.path.basename(image_file_name).replace(".fits","_Ncls{}".format(int(self._appended_clusters))))
        _rhash = hash(hash(time()) + abs(hash(self)))
        fitsname += "_{}.fits".format(_rhash)

        hdul.writeto(fitsname,overwrite=overwrite)


    def append_clusters(self,hits):
        """
        """
        killed,(new_rows,new_cols,Edep) = self._random_position(hits)
        if killed:
            return
        ## round energy to the nearest integer (ADU units is an integer variable)
        ############### APPLY QUADARTIC CALIBRATION ---- COMPTON
        if self.row_dist_fx in ["compton"]:
            # FIXME --- this is hardcoded for compton
            quadratic_cal = np.polynomial.polynomial.Polynomial([0.4088,5.177,7.85e-6])
            Edep_e = Edep/u.e2eV
            Edep_ADCs = np.round(quadratic_cal(Edep_e)).astype(int)
        else:
            #print(np.sum(Edep))
            #print(self._charge_conversion)
            Edep_ADCs = (Edep*self._charge_conversion).astype(float)
            #print(np.sum(Edep_ADCs))
            #input("press enter ...")
        
        # before appending the clusters, sum up the energy already in the image
        self._MCT_ERAW[new_rows,new_cols] += self._blank[new_rows,new_cols]

        # Append the cluster into the current blank
        self._blank[new_rows,new_cols] += Edep_ADCs
        # Update the number of clusters appendend into the current blank
        self._appended_clusters += 1

        # add event pixels also to the monte carlo true image (in units of eV)
        self._MCT[new_rows,new_cols] += Edep
        self._MCT_ID[new_rows,new_cols]   = int(self._appended_clusters)
        self._MCT_EMCT[new_rows,new_cols] = float(np.sum(Edep))
        #print(np.sum(Edep), self._MCT_EMCT[new_rows,new_cols])
        #print(self._blank[new_rows,new_cols] )

        # Apply saturation if is invoked, pix_saturation is in units of eV
        if hasattr(hits,"pix_saturation") and getattr(hits,"pix_saturation")>0:
            max_charge_ADU = hits.pix_saturation/u.e2eV * self._charge_conversion
            self._blank[self._blank>max_charge_ADU] = max_charge_ADU
            hits.Npix_sat = sum(self._blank>max_charge_ADU)

    def _random_position(self,hits,axis=0):
        """

        Parameters
        ----------

            axis : int
                Axis (0:rows, 1:columns) in which the displacement should be applied
                
        """
        # flag to control if we end-up with pixels in the new position
        killed = False

        ### to paste the cluster inside the active region
        prescan  = [u.n_rows_prescan,u.n_cols_prescan]
        overscan = [u.n_rows_overscan,u.n_cols_overscan]
#        print(prescan,overscan)

        # real active region size of the blank image
        active_region = [ self._blank.shape[0]-prescan[0]-overscan[0], self._blank.shape[1]-prescan[1]-overscan[1] ]
#        print(active_region)

        displacement = [0,0]
        # simulated positions
        positions = [hits.y_pixel,hits.x_pixel]
        # event elongation
        Dcol = int(max(hits.x_pixel) - min(hits.x_pixel))
        Drow = int(max(hits.y_pixel) - min(hits.y_pixel))
        elongation = [Drow,Dcol]
        
        if self.row_dist_fx in ['linear']:
            row_max = active_region[axis]
            nrow = np.sqrt(np.random.uniform(0.,1.0)) * row_max

        elif self.row_dist_fx in ['uniform']:
            row_max = active_region[axis]
            nrow = np.random.uniform(0.,row_max)

        else:
            # number of rows in the active region, the only place where we should paste clusters
            row_max = active_region[axis]

            # clusters during exposure time should follow a flat distriution. The number of cluster
            # per image during exposure time are given by: self.Ncls_flat
            if self.Ncls_flat>0 and self._appended_clusters<self.Ncls_flat:
                # clusters that occurs during exposure time
                nrow = np.random.uniform(0,row_max)
            else:
                # clusters during exposure time
                # ratio between first row and last row: how much frequent should be 1(last row) from 0 (first row)
                # from the normalized data this should be
                m = float(6.82)
                ### starts the random algorithm
                s_rm = 1-(1/m)**2
                val = (1-random.uniform(0,s_rm))**.5
                nrow = row_max*(m*val-1)/(m-1)

        # actualize the displacement
        displacement[axis] = int(nrow)
        # move simulated positions into the blank image (that maybe smaller than the simulated one)
        positions[axis] = (positions[axis] - positions[axis].max()) + prescan[axis] +  displacement[axis]
        

        # For axis without displacement: 
        # r = (self._blank.shape[1-axis] - prescan[1-axis] - overscan[1-axis]) / float(u.ccd_shape[1-axis])
        # displacement[1-axis] = int(positions[1-axis].min()*r) - elongation[1-axis]
        Ldim = self._blank.shape[1-axis] - prescan[1-axis] - overscan[1-axis] - elongation[1-axis]
        Lstart = prescan[1-axis]
        displacement[1-axis] = random.randint(Lstart,Ldim)
        positions[1-axis] = positions[1-axis]-positions[1].min() + displacement[1-axis]
        
        # this translation can lead to negative pixels: clusters are killed partially
        try:
            new_rows,new_cols,E = zip(*list(filter(lambda p: p if (p[0]>=0 and p[1]>=0) else None, zip(positions[0],positions[1],hits.Edep_pixel)))) 
        except ValueError:
            hits.killed = True
            new_rows,new_cols,E = None,None,None
            ## In some rare cases, all pixels will be placed in such a way that none of them will fall inside the area
            killed=True
        # remove pixels out from blank region: cluster not fully appended
        try:
            new_rows,new_cols,E = zip(*list(filter(lambda p: p if (p[0]<self._blank.shape[0] and p[1]<self._blank.shape[1]) else None,
                zip(new_rows,new_cols,E))))
        except ValueError:
            ## In some rare cases, all pixels will be placed in such a way that none of them will fall inside the area
            new_rows,new_cols,E = None,None,None
            hits.killed = True
            killed=True

        
        # list objects: convert into arrays
        E = np.array(E)
        new_rows = np.array(new_rows)
        new_cols = np.array(new_cols)

        #raise
        # Add attributes to hits to dump into the ouptut ROOT file
        setattr(hits,'rnd_rows',new_rows)
        setattr(hits,'rnd_cols',new_cols)
        setattr(hits,'rnd_Edep',E)
        
        return killed,(new_rows,new_cols,E)

