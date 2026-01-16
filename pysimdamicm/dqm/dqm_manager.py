### DATA QUALITY MONITOR PROCESS
from pysimdamicm.dqm.me import *

### SINGLETON: UNITS CLASS WITH GLOBAL PARAMETERS
from pysimdamicm.utils.units import Units
from pysimdamicm.utils.report import LaTex
u=Units()

#### other python packages externals to pysimdamicm
from matplotlib import pyplot as plt
import pandas as pd
import pickle as cpickle
import dateutil

### Send summary to server
#import socket
#import struct
import pickle

import json

### warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=Warning)

# time with respect DAY 0 (new year 2022)
DAY0 = datetime.strptime('2022-01-01 00:00:00','%Y-%m-%d %H:%M:%S')
get_deltatime_days = lambda t: (t.to_pydatetime()-DAY0).total_seconds()/60./60./24.
EVL_PLOTS = ['MEFitDCLambda','MEFitDCMu0','MEFitDCSigma','MEFitDCCalibration','MEOVSPedestalMu','MEOVSPedestalSigma']


###############################################################################################
#####       RESPONSE OF THE DETECTOR (FULL RESPONSE,. i.e. all process are account for)
#####       AND (CLUSTER) RECONSTRUCTION
###############################################################################################
class ActivateMEProperty(object):
    """Descriptor to create dynamically the activate process properties.

    For each :obj:`DigitizeProcess` added to the `ProcessManager` class, an attribute
    `activate_<DigitizeProcess.__name__>` must be included to active or
    unactivate such digitize/reconstruc process during the configuration
    of the data process manager.

    The creation of such attribute is done dynamically by using the method
    `__set__` of this class.

    Parameters
    ----------
    manager : :obj:`ProcessManager`

    process_name : str
        The name of the DigitizeProcess given by `DigitizeProcess,.__name__`
    """
    def __init__(self,dqm_manager,me_process_name):
        #### if the process do not exist, KeyError will raise
        self._me_process_class = dqm_manager.__valid_me_processes__[me_process_name]

    def __get__(self,inst_dqm_manager,class_type):
        return len(list(filter(lambda p: isinstance(p,self._me_process_class),inst_dqm_manager.__sequence__))) == 1

    def __set__(self,inst_dqm_manager,value):
        # Check if the process is in the sequence (needed in the two cases)
        # If not None is returned
        me_process = inst_dqm_manager.get_process_from_sequence(self._me_process_class)
    
        ### remove existing process, to be sure the new instanciation does not contain previous
        #       attributes
        if me_process is not None:
            inst_dqm_manager.__sequence__.remove(me_process)

        if value == True:
            # If process exists, create a new brand instance
            me_process = self._me_process_class()
            inst_dqm_manager.__sequence__.append(me_process)
        elif type(value)!=bool:
            raise RuntimeError("Not valid value (only bool accepted)")
        #### (re-)sort the remaining sequence
        inst_dqm_manager.__sequence__.sort(key=lambda p: p.__sequence_id__)

class DQMManager(object):
    """
    Class to manage (activate/deactivate, configure, ...) the full set of monitor elements
    of a CCD setup.

    """
    __sequence__ = []
    ### REMEMBER whenever you are implementing a new process:
    ### Any new ME that can be used by this manager class
    ### must be enumerated here
    __valid_me_processes__ = {
            MEQmeanDiffSingleSkips.__name__:MEQmeanDiffSingleSkips,
            MEReadoutnoiseNDCM.__name__:MEReadoutnoiseNDCM,
            MEMeanPixelChargePerRow.__name__:MEMeanPixelChargePerRow,
            MESTDPixelChargePerRow.__name__:MESTDPixelChargePerRow,
            MESlopeFromMeanPCPerRow.__name__:MESlopeFromMeanPCPerRow,
            MEInterceptFromMeanPCPerRow.__name__:MEInterceptFromMeanPCPerRow,
            MEMedianPixelChargePerRow.__name__:MEMedianPixelChargePerRow,
            MEMADPixelChargePerRow.__name__:MEMADPixelChargePerRow,
            MEMedianPixelChargePerCol.__name__:MEMedianPixelChargePerCol,
            MEMADPixelChargePerCol.__name__:MEMADPixelChargePerCol,
            MEOVSMedianPixelChargePerRow.__name__:MEOVSMedianPixelChargePerRow,
            MEOVSMADPixelChargePerRow.__name__:MEOVSMADPixelChargePerRow,
            MEOVSPCD.__name__:MEOVSPCD,
            MEPCD.__name__:MEPCD,
            MEOVSPedestalMu.__name__:MEOVSPedestalMu,
            MEOVSPedestalSigma.__name__:MEOVSPedestalSigma,
            MEOVSPedestalMuPerRow.__name__:MEOVSPedestalMuPerRow,
            MEOVSPedestalSigmaPerRow.__name__:MEOVSPedestalSigmaPerRow,
            MEDefectsMask.__name__:MEDefectsMask,
            MENpixQ1ePerCol.__name__:MENpixQ1ePerCol,
            MESTDSkipsPerCol.__name__:MESTDSkipsPerCol,
#            MEMaskedPixels.__name__:MEMaskedPixels,
#            MESinglePED.__name__:MESinglePED,
            MEFitDC.__name__:MEFitDC,
            MEFitDCMu0.__name__:MEFitDCMu0,
            MEFitDCSigma.__name__:MEFitDCSigma,
            MEFitDCLambda.__name__:MEFitDCLambda,
            MEFitDCCalibration.__name__:MEFitDCCalibration,
            MEFitMu0PerRow.__name__:MEFitMu0PerRow,
            MEFitDCPerRow.__name__:MEFitDCPerRow,
            MEFitSigmaPerRow.__name__:MEFitSigmaPerRow,
            MEFitCalibrationPerRow.__name__:MEFitCalibrationPerRow,
            MEFitMu0PerCol.__name__:MEFitMu0PerCol,
            MEFitDCPerCol.__name__:MEFitDCPerCol,
            MEFitSigmaPerCol.__name__:MEFitSigmaPerCol,
            MEFitCalibrationPerCol.__name__:MEFitCalibrationPerCol,
#            MEColTransient.__name__:MEColTransient,
#            MEColTransientMu.__name__:MEColTransientMu,
#            MEColTransientAmplitude.__name__:MEColTransientAmplitude,
#            METempShift.__name__:METempShift,
#            MEPresShift.__name__:MEPresShift,
            MESkewnessCoeff.__name__:MESkewnessCoeff,
            MECCDEqualizedImage.__name__:MECCDEqualizedImage,
            MECCDMeanImage.__name__:MECCDMeanImage,
            MECCDStdImage.__name__:MECCDStdImage,
            MECCDCalImage.__name__:MECCDCalImage,
#            MEBaselineShift.__name__:MEBaselineShift,
#            MEHorizontalClusters.__name__:MEHorizontalClusters,
#            MENoisyImages.__name__:MENoisyImages,
#            MESigmaCutoffNoise.__name__:MESigmaCutoffNoise,
#            MEPoorQualityData.__name__:MEPoorQualityData
            }

    ##### module DEBUG VARIABLES
    __verbose__ = False

    def __init__(self,level_of_debug=1):
        # INITIALIZE THE MONGO DB DICTIONARY
        setattr(self,'mongodb',{})

        # create the process activation properties (XXX Must be created to the class)
        for pname,me_process_class in self.__valid_me_processes__.items():
            setattr(DQMManager,'active_{}'.format(pname),ActivateMEProperty(DQMManager,pname))
        
        self.debug_level=level_of_debug

    def get_process_from_sequence(self,me_process_class):
        """
        Parameters
        ----------
        me_process_class : :obj:`pysimdamicm.dqm.ME`
            an instance class of such process class

        Return
        ------
        process_isnt : :obj:`pysimdamicm.dqm.ME`
            the instance of the input class from pysimdamicm.dqm_manager.__sequence__
        """
        try:
            return list(filter(lambda p: isinstance(p,me_process_class),self.__sequence__))[0]
        except IndexError:
            return None
        
    def set_configuration(self,config):
        """Function to set the parameters for the full set of processes
        """        
        for pname in config.process_to_sim:
            if hasattr(self,"active_{}".format(pname)):
                setattr(self,"active_{}".format(pname),True)
            self.set_config_process( pname, config.configuration[pname])

    def set_config_process(self, me_process_name, model_parameter):
        """
        Function to set the parameters of an specific process

        Parameters
        ----------
        me_process_name: str
            Class name for the process to be update the parameters of the model

        Raises
        ------
        KeyError: If the process name is not present in __valid_me_processes__

        RuntimeError: Process not in the current sequence
        """
        me_process_class = self.__valid_me_processes__[me_process_name]
        me_process_inst  = self.get_process_from_sequence(me_process_class)
        if me_process_inst is not None:
            me_process_inst.set_parameters(**model_parameter)
        else:
            print('InstanceError: No instance of "{}" class is found at <DQMManager>.__sequence__'.format(me_process_name))
        
        ## (re-)sort the process in the __sequence__ (in case the user changed its sequence_id value)
        self.__sequence__.sort(key=lambda p: p.__sequence_id__)

    def execute_process(self,rawdata,is_image_HR=False):
        """ This method will execute all process that have been added to self.__sequence__
        through the configuration JSON file. Any active process in this file is included in the
        process list (private data memeber self.__sequence__).
        
        Each process has an identifier (<process>.__sequence_id__) that is used to sort the 
        processes in the processing chain (self.__sequence__). Before starting the execution the 
        processes are sorted using that identifier, and then executed sequentially.  All that is
        done internally, except for processes that can be executed at any time in the processing
        chain. For this process the user can modify the data member __sequence_id__ from the
        configuration JSON file.

        """
        
        for k,me_process in enumerate(self.__sequence__):
            n_amp = -1
            for amp in sorted(rawdata.amplifier.keys()):
                rawdata.set_amplifier(amp)
                me_process.execute_process(rawdata)

                n_amp+=1
                self.append_ME_to_mongodb(me_process,rawdata,n_amp)

                if not me_process._per_amp:
                    # only one time per amplifier (like Compress image, ...)
                    break

    def execute_plot(self):
        """ This method will execute all process that have been added to self.__sequence__
        through the configuration JSON file. Any active process in this file is included in the
        process list (private data memeber self.__sequence__).
        
        Each process has an identifier (<process>.__sequence_id__) that is used to sort the 
        processes in the processing chain (self.__sequence__). Before starting the execution the 
        processes are sorted using that identifier, and then executed sequentially.  All that is
        done internally, except for processes that can be executed at any time in the processing
        chain. For this process the user can modify the data member __sequence_id__ from the
        configuration JSON file.

        """
        print("\nDQM INFO. Running plots for all ME: ")
        for k,me_process in enumerate(self.__sequence__):
            print("     * ", me_process.__name__)
            me_process.execute_plot()

    def build_report(self,run_number,me_outdir,abstract=None,dopdf=True,file_name=None):
        """
        """
        
        print("\nDQM INFO. Create latex file")
        
        if file_name is not None:
            doc_title=r"DQM Report for {}".format(file_name.replace("_","\_"))
        else:
            doc_title=r"DQM Report for run number {}".format(run_number)

        doc_author=r"{}".format(self.__module__.replace("_","\_"))
        self.report = LaTex(doc_title,doc_author)

        for k,me_process in enumerate(self.__sequence__):
            try:
                fig_name = me_process.fig_png_file
            except AttributeError:
                #### No figure to add
                continue
            #### this should be improved to add the results of the TEST and ALARMs
            try:
                fig_caption = me_process.caption
            except AttributeError:
                fig_caption = me_process.title

            if not type(fig_name) is list:
                fig_name = [fig_name]
            for fig in fig_name:
                self.report.add_figure(fig,fig_caption,short_caption=me_process.__name__)

        if file_name is not None:
            texoutfile = '{}/me_report_{}.tex'.format(me_outdir,file_name)
        else:
            texoutfile = '{}/me_report_run{:03}.tex'.format(me_outdir,int(run_number))

        if abstract is not None:
            self.report.abstract = abstract

        # if booked, create the PDF
        if dopdf:
            print("\nDQM INFO. Generate the PDF report file ")
            self.report.build_pdf(texoutfile)


    def dump_run_summary(self,run_number,me_outdir):
        """Serializing the monitor elements with cPickle (now to a file, later on into a mongoDB)
        """
        me_summary_run = {
                'run':run_number,
                'me_dc_gain':u.me_dc_gain,
                'me_dc_lambda':u.me_dc_lambda,
                'me_dc_sigma':u.me_dc_sigma,
                'me_dc_mu0':u.me_dc_mu0
                }
        for k,me_process in enumerate(self.__sequence__):
            me_summary_run[me_process.__name__] = {}
            for attr in ['x','y','fit_model','xlabel','ylabel','title','amplifier']:
                try:
                    val = getattr(me_process,attr)
                except AttributeError:
                    continue
                if attr == 'fit_model':
                    val = val.coef
                me_summary_run[me_process.__name__][attr] = val
        ### send me_summary_run to the server
        #try:
        #    self.send_to_server(me_summary_run)
        #except ConnectionRefusedError:
        #    print("\nDQM INFO. ConnectionRefusedError: Couldnt connect to Server. Is it open?")
        #### storing as pickle object
        out_file_name = '{}/me_summary_run{:03}.pkl'.format(me_outdir,int(run_number))
        file_object = open(out_file_name, 'wb')
        cpickle.dump(me_summary_run, file_object)
        file_object.close()
        print("DQM INFO. Dump all ME into a pickle file: ", out_file_name)

    @property
    def me_summary_ref(self):
        try:
            return self._me_summary_ref
        except AttributeError:
            return None
    @me_summary_ref.setter
    def me_summary_ref(self,val):
        """ Loading ME as reference and add the ref_model
        Assuming all models are poly1d!
        """
        me_summary_ref_dict = self.load_run_summary(val)
        self._me_summary_ref = MERef(me_summary_ref_dict)
        for k,me_process in enumerate(self.__sequence__):
            if me_process.__name__ in me_summary_ref_dict.keys():
                ### For the qtest
                setattr(me_process,'me_ref', me_summary_ref_dict[me_process.__name__])
                # Add the reference run number to me_ref object
                me_process.me_ref.update({'run':me_summary_ref_dict['run']})
                try:
                    val = me_summary_ref_dict[me_process.__name__]['fit_model']
                    print("     * {}: {}".format(me_process.__name__,val))
                    setattr(me_process,'ref_model', np.poly1d(val))
                except KeyError:
                    ### model was not recorded
                    continue

    def load_run_summary(self,run_ref):
        """To de-serialize the me_summary _run data stream
        """
        file_object = open(run_ref, 'rb')
        run_summary = cpickle.load(file_object)
        file_object.close()
        return run_summary

    ####################################################################################################
    #   related to mongoDB
    ####################################################################################################
    #def send_to_server(self, me_summary_run):
    #    """ Prepares the connection to the give ip and send the ME summary dict to 
    #        the server at the given ip.
    #    """
    #    ### Prepares the connection to server
    #    client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    #    ### This has to be changed to whatever ip the web will use
    #    host_ip, port = "127.0.1.1", 9999
    #    ## Connecting to the server ip
    #    client_socket.connect((host_ip, port))
    #    data = b""
    #    payload_size = struct.calcsize("Q")
    #    ### Sends the pickle to server 
    #    me_pkl = pickle.dumps(me_summary_run)
    #    message = struct.pack("Q",len(me_pkl))+me_pkl
    #    client_socket.sendall(message)
    #    client_socket.close()
    #    print( "DQM INFO. All ME sent to server (ip:port): {}:{}".format(host_ip,port))


    def append_ME_to_mongodb(self,me,rawdata,n_amp):
        """Add all information from the ME to the obj:`self.mongodb` pandas object
        """
        # if MECCDImage fill only once
        if me.__name__.count('MECCD')>0 and n_amp>0:
            return
        if me.failed:
            return

        me_data = me.get_data_as_pandas(rawdata)
        run = rawdata.n_run

        # ADD or APPEND depending on the existency of the ME to the given run
        #  self.mongodbself.mongodb is a dictionary with (run_number: DataFrame)
        if not run in self.mongodb.keys():
            self.mongodb[run] = me_data
        else:
            self.mongodb[run] = self.mongodb[run].append( me_data, ignore_index=True )

    def fill_mongoDB_documents(self):
        
        # initialize document for the mongoDB
        self.document = {}

        for run in set(self.mongodb.keys()):
            merun = self.mongodb[run]
            #print(" ADDING RUN ",)
            # initialize run document one per run
            if not run in self.document:
                #print("Initialize run document")
                self.document[run] = {
                    'number'        : run,
                    'tag'           : merun.tag.values[0],
                    'start'         : merun.start.min().to_pydatetime(),
                    'end'           : merun.end.max().to_pydatetime(),
                    'nskips'        : ",".join(map(str,list(set(merun.nskips.values.tolist())))),
                    'slow_control'  : { 'T': merun['T'].mean(), 'T_std': merun['T'].std(), 'P': merun['P'].mean(), 'P_std': merun['P'].std() },
                    'status'        : all(merun.status.values),
                    'Images'        : [],
                    'MEs'           : []
                }

            for n_image in set(merun.image.values):
                meimage = merun[merun.image==n_image]
                #fimage = ';'.join(list(set(meimage.fimage.values)))
                self.document[run]['Images'].append({
                    'number'        : n_image,
                    #'fnumber'       : fimage,
                    'run_number'    : run,
                    'start'         : meimage.start.min().to_pydatetime(),
                    'end'           : meimage.end.max().to_pydatetime(),
                    'status'        : all(meimage.status.values),
                    'nskips'        : ",".join(map(str,list(set(merun.nskips.values.tolist())))),
                    'npbin'         : meimage.npbin,
                    'nsbin'         : meimage.nsbin,
                    'texp'          : meimage.texp,
                    'tread'         : meimage.tread,
                    'nrows'         : meimage.nrows,
                    'ncols'         : meimage.ncols,
                    'slow_control'  : { 'T': meimage['T'].mean(), 'T_std': meimage['T'].std(), 'P': meimage['P'].mean(), 'P_std': meimage['P'].std() },
                    'CCDs'          : [],
                    'MEs'           : []
                })

                for n_ccd in set(meimage.ccd.values):
                    #print("for for ccd ", n_ccd)
                    meccd = meimage[meimage.ccd==n_ccd]

                    self.document[run]['Images'][-1]['CCDs'].append({
                        'number'        : n_ccd,
                        'run_number'    : run,
                        'image_number'  : n_image,
                        'fits'          : meccd.fits.values.tolist()[0],
                        'fits_id'       : meccd.fits_id.values.tolist()[0],
                        'start'         : meccd.start.min().to_pydatetime(),
                        'end'           : meccd.end.max().to_pydatetime(),
                        'status'        : all(meccd.status.values),
                        'slow_control'  : { 'T': meccd['T'].mean(), 'T_std': meccd['T'].std(), 'P': meccd['P'].mean(), 'P_std': meccd['P'].std() },
                        'wr_darkcurrent': meccd.wr_darkcurrent.values[-1].split(";"),
                        'wr_overscan'   : meccd.wr_overscan.values[-1].split(";"),
                        ### adding extra information to the CCD, for the WEB TABLE page
                        'nskips'        : meccd.nskips.values.tolist()[0],
                        'npbin'         : meccd.npbin,
                        'nsbin'         : meccd.nsbin,
                        'texp'          : meccd.texp,
                        'tread'         : meccd.tread,
                        'nrows'         : meccd.nrows,
                        'ncols'         : meccd.ncols,
                        'MEs'           : []
                    })

                    #print("Filling ME per CCD .... ")
                    self.fill_ME('ccd',meccd,self.document[run]['Images'][-1]['CCDs'][-1]['MEs'])
                    #print("CCD filled")

                #print("Filling Image MES ...... ")
                self.fill_ME('image', meimage, self.document[run]['Images'][-1]['MEs'])
                #print("IMAGE filled")

            #print("Filling RUN MEs ....... ", merun.size)
            self.fill_ME('run', merun, self.document[run]['MEs'])
            #print("RUN FILLED")

    def fill_ME(self,objtype,df,MEList):

        #### APPEND ALL MEs to the list
        for mename in sorted(set(df.name.values)):
            me = df[df.name==mename]

            ### DEFINE COMMON 'UNIQUE' parameters for the ME, like name, title, xlabel and ylabel, 
            #       ... which is only one per ME
            data = {'name':mename, 
                    'title':me.title.values[0], 
                    'info': me.caption.values[0],
                    'xlabel': me.xlabel.values[0],
                    'ylabel': me.ylabel.values[0],
                    'status': all(me.status.values.tolist()),
                    'tests':  me.tests.values.tolist()[0]}
            
            ### DEPENDING OF THE TYPE OF ME [Run/Image/CCD] we show the y-values (yaxis, amplifier, goodnes, deltatime, ... ) in different statistics
            #         CCD: raw data: one value for each amplifier,
            #       Image: the mean value per image number, and then one per (image number, amplifier)
            #         Run: the mean value per run number, and then one per (run_number, amplifer)
            #
            # So each parameter AMP, GOODNES, ... can be used on the web as color code, marker code, or even as the y-axis values, for that reasson
            # this parameters needs to be reshaped to the same size as the x-values
            if objtype == 'ccd':
                # RESHAPING AMPLIFIER AND DELTATIME
                if type(me.x.values.tolist()[0])==list:
                    Nx = [len(x) for x in me.x.values.tolist()]
                    amp = np.repeat(me.amplifier.values.tolist(),Nx).tolist()
                    dt = [get_deltatime_days(t) for t in me['start']]
                    deltatime = np.repeat(dt,Nx)
                    ccd_vals = np.repeat(me.ccd.values.astype(str).tolist(),Nx).tolist()
                else:
                    Nx = None
                    amp = me.amplifier.values.tolist()
                    deltatime = [ get_deltatime_days(t) for t in me['start']]
                    ccd_vals = me.ccd.values.astype(str).tolist()
                
                #print(data['name'])
                #print(fits_vals)

                # SPECIAL CASE THE IMAGE, WHERE THE YVALUES IS A MATRIX, AND X IS THE IMAGE NUMBER (AS THERE IS ONLY ONE IMAGE PER INPUT FILE)
                # print(" read the ccd number ", me.ccd.values.tolist()[0])
                if mename.count('MECCD')==0:
                    # any ME that is not an image
                    data.update({
                        'ccd'       : ccd_vals,
                        'amplifier' : amp,
                        'deltatime' : deltatime,
                        'x'         : sum(me.x.values.tolist(),[]) if type(me.x.values.tolist()[0])==list else me.x.values.tolist(),
                        'y'         : sum(me.y.values.tolist(),[]) if type(me.y.values.tolist()[0])==list else me.y.values.tolist(),
                        'xerror'    : sum(me.xerror.values.tolist(),[]) if type(me.xerror.values.tolist()[0])==list else me.xerror.values.tolist(),
                        'yerror'    : sum(me.yerror.values.tolist(),[]) if type(me.yerror.values.tolist()[0])==list else me.yerror.values.tolist(),
                        'yfit'      : sum(me.yfit.values.tolist(),[]) if type(me.yfit.values.tolist()[0])==list else me.yfit.values.tolist(),
                        'goodness'  : np.round(sum(me.goodness.values.tolist(),[]),2).tolist() if type(me.goodness.values.tolist()[0])==list else np.round(me.goodness.values.tolist(),2).tolist(),
                        'status'    : all(me.status.values.tolist())
                    })
                    # add fits file information
                    data.update({
                        'fits'      : '{}, ID:{}'.format(me.fits.values.tolist()[0],me.image.values.tolist()[0])
                        })

                    # add extrainfo parameter, only for CCD MEs
                    #if hasattr(me,'extrainfo'):
                    #    data.update({
                    #        'extrainfo': sum(me.extrainfo.values.tolist(),[])
                    #    })
                else:
                    # to fill the me as an image
                    data.update({
                        'ccd'       : str(me.ccd.values.tolist()[0]),
                        'amplifier' : 'UL',
                        'deltatime' : min([ get_deltatime_days(t) for t in me['start']]),
                        'x'         : None,
                        'y'         : me.y.values.tolist()[0],
                        'xerror'    : None,
                        'yerror'    : None,
                        'yfit'      : None,
                        'goodness'  : None,
                        'status'    : True
                    })
            else:
                ### Exclude ME that are not Image or RUN 
                if mename in ['MEPCD','MEOVSPCD'] or mename.count('MECCD')>0: # images : 'MECCDEqualizedImage','MECCDMeanImage','MECCDStdImage','MECCDCalImage'
                    continue

                # ADD GENERAL VALUES FOR THE MONITOR ELEMENT
                ##############################################################################################################
                values = {}
                ### convert any parameter from the plot into a flat list
                for par in ['x','y','yfit','xerror','yerror','goodness']:
                    vals = getattr(me,par).values.tolist()
                    if type(vals[0])==list:
                        values[par] = sum( vals,[] )
                    else:
                        values[par] = vals
                
                # extend the other parameters to have the same length as the individual values of x
                if type(me.x.values[0])==list:
                    values['deltatime'] = sum( [[get_deltatime_days(t)]*len(me.x.values[i]) for i,t in enumerate(me['start'])  ], [])
                    values['amplifier'] = sum( [[a]*len(me.x.values[i]) for i,a in enumerate(me['amplifier'])  ], [])
                    values['ccd']       = sum( [[str(c)]*len(me.x.values[i]) for i,c in enumerate(me['ccd']) ], [])
                    values['fits']      = sum( [['{}, url: shifter_plots/{}/{}'.format(c[0],c[1],c[2])]*len(me.x.values[i]) for i,c in enumerate(zip(me['fits'],me['run'],me['image'])) ], [])  
                else:
                    # the parameter can be floats or list of floats
                    values['deltatime'] = [get_deltatime_days(t) for i,t in enumerate(me['start'])]
                    values['amplifier'] = [a for i,a in enumerate(me['amplifier'])]
                    values['ccd'] = [str(c) for i,c in enumerate(me['ccd'])]
                    values['fits'] = ['{}, url: shifter_plots/{}/{}'.format(c[0],c[1],c[2]) for c in zip(me['fits'],me['run'],me['image'])]

                
                ####   THE MONITOR ELEMENT WILL CONTAIN: THE MEAN VALUE FOR ALL ENTRIES, A MEAN VALUE FOR EACH CCD,AMPLIFIER
                ##############################################################################################################
                # >>>> ADDING THE MEAN: data frame to groupby 'x'
                # do not include the total mean value for the DCFit as the fit has not been done with the full set pixel charge distribution
                #if mename in ['MEFitDC']:
                #    ### initialize for the MEFitDC
                for p in ['deltatime','x','y','xerror','yerror','yfit','status','amplifier','goodness','ccd','fits']:
                    data[p] = []
                
                ## ADD DATA doing the mean AMP per AMP -- only for some ME
                # >>>> ADDING THE MEAN PER AMPLIFIER     
                # DIFFERENCIATE BETWEEN RUN AND IMAGE, FROM EVOLUTION PLOTS
                if mename in EVL_PLOTS and objtype == 'run':
                    # group by file name
                    dfampmean = pd.DataFrame.from_dict(values).groupby(['ccd','amplifier','fits','x'], as_index=False).mean()
                    data['fits'].append(dfampmean.fits.values.tolist())

                else:
                    # do not group, so, ME will not contain fits file information
                    dfampmean = pd.DataFrame.from_dict(values).groupby(['ccd','amplifier','x'], as_index=False).mean()
                    if len(set(me.image.values.tolist()))>1:
                        data['fits'].append(['{}...{}'.format(me['image'].min(),me['image'].max())]*len(dfampmean.x.values.tolist()))
                    else:
                        data['fits'].append(['{}'.format(me['image'].values[0])]*len(dfampmean.x.values.tolist()))
                
                data['deltatime'].append( dfampmean.deltatime .values.tolist() )
                data['x'].append(dfampmean.x.values.tolist())
                data['y'].append(dfampmean.y.values.tolist())
                data['yfit'].append(dfampmean.yfit.values.tolist())
                data['xerror'].append((dfampmean.xerror.values).tolist())
                data['yerror'].append((dfampmean.yerror.values).tolist())
                data['status'].append(all(me.status.values.tolist()))
                data['amplifier'].append(dfampmean.amplifier.values.tolist())
                data['goodness'].append(dfampmean.goodness.round(2).values.tolist())
                data['ccd'].append(dfampmean.ccd.values.tolist())


                #### flattening the list of list to be a list, important step to get a flat list 
                # as plotly does not accept list of list of different sizes
                for p in ['deltatime','x','y','xerror','yerror','yfit','amplifier','goodness','ccd','fits']:
                    data[p] = sum(data[p],[])
                
            MEList.append(data)
        return

    def close_DB_document(self,run,me_outdir):
        if not hasattr(self,'document'):
            return

        print(" DQM INFO: save mongoDB documents as feather files")
        for run in self.document.keys():
            fname = '{}/mongoDB_document_run{:03}.pkl'.format(me_outdir,int(run))
            print("     ..document (as pkl): ", fname)
            with open(fname,'wb') as fpkl:
                pickle.dump(self.document[run], fpkl)
        return


    def save_MEs_as_dataframe(self,me_outdir,infile):
        """Save all ME of an image as an npz file
        """

        if not hasattr(self,'mongodb'):
            return

        for run in self.mongodb.keys():
            # me_outdir, run
            out_file_name = '{}/../mongoDB_document_{}_run{:03}.pkl'.format(me_outdir,infile,int(run))
            #### storing as pickle object
            print("DQM INFO: save all MEs as {}".format(infile,out_file_name))

            file_object = open(out_file_name, 'wb')
            pickle.dump(self.mongodb[run],  file_object)
            file_object.close()

        return out_file_name

#################################################################################################
#       CLASS FOR THE RUN REFERENCE
#
#################################################################################################
class ModelRef(object):
    def __init__(self,me_output):
        for key_attr in me_output:
            setattr(self,key_attr, me_output[key_attr])

class MERef(object):
    def __init__(self,me_summary_ref):
        for key_class in me_summary_ref.keys():
            if not key_class in ['run', 'me_dc_gain','me_dc_lambda','me_dc_sigma','me_dc_mu0']:
                setattr(self,key_class,ModelRef(me_summary_ref[key_class]))

def fget(_x,func,mask=None,name=None):
    """Function to compute statistical values for the ME of the Run and/or images
    """
    if mask is not None:
        # do the mean by run
        x,y,dtime = np.array(_x[0],dtype=object)[mask],np.array(_x[1],dtype=object)[mask],np.array(_x[2],dtype=object)[mask]
        xflat,yflat,dtflat = [],[],[]
        ### flattening x and y (flatten() can not be used, as different size images can be
        #       reprocessed at the same run
        for xi,yi,dti in zip(x,y,dtime):
            try:
                xflat.extend(xi)
            except TypeError:
                xflat.append(xi)
            try:
                yflat.extend(yi)
            except TypeError:
                yflat.append(yi)
            try:
                dtflat.extend(dti)
            except TypeError:
                dtflat.append(dti)

        if not name in ['MEPCD','MEOVSPCD']:
            df = pd.DataFrame.from_dict({'x':xflat,'y':yflat,'dt':dtflat})
            dfmean = df.groupby('x').mean()
            dfstd  = df.groupby('x').std()
            dfmin  = df.groupby('x').min()

            return dfmean.index.values, dfmean.y.values, dfstd.y.values, dfmin.dt.values
        else:
            return xflat,yflat,[0.0]*len(xflat),[0.0]*len(xflat)


    if func is list:
        try:
            return sum(_x,[])
        except TypeError:
            return _x
    try:
        if len(_x[0])>0:
            return func(_x,axis=0)
        else:
            return sum(_x,[])
    except TypeError:
        return _x

