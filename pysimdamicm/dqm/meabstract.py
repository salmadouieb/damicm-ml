""":mod:`pysimdamicm.dqm.me`

   Processes defining the Monitor Elements (MEs)
   to control the quality of taken skipper-CCD images.

   It only contains those MEs related to single images.

   No quality test have been added so far.

 
.. moduleauthor:: Nuria Castello-Mor
"""

from abc import ABCMeta, abstractmethod
from cmath import isnan

from pysimdamicm.dqm.qtest_library import * 
from pysimdamicm.utils.units import Units
u=Units()
from pysimdamicm.utils.plotLib import GetROOTPalette
from pysimdamicm.utils.libplot4ana import rc,plt,update_mpl_style
update_mpl_style()


import ROOT
import pandas as pd
import numpy as np
from scipy import stats
from array import array
from matplotlib import pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

###############################################################################################
#####       ABSTRACT CLASS TO DEFINE THE DIFFERENT ME
#####
###############################################################################################

class MEAbstract(metaclass=ABCMeta): 
    """
    """

    def __init__(self):
        """
        """
        self.failed = False

        # unset whanever the ME is independent of the amplifier (such us MECCDImage)
        self._per_amp = True

        ### VARIABLES TO DEFINE IF ONE MONITOR IS PER CCD/IMAGE OR PER RUN (A SET OF IMAGES)
        self.is_MECCD   = True
        self.is_MEImage = True
        self.is_MERun   = True

        ### RUN and IMAGE NUMBERS ARE MANDATORY
        self.run = []
        self.n_image = []

        ### ME is a representation of a plot where
        #       x can be column, rows, image number or run number
        #       while y is the parameter to be monitorize (i.e. ME)
        self.x = []
        self.y = []
        # errors may not exist for all ME
        self.xerr = []
        self.yerr = []
        # A CCD can be readout by several amplifiers, a silicon portion is well defined by the tuple
        #       (CCD,AMP)
        self.ccd = []
        self.amplifier = []
        self.amplifier_list = []
        # time to create evolution plots in the DQM web
        self.start = []
        
        # fitted model to data
        self.y_fit = []
        self.chi2  = np.nan

        ###### REFERENCE MODEL and Qtests
        self.qtest_results  = [{'name':'QTnone','result':True, 'x':None, 'y':None,'mean_ref':0,'n_std':0, 'std_ref':0, 'run_ref':0}]
        self.qt_valid_range = [None,None]

        ###### ME OFFLINE MODE
        #       In the offline version a PDF report is generated with all the runtime plots. This
        #       variable is a list with the full-path-file names of these plots
        self.fig_png_file = []

        ##### DEFAULT STYLES
        self.set_plotting_style()
        
        ##### DEFAULT JSON CONFIGURATION OPTIONS
        self.verbose = False
        self.saveas  = True
        # set if the plots should include errors on the y-axis
        self.plot_with_yerrors = False
        self.__units__ = {'verbose':1, 'saveas':1, 'plot_with_yerrors':1, 'qt_valid_range':1}


    def set_plotting_style(self):
        """Set default MATPLOTLIB.PYPLOT style
        """

        self.alpha = 0.6
        self.marker = {'U': 'o', 'L': '^'}
        self.linestyle ='None'
        self.markerfacecolor = {'U':'#92b9f0', 'L': '#d1bbb6'}
        self.markeredgecolor = {'U':'#0d3b7a', 'L': '#4a4544'}
        self.color = '#043bd1'
        self.color_fit = '#020fc9'
        self.color_percentile = '#45afed'
        self.lw_fit = 3
        self.color_ref = '#e61239'
        self.linewidth = 1
        self.capsize = 3
        
        self.xscale ='linear'
        self.yscale = 'linear'

        ### size for scatter plots
        self.markersize = 5
        self._set_ylim = True

        ### figure size
        self.fig_size = (10,5)
        self.dpi = 200
        self.loc_legend = 0
        self.palette = 'YlOrRd_r' #GetROOTPalette(83,dtype='palette')


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
                    setattr(self,keyword,val)
                else:
                    setattr(self,keyword,val*self.__units__[keyword])
            else:
                raise(AttributeError('"{}" invalid attribute for class {}'.format(keyword,self.__class__.__name__)))

    def info(self):
        """Show all attributes of the process
        """
        print("<{}> with sequence id {}.\n List of public data members: ".format(self.__name__,self.__sequence_id__)) 
        for attr in sorted(self.__units__.keys()):
            print("\t * {} = {} ".format(attr, getattr(self,attr)))

    
    def add_data(self,x,y,amp,**kwargs):
        """Append MANDATORY attributes to the ME

        The mandatory elements are x,y,ccd, amplifier, run, image number and datetime start, being xerr and yerr
        optional parameters

        """
        # mandatory
        mkeys = ['x','y','amplifier']
        for key,param in zip(mkeys,[x,y,amp]):
            if param is not None:
                p = getattr(self,key)
                p.append( param )

        # optionals
        for key in ['xerr','yerr','y_fit']:
            if key in kwargs:
                p = getattr(self,key)
                p.append( kwargs[key] )
        
        # add automatic linear fit if x is a list, and y_fit is not in kwargs
        # some monitor elements like histograms, has no values on the y-axis, skip this ones
        if self.__name__ in ['MEFitDC','MEPCD','MEOVSPCD','MECCDEqualizedImage','MECCDMeanImage','MECCDStdImage','MEMeanPixelChargePerRow','MECCDCalImage']:
            return
        if not 'y_fit' in kwargs and type(x)==list:
            # append fit with the method do_linear_fit
            yerr = None if not 'yerr' in kwargs else kwargs['yerr']
            self.do_linear_fit(x,y,yerr=yerr)


    def add_image(self,x,y,append=False,**kwargs):
        """Append MANDATORY attributes to the ME

        The mandatory elements are x,y 

        """
        # mandatory
        mkeys = ['x','y']
        for key,param in zip(mkeys,[x,y]):
            if param is not None:
                p = getattr(self,key)
                if append:
                    p.append(param)
                else:
                    setattr(self,key,[param])

    def get_data_as_pandas(self,rawdata):
        """Create a DataFrame with the full information needed to fill the mongoDB
        """
        
        if hasattr(rawdata,'QT_{}_{}'.format(rawdata.ccd, rawdata.execute_process_in_amp)):
            qtest_results = getattr(rawdata,'QT_{}_{}'.format(rawdata.ccd, rawdata.execute_process_in_amp))
            # delete before running another ME
            delattr(rawdata,'QT_{}_{}'.format(rawdata.ccd, rawdata.execute_process_in_amp))
        else:
            qtest_results = [{'name':'QTnone','result':True}]

        status = np.all([qt['result'] for qt in qtest_results])

        if type(self.chi2)==list:
            goodness = self.chi2
            status = False if np.max(np.abs(goodness)) > 10. else status
        else:
            chi2 = 1 if np.isnan(self.chi2) else self.chi2
            if type(self.x[-1])==list:
                goodness = [chi2]*len(self.x[-1])
                status = False if any(goodness)>10 else status
            else:
                goodness = chi2
                status = False if goodness > 10 else status
        

        nskips = int(rawdata.image_header['NSKIPS']) if 'NSKIPS' in rawdata.image_header else rawdata.nskips
        if nskips==1 and rawdata._NDCM>1:
            nskips = rawdata._NDCM

        data = {
            'name':     self.__name__,
            #### information about data
            'npbin':    rawdata.npbin,
            'nsbin':    rawdata.nsbin,
            'texp':     rawdata.exposure_time,
            'tread':    rawdata.read_time,
            'nrows':    rawdata.nrows,
            'ncols':    rawdata.ncols,
            'nskips':   nskips,
            #### metadata from rawdata
            'run':      rawdata.n_run,
            'image':    rawdata.n_image,
            'fits':     rawdata._file_name,
            'fits_id':  rawdata.n_image_from_file,
            #'fimage':   '{}:{}'.format(rawdata.ccd,rawdata.n_image_from_file),
            'ccd':      rawdata.ccd,
            'tag':      rawdata.run_tag,
            'amplifier':rawdata.execute_process_in_amp+str(rawdata.ampdir),
            'start':    rawdata.start,
            'end':      rawdata.end,
            'T':        rawdata.T[0],
            'T_err':    rawdata.T[1],
            'P':        rawdata.P[0],
            'P_err':    rawdata.P[1],
            #### status of the test for this (CCD,amp) or CCD
            'status':   status,
            'tests':    qtest_results,
            #### for ME plots
            'title':    self.title,
            'caption':  self.caption,
            'xlabel':   self.xlabel,
            'ylabel':   self._ylabel,
            'x':        self.x[-1],
            'xerror':   self.get_errors(self.x[-1],self.xerr),
            'y':        self.y[-1],
            'yerror':   self.get_errors(self.x[-1],self.yerr),
            'yfit':     self.y_fit[-1] if len(self.x)==len(self.y_fit) else self.get_errors(self.x[-1],self.y_fit),
            'goodness': goodness,
            'wr_darkcurrent': ";".join(list(set(rawdata.wr_darkcurrent))).strip(),
            'wr_overscan':    ";".join(rawdata.wr_overscan)
            }

        # Some ME will have extra information that will then be used in the web in general terms. 
        # As this infor is not common to all of them, we define a general parameter to be stored
        # This will be display by using the hover mechanism from the plotly.graph.objects
        if hasattr(self,'extrainfo'):
            data.update({'extrainfo': self.extrainfo})

        data = pd.DataFrame.from_dict(data,orient='index').transpose()

        return data
    
    def get_errors(self,x,xerr):

        if len(xerr)>0:
            return xerr[-1]
        else:
            if type(x)==list:
                return [0]*len(x)
            else:
                return 0


    @abstractmethod
    def execute_process(self,rawdata,**kwargs):
        """
        """
        raise NotImplementedError()


    def execute_qtest(self,y,ccd,amp,name,**kwargs):
        """ Execute a simple quality test: 
                any value (y-axis) must be within a given range of values
            
        Returns
        -------
            True/False if the test passed sucssesfully or not 
            (i.e. False mans that the test have been failed)

        
        Object to return should be a list of 
        [{'name':'QTnone','result':True, 'x':None, 'y':None,'mean_ref':0,'n_std':0, 'std_ref':0, 'run_ref':0}]

        For now we just return a simplier qtest object
        [{'name':'QTnone','result':True}]

        """
        print(name, self.qt_valid_range)

        if type(self.qt_valid_range)==list and self.qt_valid_range[0] is None:
            return [{ 'name':'QTnone', 'result':True }]
        
        if type(self.qt_valid_range)==dict:
            try:
                ymin,ymax = self.qt_valid_range['{},{}'.format(ccd,amp)]
            except KeyError:
                return [{ 'name':'QTnone', 'result':True }]
        else:
            ymin,ymax = self.qt_valid_range     

        if type(y)==list:
            y = np.array(y)
        else:
            y = np.array([y])

        if name=='QT:range':
            result = not (sum(np.array(y)<ymin)>0 or sum(np.array(y)>ymax)>0)
        elif name=='QT:chi2':
            result = not (kwargs['chi2']>5 or kwargs['chi2']<0)
        elif name=='QT:chi2,range':
            result = not( sum(np.array(y)<ymin)>0 or sum(np.array(y)>ymax)>0 or kwargs['chi2']>5 or kwargs['chi2']<0 )
        
        return [{'name':name,'result':result}]

        #self.qtest_results = []
        #for qtest in self.__qtest_list__:
        #    qcut = self.__qtest_cut__[qtest.__name__]
        #    if len(self.x) > 0:
        #        if is_image:
        #            x,y = self.x[-1],self.y[-1]
        #        else:
        #            x,y = self.x,self.y
        #        if qtest == qtest_chi2FromFit: 
        #            qt = qtest(x,self.chi2,None,amp,None, qcut, is_image)
        #        else:
        #            qt = qtest(x,y,self.me_ref, amp, n_std, qcut, is_image)
        #        self.qtest_results.append(qt)
    

    ############################################################################################################# PLOT FUNCTION
    def execute_plot(self,ptype='scatter',**popts):
        """Create an scatter plot with or without errors according to xerr/yerr
        """
        # Before doing nothing be sure the ME has something to plot
        if not len(self.x)>0:
            return
        # XXX Is still needed????matplotlib associate a number ID to the figure, to keep track of this and do not replicate
        # them, the singleton Units will be used to keep track of this id
        # u.n_figure+=1

        # create the canvas/frame for the figure
        self._figure = plt.Figure(figsize=self.fig_size,constrained_layout=True)
        ax = self._figure.add_subplot(111)
        if 'title' in popts:
            ax.set_title(popts['title'])
        else:
            ax.set_title(self.title+"\n [class {}]".format(self.__name__))

        
        # ADD PLOT
        if ptype=='scatter':
            self.add_scatter_plot(ax)
        elif ptype=='bar':
            self.add_bar_plot(ax)
        elif ptype=='hist':
            self.add_histogram_plot(ax,**popts)
        elif ptype=='imshow':
            self.add_imshow_plot(ax,**popts)
        
        # SET LABELS
        self.set_figure_properties(ax,**popts)

        # SAVE FIGURE
        canvas = FigureCanvasAgg(self._figure)
        outname = self.get_output_fig_file_name(**popts) 
        canvas.print_figure(outname, dpi=self.dpi)
        
    def get_output_fig_file_name(self,**popts):
        # output figure name
        suffix = self.__name__
        if 'suffix' in popts:
            suffix = "{}_{}".format(self.__name__,popts['suffix'])
        if hasattr(self,'output'):
            outname = "{}{:03}_{}.png".format(self.output,self.__sequence_id__,suffix)
        else:
            outname = "{:03}_{}.png".format(self.__sequence_id__,suffix)

        self.fig_png_file.append(outname)
        return outname

    def set_figure_properties(self,ax,**popts):
        # added as data members
        for attr in ['xlabel','ylabel','ylim','xlim','yscale','xscale']:
            if hasattr(self,attr):
                _setfunc = getattr(ax,"set_"+attr)
                _setfunc(getattr(self,attr))

        # actualize figure settings: labels, scale, limits, ... 
        for attr in popts:
            if attr in ['xlabel','ylabel','ylim','xlim','yscale','xscale']:
                _setfunc = getattr(ax,"set_"+attr)
                _setfunc(popts[attr])

        # xticks
        if 'xticks' in popts:
            _ = ax.set_xticks(popts['xticks'],minor=True)
        
        return

    def add_scatter_plot(self,ax):
        """
        """
        add_legend = True
        ampset = set(self.amplifier)

        # Use the error on the y-axis only if the user booked (or set by default)
        if self.plot_with_yerrors:
            yerr = self.yerr if len(self.yerr)>0 else [None]*len(self.x)
        else:
            yerr = [None]*len(self.x)

        # ZOOM in when outliers
        ax2 = None
        try:
            yflat = np.concatenate(self.y).flat
        except ValueError:
            yflat = self.y

        if not self.__name__ in ['MEPCD','MEOVSPCD']:
            ### only when the plot is not histogram type
            p05,p95 = np.percentile(yflat,[1,99])
            if p05/np.min(yflat)>5.0 or np.max(yflat)/p95>5.0:
                ax2 = self._figure.add_axes([0.12,0.5,0.65,0.35])
                self.loc_legend=1
       
        # plot ALL ME by (CCD,AMP)
        _amp = []
        kind = 0
        for x,y,amp,ey in zip(self.x,self.y,self.amplifier,yerr):
            if not amp in _amp and kind<4:
                _amp.append(amp)
                label = amp
            else:
                label = None

            ax.errorbar( x, y, yerr=ey,
                    # colors/shape according to the amplifier
                    marker=self.marker[amp],mfc=self.markerfacecolor[amp],mec=self.markeredgecolor[amp],
                    # rest equal for all datasets
                    ms=self.markersize,capsize=self.capsize,linestyle='None',alpha=self.alpha,
                    label=label )
            if ax2 is not None:
                ax2.errorbar( np.array(x)[np.logical_and(y>p05,y<p95)], np.array(y)[np.logical_and(y>p05,y<p95)],
                        marker=self.marker[amp],mfc=self.markerfacecolor[amp],mec=self.markeredgecolor[amp],
                        ms=self.markersize,capsize=self.capsize,linestyle='None',alpha=self.alpha)
            
            kind+=1

        if len(self.y):
            ax.legend(loc=self.loc_legend, markerscale=2, fancybox=True, framealpha=0.2)

        return

    def add_bar_plot(self,ax):
        """
        """
        _amp = [] 
        for x,y,amp in zip(self.x,self.y,self.amplifier):
            if not amp in _amp:
                _amp.append(amp)
                label = amp
            else:
                label = None

            ax.bar(x,y, edgecolor=self.markeredgecolor[amp], color=self.markerfacecolor[amp],
                    label=label )

        return
    
    def add_histogram_plot(self,ax,**kwopts):
        """
        """
        # XXX what if the RUN is different ??
        _amp = []
        if 'nbins' in kwopts:
            nbins = kwopts['nbins']
        else:
            nbins = int((np.max(np.array(self.y)) - np.min(np.array(self.y)))*3)

        if 'fill' in kwopts:
            fill = kwopts['fill']
        else:
            fill = True

        if 'histtype' in kwopts:
            histtype = kwopts['histtype']
        else:
            histtype = 'step'
        
        labels = ['{},{},{}'.format(x,y,z) for x,y,z in zip(self.n_image,self.ccd,self.amplifier)]
        # plot an stacked histogram: allow to see individual contribution, and the final plot
        ax.hist( self.x, nbins, stacked=kwopts['stacked'], fill=fill, label=labels,
                histtype=histtype )
        ax.set_yscale('log')
        
        ### ADD ADUe2 used
        ax.text(0.60,0.90, 'calibration: {} ADU/e-'.format(u.ADC2e), color='green', fontsize=14,
                transform=ax.transAxes,verticalalignment='bottom', horizontalalignment='right')

        ax.legend(loc=self.loc_legend, fancybox=True, framealpha=0.2)
        return

    def add_imshow_plot(self,ax,**kwopts):
        """Add image plot
        """
        # plot the image with the default colorcode palete
        if 'vrange' in kwopts:
            vrange = kwopts['vrange']
            fimg = ax.imshow(self.y[-1], aspect='auto', origin='lower', vmin=vrange[0], vmax=vrange[1], cmap=self.palette)
        else:
            fimg = ax.imshow(self.y[-1], aspect='auto', origin='lower', cmap=self.palette)

        # add the color bar, and reduce the space in btw
        bimg = self._figure.colorbar(fimg, pad=0.015)

        # legend for the color bar
        if 'cbar' in kwopts:
            bimg.set_label(kwopts['cbar'],rotation=270,labelpad=19)

        # add vertical lines to highlight the overscan region
        if 'vlines' in kwopts:
            for v in kwopts['vlines']:
                ax.axvline(x=v, color='#f20a19', ls='--', alpha=0.8 )

        return



    ######################################################################################################## UTIL FUNCTIONS
    def print_warning(self,msm):
        print("\x1b[93m    WARNING. {} \x1b[m".format(msm))

    def print_alert(self,msm):
        print("\x1b[93m    ALERT. {} \x1b[m".format(msm))

    def median_and_mad(self,image,axis=None,keepdims=1):
        """Median absolute debiation for masked arrays
        """

        median = np.ma.median(image,axis=axis, keepdims=keepdims)
        mad    = np.ma.median( np.ma.abs(image - median), axis=axis, keepdims=keepdims)

        return median,mad

    def do_linear_fit(self,x,y,yerr=None,zscore_max=3,**kwargs):
        """Before fitting remove outliers based on the Z-score, which is defined as

            z = (x -mean)/std
        """
        
        _xinput = x
        _yinput = y
        
        # if the zscore is too restrictive we can end up with no data to do the fit
        # ignore filter outliers in that case
        def ignore_outlier(x,y,yerr,zscore_max):
            if zscore_max is not None:
                mask_zscore = np.abs(stats.zscore(y))<zscore_max
                if len(mask_zscore)<len(x):
                    x = np.array(x)[mask_zscore]
                    y = np.array(y)[mask_zscore]
                    if yerr is not None:
                        yerr = np.array(yerr)[mask_zscore]
            return x,y,yerr
         
        x,y,yerr = ignore_outlier(x,y,yerr,zscore_max)

        if yerr is None:
            sigma = np.array([0]*len(x))
        else:
            sigma = np.array(yerr)
        
        x = np.array(x)
        y = np.array(y)
        xerr = np.array([0]*len(x))
        
        tg = ROOT.TGraphErrors(len(x), array('f',x.tolist()), array('f', y.tolist()), array('f',xerr.tolist()), array('f', sigma.tolist()) )
        
        fitfunc = ROOT.TF1('fitfunc','pol1')
        # fitting to unbinned data
        res = tg.Fit(fitfunc,'SEQ')
        res = res.Get()
        values = [res.GetParams()[1], res.GetParams()[0]]
        values_err = [res.GetErrors()[1], res.GetErrors()[0]]
        
        ### update model and errors
        self.fit_model = np.poly1d(values)
        self.fit_model_errors = values_err

        # update chi2
        self.chi2 = fitfunc.GetChisquare()/fitfunc.GetNDF()
        # update yfit values -- using the full input set of x's values
        self.y_fit.append(self.fit_model(_xinput).tolist())

        return
    


