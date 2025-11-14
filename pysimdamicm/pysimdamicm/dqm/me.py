""":mod:`pysimdamicm.dqm.me`

   Processes defining the Monitor Elements (MEs)
      to control the quality of taken skipper-CCD images.

.. moduleauthor:: Nuria Castello-Mor
"""


from pysimdamicm.dqm.meabstract import MEAbstract
from pysimdamicm.dqm.qtest_library import * #qtest_tt,qtest_residuals_per_image
from pysimdamicm.processes.skipper_analysis import FitDarkCurrentProcess, FitDarkCurrentPerRow, FitDarkCurrentPerCol, CalibrationProcess
from pysimdamicm.processes.skipper_comissioning import ChargeLossSkewnessProcess, RNvsNskipsPlot
from pysimdamicm.utils.units import Units
u=Units()

from datetime import datetime
import time
import numpy as np
import ROOT
from array import array
from astropy.io import fits

#from memory_profiler import profile

from pysimdamicm.utils.root_plot_styles import damicmStyle
style = damicmStyle()
style.cd()
ROOT.gROOT.ForceStyle()

ROOT.gErrorIgnoreLevel = ROOT.kWarning

DAY0 = datetime.strptime('2022-01-01 00:00:00','%Y-%m-%d %H:%M:%S')
datetime_to_days = lambda t: (t - DAY0).total_seconds() /60./60./24.

#####################################################################################################
#####       ME CLASSes                                                                              #
#####                                                                                               #
#####################################################################################################
##########################################################################################################################################
################ MONITOR ELEMENTS to control THE NDCM/skips
##########################################################################################################################################
class MEQmeanDiffSingleSkips(MEAbstract):
    """Mean charge of the pixel charge distribution of the difference between single skip images
    """
    __name__ = "MEQmeanDiffSingleSkips"
    __sequence_id__ = -1

    def __init__(self):
        super().__init__()
        
        # this ME will be used to find images out of the normality, i.e. is a RUN ME, but not an IMAGE ME
        self.is_MECCD = False

        self.title   = "Mean charge of the pixel charge difference distribution"
        self.caption_fmt = "It is build from the single-skip images (SSI). For each NDCM, a SSI is build with the n-th NDCM measurement of each pixel. One of these images is selected as the SSI image of reference (here {}), and pixel by pixel the difference between the SSI of reference and the n-th SSI is done. These images have as a pixel value the difference between two different skips measurements. The averaged value of all pixel values is display as a function of the NDCM ordinal. To speed up the process is done in steps of {} NDCMs."
        self.xlabel = "NDCM ordinal"
        self.ylabel = "$\\textit{mean}(\langle q - q_{ref}\\rangle) \; [ADU]$"
        self._ylabel = "mean(q-qref) [ADU]"
        
        self.skip_reference = 200
        self.skip_step = 5

        self.__units__.update({'skip_reference':1,'skip_step':1})
        
    @property
    def caption(self):
        return self.caption_fmt.format(self.skip_reference,self.skip_step)

    def execute_process(self,rawdata,**kwargs):
        
        if not len(rawdata.image.shape)==3:
            raise AttributeError("<MEQmeanDiffSingleSkips>: the input image should be an image with more than 1 skip measurements.")
        
        # set output directory 
        self.output = rawdata.me_output

        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        image_amp = rawdata.get_image_by_attribute_name("mean_compressed").copy()
        
        qmean = np.median(image_amp.data)
        
        ydiff  = []
        xskips = []
        for skip in range(0,rawdata.nskips,self.skip_step):
            ydiff.append( np.mean(np.ma.array(rawdata.image[:,:,self.skip_reference],mask=rawdata.mask_image_active_region) - np.ma.array(rawdata.image[:,:,skip],mask=rawdata.mask_image_active_region))/qmean )
            xskips.append(skip)

        # ADD DATA
        self.add_data( xskips, ydiff, rawdata.execute_process_in_amp )

        return

class MESTDSkipsPerCol(MEAbstract):
    """Mean charge of the pixel charge distribution of the difference between single skip images
    """
    __name__ = "MESTDSkipsPerCol"
    __sequence_id__ = 2

    def __init__(self):
        super().__init__()
        
        # this ME will be used to find images out of the normality, i.e. is a RUN ME, but not an IMAGE ME
        self.is_MECCD = False
        # only once per CCD
        self._per_amp = True

        self.n_std = 8

        self.title   = "Mean charge of the pixel charge difference distribution"
        self.caption_fmt = "This ME is constructed from the STANDARD DEVIATION IMAGE, where the pixel value is the STD of the acquired NDCM values. The bar means the number of pixels per column (x-axis) with an standard deviation larger than {} times the averaged STD of the STD image. Note that 'bars' on the plot are stacked per amplifier, to properly compare both amplifiers. If a column has a contribution from both amplifiers it may be due to cross-talk (or similar), so it is important to CHECK the STD IMAGE that can be found at the `Raw data` page."
        self.xlabel = "relative position to the amplifier"
        self.ylabel = "STD's NDCM outliers"
        self._ylabel = "STD's NDCM outliers"
        
        self.__units__.update({'n_std':1})
        
    @property
    def caption(self):
        return self.caption_fmt.format(self.n_std)

    def execute_process(self,rawdata,**kwargs):
        
        if not len(rawdata.image.shape)==3:
            raise AttributeError("<MESTDSkipsPerCol>: the input image should be an image with more than 1 skip measurements.")        
        if not hasattr(rawdata,'image_std_compressed'):
            raise AttributeError("<MESTDSkipsPerCol>: the STD must be booked, add 'std' to CompressSkipperProcess.func_to_compress in the config JSON file.")

        # set output directory 
        self.output = rawdata.me_output

        # get the STD image for the given ampl
        image_amp = rawdata.get_image_by_attribute_name("std_compressed")
        image_amp = np.ma.array(image_amp, mask=rawdata.mask_image_active_region)

        mean,std = np.ma.mean(image_amp),np.ma.std(image_amp)
        #p95      = np.ma.percentile(image_amp,95)
        Nstd_outliers = (image_amp > mean+self.n_std*std).sum(axis=0).compressed()

        xcols = np.linspace(1,len(Nstd_outliers),len(Nstd_outliers)).astype(int)

        # flip x-cols values for amplifier U, as we want to find companions with the same distance to the amplifier
        if rawdata.execute_process_in_amp=='U':
            Nstd_outliers = np.flip(Nstd_outliers)

        # ADD DATA
        self.add_data( xcols.tolist(), Nstd_outliers.tolist(), rawdata.execute_process_in_amp )
        
        return


class MEReadoutnoiseNDCM(MEAbstract):
    """Mean charge of the pixel charge distribution of the difference between single skip images
    """
    __name__ = "MEReadoutnoiseNDCM"
    __sequence_id__ = 3

    def __init__(self):
        super().__init__()
        
        # this ME will be used to find images out of the normality, i.e. is a RUN ME, but not an IMAGE ME
        self.is_MECCD = False

        self.title   = "Measured readout as a function of the number of skips"
        self.caption_fmt = "Measured readout noise as a function of the number of non-destructive readout samples per pixel for the skipper-CCD. Points are the standard deviation of the empty pixels distribution as a function of the number of averaged samples. The line is the theoretical expectation assuming independent, uncorrelated samples. The first {} measurements have been excluded."
        self.xlabel  = "Samples per pixel"
        self.ylabel  = "$\\textit{Readout noise} \;[e^{-}\\textit{rms/pix}]$"
        self._ylabel = "Readout Noise [ADU]"

        self.id_skip_start = 1
        self.n_skips_per_block = -50
        #self.yscale = 'log'


        self.__units__.update({'n_skips_per_block':1,'id_skip_start':1})
    
    @property
    def caption(self):
        return self.caption_fmt.format(self.id_skip_start)

    def execute_process(self,rawdata,**kwargs):
        
        if not len(rawdata.image.shape)==3:
            raise AttributeError("<MEReadoutnoiseNDCM>: the input image should be an image with more than 1 skip measurements.")

        # set output directory 
        self.output = rawdata.me_output

        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        image_amp = rawdata.get_image_by_attribute_name("mean_compressed").copy()
        
        pnoise = RNvsNskipsPlot()
        pnoise.n_skips_per_block = self.n_skips_per_block
        pnoise.id_skip_start     = self.id_skip_start
        pnoise.__silent__ = True        
        pnoise.execute_process(rawdata)

        # ADD DATA
        x = pnoise.results['NDCM']
        y = pnoise.results['noise']
        sigma_last = y[-1] * np.sqrt(x[-1])
        y_fit = [sigma_last/np.sqrt(xi) for xi in x]

        self.add_data( x, y, rawdata.execute_process_in_amp , y_fit=y_fit)
        
        # execute plot image by image, and skip the plot for run
        infile = "{}{}".format(rawdata.me_img,rawdata._file_name.split(".")[0])
        _ = self.execute_me_plot_per_file(rawdata.execute_process_in_amp,infile)

        return

    def execute_me_plot_per_file(self,amp,infile,**kwargs):
        """
        """
        ROOT.gROOT.SetBatch(1)

        c = ROOT.TCanvas()
        leg = ROOT.TLegend(0.85,0.85,0.9,0.9)
        
        x = array('d',self.x[-1])
        
        tg = ROOT.TGraph(len(x), x, array('d',self.y[-1]))
        tg.SetMarkerStyle(20)
        tg.SetMarkerSize(0.9)
        tg.SetMarkerColor(ROOT.kBlack)
        tg.GetHistogram().GetXaxis().SetTitle("Samples per pixel")
        tg.GetHistogram().GetYaxis().SetTitle("Readout noise (e^{-} rms/pix)")
        tg.GetHistogram().SetTitle(infile.split("/")[-1])
        leg.AddEntry(tg," amplifier {}".format(amp))
        tg.SetTitle(infile.split("/")[-1])
        tg.Draw("A P")

        # expected: 1 over square N
        #xm = np.arange(min(x),max(x),1)
        tgm = ROOT.TGraph(len(x), x, array('d',self.y_fit[-1]))
        tgm.SetLineColor(ROOT.kRed)
        tgm.SetLineStyle(9)
        tgm.SetLineWidth(2)
        leg.AddEntry(tgm,"1/\sqrt{N}","L")
        tgm.Draw("L SAME")
        
        leg.Draw("SAME")
        c.SetLogy()
        c.SetLogx()
        c.Update()
        c.Draw()

        # output file name
        outfile = "{}_amp{}_{}.png".format(infile,amp,self.__name__)
        # save figure
        c.SaveAs(outfile)
        # add output file png to the list of files to be included in the PDF report
        self.fig_png_file.append(outfile)

        return outfile.split('/')[-1]
    
    def execute_plot(self):
        """
        """
        # nothing to do, all plots have been done previously, one by one with the single images

        return




#####################################################################################################
#                                                                                                   #
#       DARK CURRENT AS THE MEDIAN OF AN IMAGE AFTER MASK PIXELS ABOVE THRESHOLD                    #
#                                                                                                   #
#####################################################################################################
class MEMeanPixelChargePerRow(MEAbstract):
    """ Selects those pixels of the pedestal subtracted image active area that satisfies :math:`q_{i,j}< q^{th}_{i}`. Where
    :math:`i,j` are the row and column number and :math:`q^{th}_{i} = median(q_{j})_{i}+n_{mad}*MAD(q_{j})_{i}`.

    The resulting median of the selected pixels is plotted vs the row number. This is done for every image in the RUN.

    A linear fit of the median per row of all images is performed. The slope and intersection values allow us to
    see the correlation between rows and charge.

    Attributes
    ----------
    n_mad : int
        It's the number of MAD the threshold is set above the median pixel value of the active
        area.
    """

    __name__ = "MEMeanPixelChargePerRow"
    __sequence_id__ = 4

    
    def __init__(self):
        super().__init__()
        
        # this ME acts on the pedestal subtracted image
        self.image = 'mean_compressed_pedestal_subtracted'

        self.title   = "Mean value of the pixel charge in a row (only pixels in the active region)."
        self.caption_fmt = "Average value of the pixel charge in a row, considering only values below {}xMAD. This ME should spot any issues in the level of background when compared with the reference (hot rows or regions). Note that only pixels in the sensitive region are included. Calibration used: {} ADU/e-"
        self.xlabel  = "row number"
        self.ylabel  = "$\\textit{Mean Charge} \; [e^{-}/pix]$"
        self._ylabel = "Mean Pixel Charge [e-/pix]"

        self.n_mad    = 5

        self.__units__.update({'image':1, 'n_mad':1})

    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(self.n_mad,ADC2e)

    def execute_process(self,rawdata,**kwargs):
        """Estimate the dark current per row as the median of those pixels above a given threshold,
        which is given by
            - q_median_row + N q_mad_row

        """
        # set output directory 
        self.output = rawdata.me_output

        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        image_amp = rawdata.get_image_by_attribute_name(self.image)

        # ONLY ACTIVE REGION, AND CHARGE IN UNITS OF ELECTRONS
        image_amp_e = np.ma.array(image_amp, mask=rawdata.mask_image_active_region) / u._gain[rawdata.execute_process_in_amp]
        
        # GET THE PIXEL TRHESHOLD CHARGE PER ROW
        q_median_row, q_mad_row = self.median_and_mad(image_amp_e, axis=1, keepdims=1)
        q_threshold_row = q_median_row + self.n_mad * q_mad_row
        
        # AVERAGE AND STANDARD DEVIATION OF THOSE PIXELS WITH CHARGE BELOW CERTAIN VALUE
        ma_zeros = np.ma.zeros( q_mad_row.shape, dtype=np.float64)
        ma_zeros.mask = q_mad_row.mask
        dc_mean_row = ma_zeros.copy()
        dc_std_row  = ma_zeros.copy()
        rows = np.arange(rawdata.rows[0],rawdata.rows[1])
        
        for r in rows:
            q_row = image_amp_e[r,:]
            dc_mean_row[r] = np.ma.average( q_row[q_row < q_threshold_row[r]] )
            dc_std_row[r]  = np.ma.std( q_row[q_row < q_threshold_row[r]] )

        ### ADD DATA: dc per row
        rows = rows + 1
        self.add_data( rows.tolist(), dc_mean_row.compressed().tolist(), rawdata.execute_process_in_amp, yerr=dc_std_row.compressed().tolist() )
        # XXX STD can be assumed as an error of the mean DC value, however we do not add this value,
        # as sometimes is too large and do not allow us to see the 'local' behaviour of the mean
        # values XXX by default the error is not plotted in the final plot
        
        ### SOME OF THE PARAMETERS COMPUTED IN THIS ME WILL BECAME ME TOO         ################################################################
        ##########################################################################################################################################
        setattr(rawdata,'me_{}_sensor_median_row'.format(rawdata.execute_process_in_amp), q_median_row.compressed().tolist())
        setattr(rawdata,'me_{}_sensor_mad_row'.format(rawdata.execute_process_in_amp), q_mad_row.compressed().tolist())
        setattr(rawdata,'me_{}_sensor_rows'.format(rawdata.execute_process_in_amp),rows.tolist())
        setattr(rawdata,'me_{}_dc_std_row'.format(rawdata.execute_process_in_amp), dc_std_row.compressed().tolist())

        ### FIT DARK CURRENT PER ROW TO A LINEAR MODEL, AND STORE THE BEST FIT VALUES, THIS WILL BE
        # THEN SET AS ME
        self.do_linear_fit(self.x[-1], self.y[-1], dc_std_row.compressed().tolist() )
        setattr(rawdata,'me_{}_dc_mean_row_slope'.format(rawdata.execute_process_in_amp),(self.fit_model.coef[0],self.fit_model_errors[0]))
        setattr(rawdata,'me_{}_dc_mean_row_interception'.format(rawdata.execute_process_in_amp),(self.fit_model.coef[1],self.fit_model_errors[1]))

        return

class MESTDPixelChargePerRow(MEAbstract):
    """Standard deviation of the dark current per row as the STD of the column pixel charge after removing pixels above certain value
    """
    __name__ = "MESTDPixelChargePerRow"
    __sequence_id__ = 5

    def __init__(self):
        super().__init__()

        self.title   = "Standard deviation of DC per row [see MEMeanPixelChargePerRow]"
        self.caption_fmt = "STD of the pixel charge values in a row, after excluding some pixels [see MEMeanPixelChargePerRow]. Calibration used: {} ADU/e-."
        self.xlabel  = "row number"
        self.ylabel  = "$\\textit{STD} \;[e^{-}/pix]$"
        self._ylabel = "STD [e-/pix]"
       
    
    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):
        
        # Check if the ME has been defined previously
        try:
            yvalues = getattr(rawdata,'me_{}_dc_std_row'.format(rawdata.execute_process_in_amp))
            xvalues = getattr(rawdata,'me_{}_sensor_rows'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            self.print_warning("MEMeanPixelChargePerRow is MANDATORY if MESTDPixelChargePerRow is booked.")
            return

        # set output directory
        self.output = rawdata.me_output

        self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp )
                #, rawdata.ccd,
                # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start)  )

        return


##########################################################################################################################################
################ ADD ALL ME THAT ARISE FROM THIS ONE
##########################################################################################################################################
class MESlopeFromMeanPCPerRow(MEAbstract):
    """Linear dependence of the median pixel charge per row (after removing pixels above certain
    threshold)
    """
    __name__ = "MESlopeFromMeanPCPerRow"
    __sequence_id__ = 6

    def __init__(self):
        super().__init__()
        
        # this ME will be used to find images out of the normality, i.e. is a RUN ME, but not an IMAGE ME
        self.is_MECCD = False

        self.title   = "Slope from the DC(row) fit [see MEMeanPixelChargePerRow]"
        self.caption_fmt = "Is the slope (charge variation with the rows) from the linear regresion of MEMeanPixelChargePerRow (mean charge versus row). Note that not all pixel values are included in MEMeanPixelChargePerRow (find information on his plot). Calibration used: {} ADU/e-."
        self.xlabel = "image file number"
        self.ylabel = "$\\textit{slope} \;[e^{-}/pix]$"
        self._ylabel = "Slope of DC(row) [e-/pix]"
        
        
    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):
        
        # Check if the ME has been defined previously
        try:
            yvalues = getattr(rawdata,'me_{}_dc_mean_row_slope'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            self.print_warning("MEMeanPixelChargePerRow is MANDATORY if MESlopeFromMeanPCPerRow is booked.")
            return

        # set output directory
        self.output = rawdata.me_output
        
        # ADD DATA
        self.add_data( int(rawdata.n_image), yvalues[0], rawdata.execute_process_in_amp, yerr=yvalues[1] )
                #, rawdata.ccd,
                # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start),
                
        return

class MEInterceptFromMeanPCPerRow(MEAbstract):
    """Dark current as the interception of the linear dependency of the averaged pixel charge versus
    rows
    """
    __name__ = "MEInterceptFromMeanPCPerRow"
    __sequence_id__ = 7

    def __init__(self):
        super().__init__()

        # this ME will be used to find images out of the normally, i.e. is a RUN ME, an IMAGE ME,
        # but not an CCD ME
        self.is_MECCD = False

        self.title   = "Interception point of the linear fit of MEMeanPixelChargePerRow"
        self.caption_fmt = "Is the interception point from the linear regresion of MEMeanPixelChargePerRow (mean charge versus row). Note that not all pixel values are included in MEMeanPixelChargePerRow (find information on his plot). Calibration used: {} ADU/e-."
        self.xlabel  = "image file number"
        self.ylabel  = "$\\textit{intercept} \; [e^{-}/pix]$"
        self._ylabel  = "Intercept of DC(row) [e-/pix]"

        # redefine the QTest for this ME
    
    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):

        # Check if the ME has been defined previously
        try:
            yvalues = getattr(rawdata,'me_{}_dc_mean_row_interception'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            self.print_warning("MEMeanPixelChargePerRow is MANDATORY if MEInterceptFromMeanPCPerRow is booked.")
            return

        # set output directory
        self.output = rawdata.me_output
        
        # ADD DATA
        self.add_data( int(rawdata.n_image), yvalues[0], rawdata.execute_process_in_amp, yerr=yvalues[1] )
                #, rawdata.ccd,
                #rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start),
                

        return

class MEMedianPixelChargePerRow(MEAbstract):
    """Mean pixel charge per row
    """

    __name__ = "MEMedianPixelChargePerRow"
    __sequence_id__ = 8

    def __init__(self):
        super().__init__()
        
        self.title   = "Median pixel charge (only active region)"
        self.caption_fmt = "Median value of the pixel charge belonging to the same row. Note that pixels above certain MAD have been excluded, mostly to ignore any possible cluster. Calibration used: {} ADU/e-."
        self.xlabel  = "row number"
        self.ylabel  = "$\\textit{Mean Charge} \; [e^{-}/pix]$"
        self._ylabel = "Mean Pixel Charge [e-/pix]"


    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):

        # Check if the ME has been defined previously
        try:
            yvalues = getattr(rawdata,'me_{}_sensor_median_row'.format(rawdata.execute_process_in_amp))
            xvalues = getattr(rawdata,'me_{}_sensor_rows'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            self.print_warning("MEMeanPixelChargePerRow is MANDATORY if MEMedianPixelChargePerRow is booked.")
            return
        
        # set output directory
        self.output = rawdata.me_output

        # ADD DATA
        self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp )
                #, rawdata.ccd,
                # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start))

        # fit data to linear fit 
        # self.do_linear_fit(self.x[-1], self.y[-1])

class MEMADPixelChargePerRow(MEAbstract):
    """Median Absolute Deviation of the column pixel charge as a funciton of the rows
    """

    __name__ = "MEMADPixelChargePerRow"
    __sequence_id__ = 9

    def __init__(self):
        super().__init__()
        
        self.title   = "Median pixel charge (only active region)"
        self.caption_fmt = "MAD of the pixel charge values in a row (similar to MESTDPixelChargePerRow, but considering all pixels in the sensitive region!). Calibration used: {} ADU/e-."
        self.xlabel  = "row number"
        self.ylabel  = "$\\textit{MAD(Pixel Charge)} \; [e^{-}/pix]$"
        self._ylabel = "MAD [e-/pix]"


    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):

        # Check if the ME has been defined previously
        try:
            yvalues = getattr(rawdata,'me_{}_sensor_mad_row'.format(rawdata.execute_process_in_amp))
            xvalues = getattr(rawdata,'me_{}_sensor_rows'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            self.print_warning("MEMeanPixelChargePerRow is MANDATORY if MEMADPixelChargePerRow is booked.")
            return
        
        # set output directory
        self.output = rawdata.me_output

        # ADD DATA
        self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp )
                # , rawdata.ccd,
                # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start))

        # fit data to linear fit 
        # self.do_linear_fit(self.x[-1], self.y[-1])


###############################################################################################################
#
#   PIXEL CHARGE AS A FUNCTION OF COLUMN
#
###############################################################################################################
class MEMedianPixelChargePerCol(MEAbstract):
    """Median Charge as the averaged rows pixel charge as a function of columns
    """
    
    __name__ = "MEMedianPixelChargePerCol"
    __sequence_id__ = 10

    def __init__(self):
        super().__init__()

        # this ME acts on the pedestal subtracted image
        self.image = 'mean_compressed_pedestal_subtracted'

        self.title  = "Median pixel charge (in the active region)"
        self.caption_fmt = "Median value of the pixel charge belonging to the same column. Note that pixels above certain MAD have been excluded, mostly to ignore any possible cluster. Calibration used: {} ADU/e-"
        self.xlabel = "column number"
        self.ylabel = "$\\textit{Median Charge} [e^{-}/pix]$"
        self._ylabel = "Median Pixel Charge [e-/pix]"

    
    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):
        
        # set output directory
        self.output = rawdata.me_output

        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        image_amp = rawdata.get_image_by_attribute_name(self.image)
        
        # ONLY ACTIVE REGION, AND CHARGE IN UNITS OF ELECTRONS
        image_amp_e = np.ma.array(image_amp, mask=rawdata.mask_image_active_region) / u._gain[rawdata.execute_process_in_amp]
        
        # GET THE PIXEL TRHESHOLD CHARGE PER COLUMN
        q_median_col, q_mad_col = self.median_and_mad(image_amp_e, axis=0, keepdims=1)
        
        # ADD DATA
        cols = np.arange(rawdata.cols[0],rawdata.cols[1]) + 1

        self.add_data( cols.tolist(), q_median_col.compressed().tolist(), rawdata.execute_process_in_amp ) #, yerr=q_mad_col.compressed().tolist() )
                # , rawdata.ccd,
                # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start),
                # 

        ### SOME OF THE PARAMETERS COMPUTED IN THIS ME WILL BECAME ME TOO         ################################################################
        ##########################################################################################################################################
        setattr(rawdata,'me_{}_sensor_mad_col'.format(rawdata.execute_process_in_amp), q_mad_col.compressed().tolist())
        setattr(rawdata,'me_{}_sensor_cols'.format(rawdata.execute_process_in_amp),cols.tolist())

        
class MEMADPixelChargePerCol(MEAbstract):
    """Median absolute deviation as the mad of the row pixel charge as a funciton of columns
    """

    __name__ = "MEMADPixelChargePerCol"
    __sequence_id__ = 11

    def __init__(self):
        super().__init__()

        self.title  = "Median Absolute Deviation (in the active region)"
        self.caption_fmt = "MAD of pixel charge values in a column (all pixels in the active region). Calibration used: {} ADU/e-."
        self.xlabel = "column number"
        self.ylabel = "$\\textit{MAD(Pixel Charge)} \; [e^{-}/pix]$"
        self._ylabel = "MAD [e-/pix]"

    
    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):
        
         # Check if the ME has been defined previously
        try:
            yvalues = getattr(rawdata,'me_{}_sensor_mad_col'.format(rawdata.execute_process_in_amp))
            xvalues = getattr(rawdata,'me_{}_sensor_cols'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            self.print_warning("MEMedianPixelChargePerCol is MANDATORY if MEMADPixelChargePerCol is booked.")
            return
        
        # set output directory
        self.output = rawdata.me_output

        # ADD DATA
        self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp )
                # rawdata.ccd,
                # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start) )


##### ME TO CHECK HOT COLUMNS
class MENpixQ1ePerCol(MEAbstract):
    """ME distribution of pixels with at least 1 electron (mask with a box any pixel with q>20e)
    """
    
    __name__ = "MENpixQ1ePerCol"
    __sequence_id__ = 12

    def __init__(self):
        super().__init__()
        self._per_amp = False

        # this ME acts on the pedestal subtracted image
        self.image = 'mean_compressed_pedestal_subtracted_e'
        self.ne = 1
        self.mask_clusters = False
        
        # in electrons
        self.Qhigh_pix = 20.5
        self.DCOLS = 5
        self.DROWS = 5

        self.title  = "Distribution of pixels with 1 electron (or more) after excluding any pixel in a cluster E>20e"
        self.caption= "Number of pixels with more than 1 electron per column. Note that pixels with loading values above 20 electrons have been masked with a vertical box: 1 col x 10 rows. The main goal of this plot is to highlight the hot columns, and clusters if CTI is present may contribute, reason for masking."
        self.xlabel = "column number"
        self.ylabel = "Npix(q>1e)"
        self._ylabel = "n. pix with q>0 e-"

        
        self.__units__.update({'image':1,'ne':1,'mask_clusters':1}) 

    def execute_process(self,rawdata,**kwargs):
        
        # set output directory
        self.output = rawdata.me_output

        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        image_amp = rawdata.get_image_by_attribute_name(self.image)
        
        # build the mask for the pseudo-clusters
        mask = np.zeros_like(image_amp.data)
        if self.mask_clusters:
            LY,LX = image_amp.data.shape
            rows,cols = np.where(image_amp.data>self.Qhigh_pix)
            for r,c in zip(rows,cols):
                mask[r:min(r+self.DROWS,LY),c:min(c+self.DCOLS,LX)] = True

        # masked data with charge > ne
        npix_q1e = (np.ma.round(np.ma.array(image_amp.data, mask=mask)) > self.ne).sum(axis=0)        
        # column density of events with q>1e
        x = np.arange(0,len(npix_q1e))
        
        self.add_data(x.tolist(), npix_q1e.tolist(), 'UL')
        
###############################################################################################################
#
#   MONITOR ELEMENTS TO CHECK THE OVERSCAN
#
###############################################################################################################
class MEOVSMedianPixelChargePerRow(MEAbstract):
    """Calculates the median of each row in the overscan as a function of the row.
    """

    __name__ = "MEOVSMedianPixelChargePerRow"
    __sequence_id__ = 20

    def __init__(self):
        super().__init__()

        # this ME acts on the pedestal subtracted image
        self.image = 'mean_compressed_pedestal_subtracted'

        self.title   = "Median Charge in the overscan."
        self.caption_fmt = "ONLY OVERSCAN REGION: Median pixel charge values in a row in units of e-. Calibration used: {} ADU/e-"
        self.xlabel  = "row number"
        self.ylabel  = "$\\textit{Mean Charge} \; [e^{-}/pix]$"
        self._ylabel = "overscan Mean Pixel Charge [e-/pix] "


        self.__units__.update({'image':1}) 
    
    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):

        # Do we have OVERSCAN REGION?
        if np.logical_not(rawdata.mask_image_overscan_cols).sum()==0:
            return

        # set output directory 
        self.output = rawdata.me_output

        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        image_amp = rawdata.get_image_by_attribute_name(self.image)

        # ONLY ACTIVE REGION, AND CHARGE IN UNITS OF ELECTRONS
        image_amp_e = np.ma.array(image_amp, mask=rawdata.mask_image_overscan_cols) / u._gain[rawdata.execute_process_in_amp]
        
        # median and mad pixel charge in the overscan region
        q_ovs_median_row, q_ovs_mad_row = self.median_and_mad(image_amp_e, axis=1, keepdims=1)
        
        ### ADD DATA: dc per row
        rows = np.arange(rawdata.rows[0],rawdata.rows[1]) + 1
        self.add_data( rows.tolist(), q_ovs_median_row.compressed().tolist(), rawdata.execute_process_in_amp, yerr=q_ovs_mad_row.compressed().tolist() )
                # , rawdata.ccd,
                # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start),
                

        
        ### SOME OF THE PARAMETERS COMPUTED IN THIS ME WILL BECAME ME TOO         ################################################################
        ##########################################################################################################################################
        setattr(rawdata,'me_{}_ovs_mad_row'.format(rawdata.execute_process_in_amp), q_ovs_mad_row.compressed().tolist())
        setattr(rawdata,'me_{}_ovs_rows'.format(rawdata.execute_process_in_amp),rows.tolist())

        ### FIT DARK CURRENT PER ROW TO A LINEAR MODEL, AND STORE THE BEST FIT VALUES, THIS WILL BE
        # THEN SET AS ME
        # self.do_linear_fit(self.x[-1], self.y[-1], q_ovs_mad_row.compressed().tolist())
        # now the fit is done automatically with method add_data ---- 
        setattr(rawdata,'me_{}_dc_ovs_row_slope'.format(rawdata.execute_process_in_amp),(self.fit_model.coef[0],self.fit_model_errors[0]))
        setattr(rawdata,'me_{}_dc_ovs_row_interception'.format(rawdata.execute_process_in_amp),(self.fit_model.coef[1],self.fit_model_errors[1]))
       
        return


class MEOVSMADPixelChargePerRow(MEAbstract): 
    """Median absolute deviation of the column pixel charge in the overscan as a function of the row
    """

    __name__ = "MEOVSMADPixelChargePerRow"
    __sequence_id__ = 21

    def __init__(self):
        super().__init__()

        self.title   = "Median Absolute Deviation of the pixel charge in the overscan."
        self.caption_fmt = "ONLY OVERSCAN REGION: MAD of the pixel charge values in a row in units of e-. Calibration used: {} ADU/e-"
        self.xlabel  = "row number"
        self.ylabel  = "$\\textit{MAD(Pixel Charge)} \; [e^{-}/pix]$"
        self._ylabel = "overscan MAD [e-/pix]"
        
    
    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):

         # Check if the ME has been defined previously
        try:
            yvalues = getattr(rawdata,'me_{}_ovs_mad_row'.format(rawdata.execute_process_in_amp))
            xvalues = getattr(rawdata,'me_{}_ovs_rows'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            self.print_warning("MEOVSMedianPixelChargePerRow is MANDATORY if MEOVSMADPixelChargePerRow is booked.")
            return
        
        # set output directory
        self.output = rawdata.me_output

        # ADD DATA
        self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp )
                # rawdata.ccd,
                # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start))

        return



###############################################################################################################
#
#   MONITOR ELEMENTS CREATED DURING PEDESTAL SUBTRACTION
#
###############################################################################################################
class MEOVSPedestalMuPerRow(MEAbstract):
    """Pedestal or baseline estimated by a Gaussian fit to the 0e- peak from the pixel charge distribution of the overscan region
    """

    __name__ = "MEOVSPedestalMuPerRow"
    __sequence_id__ = 30

    def __init__(self):
        super().__init__()

        self.title   = "Pedestal/baseline as the mean gaussian fit (Overscan PCD)"
        self.caption = "It is the mu parameter of the the Gaussian fit perfomed on the pixel charge distribution of a row. Possible horitzontal clusters are masked by applying a charge threshold. This pedestal is then subtracted row-by-row to the mean compressed image."
        self.xlabel = "row number"
        self.ylabel = "$\mu_{row}^{OVS} \;[ADU]$"
        self._ylabel = "Pedestal,mu [ADU]"

    
    def execute_process(self,rawdata,**kwargs):
        
        # check if the attribute has been calculated 
        try:
            yvalues = getattr(rawdata, 'pedestal_mu')
            yvalues = yvalues[rawdata.execute_process_in_amp]
        except AttributeError:
            self.print_warning("pedestal_mu is not found! Is PedestalSubtractedProcess done? ")
            return

        # pedestal can be subtracted in both axis (row:0, cols:1), 
        # this ME monitorize the pedestal estimated for each row
        try:
            yvalues = yvalues[0]
        except KeyError:
            self.print_warning("pedestal_mu is not found! Is PedestalSubtractedProcess done row by row? ")

        # Errors and chi2 information ony exists if pedestal subtraction was done with 'gaus_fit' method, if exists, include it
        try:
            yerr = np.concatenate(rawdata.pedestal_mu_error[rawdata.execute_process_in_amp][0]).ravel().tolist()
        except (AttributeError,KeyError):
            yerr = None

        xvalues = getattr(rawdata,'me_{}_ovs_rows'.format(rawdata.execute_process_in_amp))
        # mu is a list of numpys of 1 dimension, flattening the arrays
        yvalues = np.concatenate(yvalues).ravel().tolist() 
        apply_nan_mask = False
        if len(xvalues)!=len(yvalues):
            nan_mask = ~np.isnan(yvalues)
            apply_nan_mask = True
            yvalues = np.array(yvalues)[nan_mask]
            if len(xvalues)!=len(yvalues):
                raise RuntimeError(f"<self.__name__>: Length of y differs from x even ignoring nan parameters")
            yvalues = yvalues.tolist()
        
        # set output directory
        self.output = rawdata.me_output
         
        # ADD DATA
        if yerr is None:
            self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp ) 
            #, rawdata.ccd,
            # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start))
        else:
            if apply_nan_mask:
                self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp, yerr=np.array(yerr)[nan_mask].tolist() ) 
            else:
                self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp, yerr=yerr )

        # If the pedestal is given by the gaussian fit, associate to each single y-value the chi2 from the gaussian fit, 
        # i.e. overrride the chi2 given by the linear fit, done by the add_data method
        if yerr is not None:
            if apply_nan_mask:
                self.chi2 = rawdata.pedestal_chi2[rawdata.execute_process_in_amp][0][nan_mask].tolist()
            else:
                self.chi2 = rawdata.pedestal_chi2[rawdata.execute_process_in_amp][0].tolist()

        return


class MEOVSPedestalSigmaPerRow(MEAbstract):
    """Electron resolution as the sigma of the overscan pixel charge distribution estimated by a Gaussian fit to the 0e=- peak over the ovs PCD
    """

    __name__ = "MEOVSPedestalSigmaPerRow"
    __sequence_id__ = 31

    def __init__(self):
        super().__init__()

        self.title   = "Electronic noise characterization from the Gaussian fit over the overscan region"
        self.caption_fmt = "It is the sigma parameter of the the Gaussian fit perfomed on the pixel charge distribution of a row. Possible horitzontal clusters are masked by applying a charge threshold. Calibration used: {} ADU/e-."
        self.xlabel = "row number"
        self.ylabel = "$\sigma_{row}^{OVS} \;[ADU]$"
        self._ylabel = "Pedestal,STD [ADU]"
    

    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):
        
        # check if the attribute has been calculated 
        try:
            yvalues = getattr(rawdata, 'pedestal_sigma')
            yvalues = yvalues[rawdata.execute_process_in_amp]
        except AttributeError:
            self.print_warning("pedestal_sigma is not found! Is PedestalSubtractedProcess done? ")
            return

        # pedestal can be subtracted in both axis (row:0, cols:1), 
        # this ME monitorize the pedestal estimated for each row
        try:
            yvalues = yvalues[0]
        except KeyError:
            self.print_warning("pedestal_sigma is not found! Is PedestalSubtractedProcess done row by row? ")

        # Errors and chi2 information ony exists if pedestal subtraction was done with 'gaus_fit' method, if exists, include it
        try:
            yerr = np.concatenate(rawdata.pedestal_sigma_error[rawdata.execute_process_in_amp][0]).ravel().tolist()
        except (AttributeError,KeyError):
            yerr = None

        xvalues = getattr(rawdata,'me_{}_ovs_rows'.format(rawdata.execute_process_in_amp))
        # mu is a list of numpys of 1 dimension, flattening the arrays
        yvalues = np.concatenate(yvalues).ravel().tolist() 
        apply_nan_mask = False
        if len(xvalues)!=len(yvalues):
            nan_mask = ~np.isnan(yvalues)
            apply_nan_mask = True
            yvalues = np.array(yvalues)[nan_mask]
            if len(xvalues)!=len(yvalues):
                raise RuntimeError(f"<self.__name__>: Length of y differs from x even ignoring nan parameters")
            yvalues = yvalues.tolist()
       
        # set output directory
        self.output = rawdata.me_output
         
        # ADD DATA
        if yerr is None:
            self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp ) 
        else:
            if apply_nan_mask:
                self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp, yerr=np.array(yerr)[nan_mask].tolist() ) 
            else:
                self.add_data( xvalues, yvalues, rawdata.execute_process_in_amp, yerr=yerr )

        # If the pedestal is given by the gaussian fit, associate to each single y-value the chi2 from the gaussian fit, 
        # i.e. overrride the chi2 given by the linear fit, done by the add_data method
        if yerr is not None:
            if apply_nan_mask:
                self.chi2 = rawdata.pedestal_chi2[rawdata.execute_process_in_amp][0][nan_mask].tolist()
            else:
                self.chi2 = rawdata.pedestal_chi2[rawdata.execute_process_in_amp][0].tolist()
        
        return

class MEOVSPedestalMu(MEAbstract):
    """Pedestal or baseline estimated by a Gaussian fit to the 0e=- peak from the pixel charge distribution of the overscan region
    """

    __name__ = "MEOVSPedestalMu"
    __sequence_id__ = 32

    def __init__(self):
        super().__init__()

        self.title   = "Pedestal/baseline as the mean gaussian fit (Overscan PCD)"
        self.caption = "The pedestal has been estimated row by row by fitting a gaussian to the PCD of each row, being the best fitted mu. This ME is the average of all row-based pedestal values (per CCD, per amp) in units of ADU."
        self.xlabel = "time [days]"
        self.ylabel = "$\\textit{median}(\mu_{row}^{OVS}) \;[ADU]$"
        self._ylabel = "Mean(Pedestal,mu per row) [ADU]"


    def execute_process(self,rawdata,**kwargs):
        
        # check if the attribute has been calculated 
        try:
            yvalues = getattr(rawdata, 'pedestal_mu')
            yvalues = yvalues[rawdata.execute_process_in_amp]
        except AttributeError:
            self.print_warning("pedestal_mu is not found! Is PedestalSubtractedProcess done? ")
            return

        # pedestal can be subtracted in both axis (row:0, cols:1), 
        # this ME monitorize the pedestal estimated for each row
        try:
            yvalues = yvalues[0]
        except KeyError:
            self.print_warning("pedestal_mu is not found! Is PedestalSubtractedProcess done row by row? ")
        
        # mu is a list of numpys of 1 dimension, flattening the arrays
        yvalues = np.concatenate(yvalues).ravel().tolist() 
        
        # set output directory
        self.output = rawdata.me_output

        # ADD DATA
        self.add_data( rawdata.start, np.nanmean(yvalues), rawdata.execute_process_in_amp ) 

        setattr(rawdata,'QT_{}_{}'.format(rawdata.ccd, rawdata.execute_process_in_amp), self.execute_qtest(yvalues, rawdata.ccd, rawdata.execute_process_in_amp, 'QT:range') )

        return

class MEOVSPedestalSigma(MEAbstract):
    """Electron resolution as the sigma of the overscan pixel charge distribution estimated by a Gaussian fit to the 0e=- peak over the ovs PCD
    """

    __name__ = "MEOVSPedestalSigma"
    __sequence_id__ = 33

    def __init__(self):
        super().__init__()

        self.title   = "Resolution as the sigma gaussian fit (Overscan PCD)"
        self.caption = "The pedestal has been estimated row by row by fitting a gaussian to the PCD of each row. This ME refers to the best fitted sigma values, being the average of the row-based sigma values (per CCD, per amp) in units of ADU."
        self.xlabel  = "time [days]"
        self.ylabel  = "$\\textit{median}(\sigma_{row}^{OVS}) \;[ADU]$"
        self._ylabel = "Mean(Pedestal,sigma per row) [ADU]"


    def execute_process(self,rawdata,**kwargs):
        
        # check if the attribute has been calculated 
        try:
            yvalues = getattr(rawdata, 'pedestal_sigma')
            yvalues = yvalues[rawdata.execute_process_in_amp]
        except AttributeError:
            self.print_warning("pedestal_sigma is not found! Is PedestalSubtractedProcess done? ")
            return

        # pedestal can be subtracted in both axis (row:0, cols:1), 
        # this ME monitorize the pedestal estimated for each row
        try:
            yvalues = yvalues[0]
        except KeyError:
            self.print_warning("pedestal_sigma is not found! Is PedestalSubtractedProcess done row by row? ")

        # mu is a list of numpys of 1 dimension, flattening the arrays
        yvalues = np.concatenate(yvalues).ravel().tolist() 
        
        # set output directory
        self.output = rawdata.me_output
         
        # ADD DATA
        self.add_data( rawdata.start, np.nanmean(yvalues), rawdata.execute_process_in_amp) 

        setattr(rawdata,'QT_{}_{}'.format(rawdata.ccd, rawdata.execute_process_in_amp), self.execute_qtest(yvalues, rawdata.ccd, rawdata.execute_process_in_amp, 'QT:range') )
       
        return

class MEOVSPCD(MEAbstract):
    """
    """
    __name__ = "MEOVSPCD"
    __sequence_id__ = 40

    def __init__(self):
        super().__init__()

        self.image = "mean_compressed"
        
        self.do_histogram = False
        self.do_stacked_hist = True

        self.title   = "Pixel Charge distribution (overscan)"
        self.caption_fmt = "Distribution of the pixel charge in the overscan region. Calibration used: {} ADU/e-."
        self.xlabel  = "Pixel Charge [e-]"
        self.ylabel  = "counts"
        self._ylabel = self.ylabel


        self.__units__.update({'image':1,'do_histogram':1, 'do_stacked_hist':1})

    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):

        # Do we have OVERSCAN REGION?
        if np.logical_not(rawdata.mask_image_overscan_cols).sum()==0:
            return

        # set output directory 
        self.output = rawdata.me_output

        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        image_amp = rawdata.get_image_by_attribute_name(self.image)

        # ONLY ACTIVE REGION, AND CHARGE IN UNITS OF ELECTRONS
        image_ovs_amp = np.ma.array(image_amp, mask=rawdata.mask_image_overscan_cols) / u._gain[rawdata.execute_process_in_amp]

        # ADD DATA
        # THE ME CAN CREATE THE HISTOGRAM AND APPEND THE FREQUENCIES TO BE PLOTTED, OR ACCUMULATE
        # PIXEL CHARGE TO FINALLY DO AN HISTOGRAM DISTRIBUTION
        if self.do_histogram:
            # histogram
            yvalues, xvalues = np.histogram( image_ovs_amp.compressed(), int(np.sqrt(image_ovs_amp.size)) )
            xcenters = xvalues[1:] - np.diff(xvalues)[0]
            
            self.add_data( xcenters.tolist(), yvalues.tolist(), rawdata.execute_process_in_amp ) # , rawdata.ccd,
                    # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image) )
        else:
            xvalues = image_ovs_amp.compressed().tolist()
            self.add_data( xvalues, [None]*len(xvalues), rawdata.execute_process_in_amp ) #, rawdata.ccd,
                    #rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start))

        return

    def execute_plot(self,**kwargs):
        if self.do_histogram:
            super().execute_plot(ptype='bar')    
        else:
            xmax,xmin = np.array(self.x,dtype=object).ravel().max(),np.array(self.x,dtype=object).ravel().min()
            super().execute_plot( ptype='hist',
                    **{'nbins':int((np.max(xmax)-np.min(xmin))*50),
                        'fill':False, 'histtype':"step", 'stacked':self.do_stacked_hist })



class MEPCD(MEAbstract):
    """
    """
    __name__ = "MEPCD"
    __sequence_id__ = 41

    def __init__(self):
        super().__init__()
        
        self.image = "mean_compressed_pedestal_subtracted"
        
        self.threshold = 500*u.e
        # type of plot
        self.do_histogram    = False
        self.do_stacked_hist = True

        self.title   = "Pixel Charge distribution (active region)"
        self.caption_fmt = "Distribution of the pixel charge in the active region. Calibration used: {} ADU/e-."
        self.xlabel = "Pixel Charge [e-]"
        self.ylabel  = "counts"
        self._ylabel = self.ylabel


        self.__units__.update({'image':1,'threshold':u.e,'do_histogram':1, 'do_stacked_hist':1})

    @property
    def caption(self):
        ADC2e = ','.join('='.join((str(key),str(val))) for (key,val) in u._gain.items())
        return self.caption_fmt.format(ADC2e)

    def execute_process(self,rawdata,**kwargs):
        """Append pixel charge to the x-values. The histogram distribution will be done, during the
        execute_plot
        """

        # Do we have OVERSCAN REGION?
        if np.logical_not(rawdata.mask_image_overscan_cols).sum()==0:
            return

        # set output directory 
        self.output = rawdata.me_output

        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        image_amp = rawdata.get_image_by_attribute_name(self.image)

        # ONLY ACTIVE REGION, AND CHARGE IN UNITS OF ELECTRONS
        image_amp = np.ma.array(image_amp, mask=rawdata.mask_image_active_region) / u._gain[rawdata.execute_process_in_amp]
        
        # ADD DATA
        qmax = np.ma.min(image_amp) + self.threshold
        qvalues = image_amp[image_amp<qmax].compressed().tolist()
        # THE ME CAN CREATE THE HISTOGRAM AND APPEND THE FREQUENCIES TO BE PLOTTED, OR ACCUMULATE
        # PIXEL CHARGE TO FINALLY DO AN HISTOGRAM DISTRIBUTION
        if self.do_histogram:
            yvalues, xvalues = np.histogram(qvalues, int(np.sqrt(image_amp.size)) )
            xcenters = xvalues[1:] - np.diff(xvalues)[0]
            self.add_data( xcenters.tolist(), yvalues.tolist(), rawdata.execute_process_in_amp ) # , rawdata.ccd,
                    # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start))
        else:
            self.add_data( qvalues, [None]*len(qvalues), rawdata.execute_process_in_amp ) #, rawdata.ccd,
                    # rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start))

        return

    def execute_plot(self,**kwargs):
        if self.do_histogram:
            super().execute_plot(ptype='bar')
        else:
            xmax,xmin = np.array(self.x,dtype=object).ravel().max(),np.array(self.x,dtype=object).ravel().min()
            super().execute_plot( ptype='hist',
                    **{'nbins':int((np.max(xmax)-np.min(xmin))*50),
                        'fill':False,'stacked':self.do_stacked_hist})


#####################################################################################################
#
#   FIT DARK CURRENT ME
#
#####################################################################################################
class MEFitDCPerRow(MEAbstract):
    """
    """

    __name__ = "MEFitDCPerRow"
    __sequence_id__ = 50

    def __init__(self):
        super().__init__() 
        
        self.is_MERun = False

        self.title  = "Dark Current Fit row by row"
        self.caption= "Dark current from the DCFit done row by row (in units of electrons/pix/img)"
        self.xlabel = "rows"
        self.ylabel = "$\lambda \; [e^{-}/pix/img]$"
        self._ylabel = "Poisson DC [e-/pix/img]"

        self.image = "mean_compressed_pedestal_subtracted"
        
        # Attributes for the fitting algorithm
        self.row_start = 2
        self.row_end   = -1
        self.row_step  = 10
        
        # parameters for the DC fit
        self.cols_to_mask = []
        self.rows_to_mask = []      
        self.n_peaks                    = 3
        self.n_sigma_fit                = 10
        self.calibration                = 10.0
        self.binning_size               = 0.5
        self.x_min                      = -15
        self.x_max                      = 60
        self.mu_min,self.mu_max         = -1,1
        self.sigma_min,self.sigma_max   = 0.01,3.0
        self.dc_min,self.dc_max         = 0.0001,50
        self.gain_min,self.gain_max     = 8.0,15.0


        self.__units__.update({'image':1,
            'rows_to_mask':1, 'cols_to_mask':1,'row_start':1, 'row_end':1, 'row_step':1,
            'calibration':1,'n_peaks':1,'n_sigma_fit':1,'max_img_pcd':1,
            'binning_size':1,
            'sigma_min':1,'sigma_max':1,
            'mu_min':1,'mu_max':1,
            'dc_min':1,'dc_max':1,
            'gain_min':1,'gain_max':1,
            'x_min':1,'x_max':1})

    
    def execute_process(self,rawdata,**kwargs):
        """
        """
        # set output directory 
        self.output = rawdata.me_dcfit

        dcfit_row = FitDarkCurrentPerRow()
        dcfit_row.__silent__    = True
        for par in ['image', 'n_peaks', 'calibration', 'x_min', 'x_max', 'mu_min', 'mu_max', 'sigma_min','sigma_max', 'dc_min',
                'dc_max', 'gain_min', 'gain_max','rows_to_mask', 'cols_to_mask','row_start', 'row_end', 'row_step']:
            setattr(dcfit_row, par, getattr(self,par))

        results = dcfit_row.execute_process(rawdata,**{'results':True})
        
        x = results['xrows']
        self.add_data( results['xrows'][0], results['dc'][0],
                rawdata.execute_process_in_amp, yerr=results['dc'][1],
                xerr=results['xrows'][1]) 
 
        setattr(rawdata,'me_DCfit_row_{}'.format(rawdata.execute_process_in_amp),results)


class MEFitSigmaPerRow(MEAbstract):
    """
    """

    __name__ = "MEFitSigmaPerRow"
    __sequence_id__ = 51

    def __init__(self):
        super().__init__() 
        
        self.is_MERun = False

        self.title  = "Resolution fitted row by row"
        self.caption= "Resolution as the sigma from the gaussian fit to the zero-electron peak in the row by row dark current fit."
        self.xlabel = "rows"
        self.ylabel = "$\\textit{resolution}, \; \sigma_0 [e^{-}]$"
        self._ylabel = "zero-e Gaussian Sigma [e-]"


        self.__units__.update({})

    
    def execute_process(self,rawdata,**kwargs):
        """
        """
        # set output directory 
        self.output = rawdata.me_dcfit

        try:
            results = getattr(rawdata,'me_DCfit_row_{}'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            return

        self.add_data( results['xrows'][0], results['sigma'][0],
                rawdata.execute_process_in_amp, yerr=results['sigma'][1],
                xerr=results['xrows'][1]) 
 

class MEFitCalibrationPerRow(MEAbstract):
    """
    """

    __name__ = "MEFitCalibrationPerRow"
    __sequence_id__ = 52

    def __init__(self):
        super().__init__() 
        
        self.is_MERun = False

        self.title  = "Calibration constant fitted row by row"
        self.caption= "Calibration as the distance between peaks in the row by row dark current fit."
        self.xlabel = "rows"
        self.ylabel = "$\\textit{calibration} \; [ADU/e^{-}]$"
        self._ylabel = "Calibration [ADU/e-]"


        self.__units__.update({})

    
    def execute_process(self,rawdata,**kwargs):
        """
        """
        # set output directory 
        self.output = rawdata.me_dcfit

        try:
            results = getattr(rawdata,'me_DCfit_row_{}'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            return

        self.add_data( results['xrows'][0], results['calibration'][0],
                rawdata.execute_process_in_amp, yerr=results['calibration'][1],
                xerr=results['xrows'][1]) 
 

class MEFitMu0PerRow(MEAbstract):
    """
    """

    __name__ = "MEFitMu0PerRow"
    __sequence_id__ = 53

    def __init__(self):
        super().__init__() 
        
        self.is_MERun = False

        self.title  = "Center of the zero-electron peak fitted row by row"
        self.caption= "Center of the zero-electron peak as the mu from the gaussian on the row by row dark current fit."
        self.xlabel = "rows"
        self.ylabel = "$\mu_0 \; [e^{-}]$"
        self._ylabel = "zero-e Peak Center [e-]"


        self.__units__.update({})

    
    def execute_process(self,rawdata,**kwargs):
        """
        """
        # set output directory 
        self.output = rawdata.me_dcfit

        try:
            results = getattr(rawdata,'me_DCfit_row_{}'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            return

        self.add_data( results['xrows'][0], results['mu'][0],
                rawdata.execute_process_in_amp, yerr=results['mu'][1],
                xerr=results['xrows'][1]) 
 


class MEFitDCPerCol(MEAbstract):
    """
    """

    __name__ = "MEFitDCPerCol"
    __sequence_id__ = 54

    def __init__(self):
        super().__init__() 
        
        self.is_MERun = False

        self.title  = "Dark Current Fit column by column"
        self.caption= "Dark current from the DCFit done column by column (in units of electrons/pix/img)"
        self.xlabel = "columns"
        self.ylabel = "$\lambda \; [e^{-}/pix/img]$"
        self._ylabel = "Poisson DC [e-/pix/img]"

        self.col_start = 2
        self.col_end   = -1
        self.col_step  = 5

        # parameters for the DC fit
        self.image = "mean_compressed_pedestal_subtracted"
        self.cols_to_mask = []
        self.rows_to_mask = []      
        self.n_peaks                    = 3
        self.n_sigma_fit                = 10
        self.calibration                = 10.0
        self.binning_size               = 0.5
        self.x_min                      = -15
        self.x_max                      = 60
        self.mu_min,self.mu_max         = -1,1
        self.sigma_min,self.sigma_max   = 0.01,3.0
        self.dc_min,self.dc_max         = 0.0001,50
        self.gain_min,self.gain_max     = 8.0,15.0


        self.__units__.update({'image':1,
            'rows_to_mask':1, 'cols_to_mask':1,'col_start':1, 'col_end':1, 'col_step':1,
            'calibration':1,'n_peaks':1,'n_sigma_fit':1,'max_img_pcd':1,
            'binning_size':1,
            'sigma_min':1,'sigma_max':1,
            'mu_min':1,'mu_max':1,
            'dc_min':1,'dc_max':1,
            'gain_min':1,'gain_max':1,
            'x_min':1,'x_max':1})

    
    def execute_process(self,rawdata,**kwargs):
        """
        """
        # set output directory 
        self.output = rawdata.me_dcfit

        dcfit_col = FitDarkCurrentPerCol()
        dcfit_col.__silent__    = True
        for par in ['image', 'n_peaks', 'calibration', 'x_min', 'x_max', 'mu_min', 'mu_max', 'sigma_min','sigma_max', 'dc_min',
                'dc_max', 'gain_min', 'gain_max','rows_to_mask', 'cols_to_mask','col_start',
                'col_end', 'col_step']:
            setattr(dcfit_col, par, getattr(self,par))


        results = dcfit_col.execute_process(rawdata,**{'results':True})
        
        x = results['xcols']
        self.add_data( results['xcols'][0], results['dc'][0],
                rawdata.execute_process_in_amp, yerr=results['dc'][1],
                xerr=results['xcols'][1]) 
 
        setattr(rawdata,'me_DCfit_col_{}'.format(rawdata.execute_process_in_amp),results)


class MEFitSigmaPerCol(MEAbstract):
    """
    """

    __name__ = "MEFitSigmaPerCol"
    __sequence_id__ = 55

    def __init__(self):
        super().__init__() 
        
        self.is_MERun = False

        self.title  = "Resolution fitted column by column"
        self.caption= "Resolution as the sigma from the gaussian fit to the zero-electron peak in the column by column dark current fit."
        self.xlabel = "columns"
        self.ylabel = "$\\textit{resolution},\; \sigma_0 [e^{-}]$"
        self._ylabel = "zero-e Gaussian Sigma [e-]"


        self.__units__.update({})

    
    def execute_process(self,rawdata,**kwargs):
        """
        """
        # set output directory 
        self.output = rawdata.me_dcfit

        try:
            results = getattr(rawdata,'me_DCfit_col_{}'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            return

        self.add_data( results['xcols'][0], results['sigma'][0],
                rawdata.execute_process_in_amp, yerr=results['sigma'][1],
                xerr=results['xcols'][1]) 
 

class MEFitCalibrationPerCol(MEAbstract):
    """
    """

    __name__ = "MEFitCalibrationPerCol"
    __sequence_id__ = 56

    def __init__(self):
        super().__init__() 
        
        self.is_MERun = False

        self.title  = "Calibration constant fitted column by column"
        self.caption= "Calibration as the distance between peaks in the column by column dark current fit."
        self.xlabel = "columns"
        self.ylabel = "$\\textit{calibration} \; [ADU/e^{-}]$"
        self._ylabel = "Calibration [ADU/e-]"


        self.__units__.update({})

    
    def execute_process(self,rawdata,**kwargs):
        """
        """
        # set output directory 
        self.output = rawdata.me_dcfit

        try:
            results = getattr(rawdata,'me_DCfit_col_{}'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            return

        self.add_data( results['xcols'][0], results['calibration'][0],
                rawdata.execute_process_in_amp, yerr=results['calibration'][1],
                xerr=results['xcols'][1]) 
 

class MEFitMu0PerCol(MEAbstract):
    """
    """

    __name__ = "MEFitMu0PerCol"
    __sequence_id__ = 57

    def __init__(self):
        super().__init__() 
        
        self.is_MERun = False

        self.title  = "Center of the zero-electron peak fitted column by column"
        self.caption= "Center of the zero-electron peak as the mu from the gaussian on the column by column dark current fit."
        self.xlabel = "columns"
        self.ylabel = "$\mu_0 \; [e^{-}]$"
        self._ylabel = "zero-e- Peak Center [e-]"

        self.__units__.update({})

    
    def execute_process(self,rawdata,**kwargs):
        """
        """
        # set output directory 
        self.output = rawdata.me_dcfit

        try:
            results = getattr(rawdata,'me_DCfit_col_{}'.format(rawdata.execute_process_in_amp))
        except AttributeError:
            return

        self.add_data( results['xcols'][0], results['mu'][0],
                rawdata.execute_process_in_amp, yerr=results['mu'][1],
                xerr=results['xcols'][1]) 


#############################################################################################################
#
#
#############################################################################################################
class MEFitDC(MEAbstract):
    """
    """

    __name__ = "MEFitDC"
    __sequence_id__ = 60

    def __init__(self):
        super().__init__() 
        
        self.is_MERun = False

        self.title  = "Dark Current and calibration Fit"
        self.caption= "The distribution of pixel values (PCD) in the active region is fitted by the convolution of a Poisson (DC) and a Gaussian (single e- resolution, zero-e peak position). The fit is done in units of ADU being the calibration (distance between peaks) a free parameter during the fit, and assumes linearity with the number of electrons."
        self.xlabel = "Pixel Charge [ADU]"
        self.ylabel = "counts"
        self._ylabel = self.ylabel

        self.image = "mean_compressed_pedestal_subtracted"

        # Attributes for the fitting algorithm
        self.n_peaks = 3
        self.n_sigma_fit = 10
        self.calibration = 10.0
        self.binning_size = 0.5
        self.x_min = -20
        self.x_max = 60
        self.mu_min,self.mu_max         = -1,1
        self.sigma_min,self.sigma_max   = 0.01,3.0
        self.dc_min,self.dc_max         = 0.001,50
        self.gain_min,self.gain_max     = 8.0,15.0

        #  mask rows and columns
        self.rows_to_mask = []
        self.cols_to_mask = []

        self.mask_clusters = False
        
        # add parameter to update calibrated images after this fitting
        self.recalibrate = False

        self.__units__.update({'image':1,
            'calibration':1,'n_peaks':1,'n_sigma_fit':1,'max_img_pcd':1,
            'binning_size':1,
            'sigma_min':1,'sigma_max':1,
            'mu_min':1,'mu_max':1,
            'dc_min':1,'dc_max':1,
            'gain_min':1,'gain_max':1,
            'x_min':1,'x_max':1,
            'rows_to_mask':1, 'cols_to_mask':1,
            'mask_clusters':1,
            'recalibrate':1
            })

    
    def execute_process(self,rawdata,**kwargs):
        """
        """
        # set output directory 
        self.output = rawdata.me_dcfit

        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        image_amp = rawdata.get_image_by_attribute_name(self.image)

        # get DCFit object and Fit
        dcfit = self.get_dcfit_process()
        dcfit.execute_process(rawdata)

        # ADD ME THAT CAME OUT FROM THIS PROCESS (all in electron units!)
        results = dcfit.get_results_in_electrons(rawdata)
        setattr(rawdata,'me_dcfit_results_{}'.format(rawdata.execute_process_in_amp),results )
        # from results get chi2 for a binned L fit: Baker-Cousins Chi2
        self.chi2 = results['chi2']

        # ADD DATA
        x,y,y_fit = dcfit.get_dcfit_data()
        self.add_data( x, y, rawdata.execute_process_in_amp, y_fit=y_fit )

        # execute plot needs the dcfit 
        self.dcfit = dcfit.pcd_hist
        self.output += rawdata._file_name.replace(".{}".format(rawdata.__type__),"")
        # text to inclue in the plot
        nskips = int(rawdata.image_header['NSKIPS']) if 'NSKIPS' in rawdata.image_header else rawdata.nskips
        if nskips==1 and rawdata._NDCM>1:
            nskips = rawdata._NDCM
        # text_res=[
        #         "  #sigma="+str(round(results['sigma'][0],6))+"e^{-}/pix" ,
        #         "  #lambda="+str(round(results['dc'][0],6))+" e^{-}/pix/img",
        #         "  #chi^{2}_{red}[LL]="+str(round(self.chi2,3)),
        #         "bin: {}x{}".format(rawdata.npbin,rawdata.nsbin)+", "+"t_{exp}="+str(round(rawdata.exposure_time,1))+"s , t_{read}="+str(round(rawdata.read_time,1))+"s"]
        text_res = [
                f"image ID: {rawdata.n_image_from_file}",
                f"CCD-{rawdata.execute_process_in_amp}",
                f"frame: {rawdata.nrows}x{rawdata.ncols}x{nskips}",
                f"bin: {rawdata.npbin}x{rawdata.nsbin}x{nskips}",
                "t_{exp}="+str(round(rawdata.exposure_time,1))+"s",
                "t_{read}="+str(round(rawdata.read_time,1))+"s"
                ]
        wr_darkcurrent = self.execute_fit_plot(text=text_res,amp=rawdata.execute_process_in_amp)
        rawdata.wr_darkcurrent.append( wr_darkcurrent.strip() )

        # set attributes to write down 
        self.attributes_to_update_fits_header(dcfit,rawdata.execute_process_in_amp)
        
        # recalibrate the images
        cal = CalibrationProcess()
        cal.__silent__ = True
        cal.image = self.image.replace("image_","").replace("_e","")
        print("DQM INFO: Image {} has been recalibrated using the best fitted gain {} [amplifier {}]".format(cal.image,dcfit.dc_gain,rawdata.execute_process_in_amp))
        #if type(cal.gain) is not dict:
        #    # create a dictionary in case is not 
        #    cal.gain = {}
        cal.gain = {}
        cal.gain[rawdata.execute_process_in_amp] = dcfit.dc_gain
        cal.execute_process(rawdata)


    def get_dcfit_process(self):
        # function to call the FitDarkCurrentProcess
        dcfit = FitDarkCurrentProcess()
        # reduce displays
        dcfit.__silent__    = True
        dcfit.__display__   = False
        dcfit.__DEBUG__     = False
        dcfit.__verbose__   = False

        dcfit.save_plots    = False
        dcfit.do_calibration= True
        dcfit.image         = self.image
        dcfit.rows_to_mask  = self.rows_to_mask
        dcfit.cols_to_mask  = self.cols_to_mask
        
        # Q:quite, L:binned data with low-stats mostly for second peak of the poisson, S: to use the chi-2 for the goodness of the fit
        # dcfit.fit_options = "QLES"

        for par in ['n_peaks','calibration','binning_size', 'x_min', 'x_max', 'mu_min', 'mu_max', 'sigma_min','sigma_max', 'dc_min',
                'dc_max', 'gain_min', 'gain_max','rows_to_mask', 'cols_to_mask','mask_clusters']:
            setattr(dcfit,par,getattr(self,par))


        return dcfit

    def attributes_to_update_fits_header(self,dcfit,amp):
        """
        """
        u._tohdr.update({
            'MEGAIN{}'.format(amp)  :(dcfit.dc_gain,'DQM: fitted gain ADU/e-'),
            'MEDC{}'.format(amp)    :(dcfit.dc_lambda,'DQM: fitted DC ADU/bin/img'),
            'MESIG{}'.format(amp)   :(dcfit.dc_sigma,'DQM: fitted sigma ADU'),
            'MEMU0{}'.format(amp)   :(dcfit.dc_mu0,'DQM: fitted zero-e peak center in ADU'),
            'EMEGAIN{}'.format(amp) :(dcfit.dc_gain_err,'DQM: error of MEGAIN'),
            'EMEDC{}'.format(amp)   :(dcfit.dc_lambda_err,'DQM: error of MELAMBDA'),
            'EMESIG{}'.format(amp)  :(dcfit.dc_sigma_err,'DQM: error of MESIGMA'),
            'EMEMU0{}'.format(amp)  :(dcfit.dc_mu0_err,'DQM: error of MEMU0')
            })

    def execute_fit_plot(self,text=None,amp=None,**kwargs):
        """
        """
        ROOT.gROOT.SetBatch(1)

        c = ROOT.TCanvas("DCFit_FinalFit")
        self.dcfit.GetXaxis().SetTitle(" pixel charge [ADC]")
        self.dcfit.GetYaxis().SetTitle(" counts ")
        self.dcfit.Draw()
        c.Update()

        ### re-place the TPaveStats
        stat_box = self.dcfit.FindObject("stats")
        stat_box.SetX1NDC(0.62)
        stat_box.SetY1NDC(0.63)
        stat_box.SetX2NDC(0.96)
        stat_box.SetY2NDC(0.91)
        
        if text is not None:
            for i,t in enumerate(text):
                ttext = ROOT.TLatex()
                ttext.SetTextSize(0.03)
                ttext.DrawLatexNDC(0.55,0.33+0.06*i, "#color[810]{"+t+"}")

        c.SetLogy()
        c.Update()
        c.Draw()

        outfile = "{}_amp_{}_{}.png".format(self.output,amp,self.__name__)
        c.SaveAs(outfile)
        self.fig_png_file.append(outfile)

        return outfile.split('/')[-1]
    
    def execute_plot(self):
        """
        """
        # nothing to do, all plots have been done previously, one by one with the single images

        return


class MEFitDCMu0(MEAbstract): 
    """Median absolute deviation of the column pixel charge in the overscan as a function of the row
    """

    __name__ = "MEFitDCMu0"
    __sequence_id__ = 61

    def __init__(self):
        super().__init__()
        
        self.is_MECCD = False

        self.title   = "Position of the 0-electron peak (from DC fit)"
        self.caption = "This is the center of the first gaussian (i.e. the zero-electron peak) of the pixel charge distribution (PCD) that is fitted to a set of gaussians convolved by a poisson. A systematic offset of the position of this peak may indicate a problem on the pedestal estimation. "
        self.xlabel  = "image number"
        self.ylabel  = "$\\textit{zero-electron peak}, \; \mu_{0} \; [e-]$"
        self._ylabel = "zero-e- Peak Center [e-]"
        
        # This does not depend on the CCD, or science run, so rnages are defined here 
        self.plot_with_yerrors = [-0.01,0.01]


    def execute_process(self,rawdata,**kwargs):

         # Check if the ME has been defined previously
        try:
            results = getattr(rawdata,'me_dcfit_results_{}'.format(rawdata.execute_process_in_amp))
            yvalue  = results['mu0']
        except AttributeError:
            self.print_warning("MEFitDC is MANDATORY if MEFitDCMu0 is booked.")
            return
        
        # set output directory
        self.output = rawdata.me_output
        self.chi2 = results['chi2']
        # ADD DATA
        self.add_data( int(rawdata.n_image), yvalue[0], rawdata.execute_process_in_amp, yerr=yvalue[1] ) #, rawdata.ccd,
                #rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image),
                #datetime_to_days(rawdata.start), 

        setattr(rawdata,'QT_{}_{}'.format(rawdata.ccd, rawdata.execute_process_in_amp), self.execute_qtest(yvalue[0], rawdata.ccd, rawdata.execute_process_in_amp,'QT:chi2,range',**{'chi2':self.chi2}) )

        return


class MEFitDCSigma(MEAbstract): 
    """Median absolute deviation of the column pixel charge in the overscan as a function of the row
    """

    __name__ = "MEFitDCSigma"
    __sequence_id__ = 62

    def __init__(self):
        super().__init__()
        
        self.is_MECCD = False

        self.title   = "Single electron resolution (from DC fit)"
        self.caption = "This is the readout noise obtained from the fit of the pixel charge distribution (PCD) to a set of gaussians convolved by a poisson. The readout noise is the sigma from the Gaussian distribution, and is assumed to be independent of the number of electrons, so the zero-electron peak is the one that contributes the most, just for statistics. It is given in units of e-."
        self.xlabel  = "image number"
        self.ylabel  = "$\\textit{single-e resolution}, \; \sigma \; [e^{-}]$"
        self._ylabel = "zero-e Gaussian Sigma [e-]"
        

    def execute_process(self,rawdata,**kwargs):

         # Check if the ME has been defined previously
        try:
            results = getattr(rawdata,'me_dcfit_results_{}'.format(rawdata.execute_process_in_amp))
            yvalue  = results['sigma']
        except AttributeError:
            self.print_warning("MEFitDC is MANDATORY if MEFitDCMu0 is booked.")
            return
        
        # set output directory
        self.output = rawdata.me_output
        self.chi2 = results['chi2']
        # ADD DATA
        self.add_data( int(rawdata.n_image), yvalue[0], rawdata.execute_process_in_amp, yerr=yvalue[1] ) #, rawdata.ccd,
                #rawdata.execute_process_in_amp, int(rawdata.n_run), int(rawdata.n_image), datetime_to_days(rawdata.start),
                
        setattr(rawdata,'QT_{}_{}'.format(rawdata.ccd, rawdata.execute_process_in_amp), self.execute_qtest(yvalue[0], rawdata.ccd, rawdata.execute_process_in_amp,'QT:chi2,range',**{'chi2':self.chi2}) )

        return


class MEFitDCLambda(MEAbstract): 
    """Median absolute deviation of the column pixel charge in the overscan as a function of the row
    """

    __name__ = "MEFitDCLambda"
    __sequence_id__ = 63

    def __init__(self):
        super().__init__()
        
        self.is_MECCD = False

        self.title   = "Dark current (from DC fit)"
        self.caption = "This is the best fitted dark current, given by the fit of the pixel charge distribution (PCD) to a set of gaussians convolved by a poisson. The fit is done in the ADUs parameter space to be able to fit also the calibration contant (distance between two consecutive single-electron peaks). The dark current is given in units of electrons per non-binned pixel and per image exposure time, i.e. e-/pix/img."
        self.xlabel  = "image number"
        self.ylabel  = "$\\textit{dark current}, \; \lambda \; [e^{-}/pix/img]$"
        self._ylabel = "Poisson DC [e-/pix/img]"
        
    def execute_process(self,rawdata,**kwargs):

         # Check if the ME has been defined previously
        try:
            results = getattr(rawdata,'me_dcfit_results_{}'.format(rawdata.execute_process_in_amp))
            yvalue  = results['dc']
        except AttributeError:
            self.print_warning("MEFitDC is MANDATORY if MEFitDCLambda is booked.")
            return
        
        # set output directory
        self.output = rawdata.me_output
        self.chi2 = results['chi2']
        # ADD DATA
        self.add_data( int(rawdata.n_image), yvalue[0], rawdata.execute_process_in_amp, yerr=yvalue[1] ) 
        
        # PASS QUALITY TEST TO THE LAST INPUT FILE
        setattr(rawdata,'QT_{}_{}'.format(rawdata.ccd, rawdata.execute_process_in_amp), self.execute_qtest(yvalue[0], rawdata.ccd, rawdata.execute_process_in_amp,'QT:chi2,range',**{'chi2':self.chi2}) )
        return

class MEFitDCCalibration(MEAbstract): 
    """Median absolute deviation of the column pixel charge in the overscan as a function of the row
    """

    __name__ = "MEFitDCCalibration"
    __sequence_id__ = 64

    def __init__(self):
        super().__init__()
        
        self.is_MECCD = False

        self.title   = "Calibration (from DC fit)"
        self.caption = "This is the calibration constant from the fit of the pixel charge distribution (PCD) to a set of gaussians convolved by a poisson. The fit is done in the ADUs parameter space, so the calibration is then the distance between two consecutive single-electron peaks, and is assumed to be the same in the fitted region. The units are ADU/e-."
        self.xlabel  = "image number"
        self.ylabel  = "$\\textit{calibration}, \; k \; [ADU/e^{-}]$"
        self._ylabel = "Calibration [ADU/e-]"



    def execute_process(self,rawdata,**kwargs):

         # Check if the ME has been defined previously
        try:
            results = getattr(rawdata,'me_dcfit_results_{}'.format(rawdata.execute_process_in_amp))
            yvalue  = results['calibration']
        except AttributeError:
            self.print_warning("MEFitDC is MANDATORY if MEFitDCLambda is booked.")
            return
        
        # set output directory
        self.output = rawdata.me_output
        self.chi2 = results['chi2']
        # ADD DATA
        self.add_data( int(rawdata.n_image), yvalue[0], rawdata.execute_process_in_amp,yerr=yvalue[1] ) 

        setattr(rawdata,'QT_{}_{}'.format(rawdata.ccd, rawdata.execute_process_in_amp), self.execute_qtest(yvalue[0], rawdata.ccd, rawdata.execute_process_in_amp,'QT:chi2,range',**{'chi2':self.chi2}) )

        return


##########################################################################################################################################
##
##  images to be included in the Raw Data WEB PAGE
##
##########################################################################################################################################

class MECCDEqualizedImage(MEAbstract):
    """CCD Image, equalized histogram
    """
    __name__= "MECCDEqualizedImage"
    __sequence_id__= 70

    def __init__(self):
        super().__init__()

        # only once per CCD
        self._per_amp = False

        # is not a Run or Image ME
        self.is_MEImage = False
        self.is_MERun = False

        # image can be defined by the user
        self.image = "mean_compressed_pedestal_subtracted"

        self.title   = "Equalized CCD Image "
        self.caption = "Equalized image after pedestal subtraction (see pedestal subtracted in MEOVSPedestalMuPerRow)."
        self.xlabel  = "columns"
        self.ylabel  = "rows"
        self._ylabel = self.ylabel

        self.sigma = 0.1


        self.__units__.update({'image':1, 'sigma':1})
    
    def execute_process(self, rawdata, **kwargs):
        """
        """
        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        image_amp = rawdata.get_image_by_attribute_name(self.image).copy()
        image_amp = image_amp.data

        # make image to be positive
        image_amp = image_amp.data + abs(image_amp.min()) + 1
        image_amp = np.random.normal(np.log10(image_amp), self.sigma)

        # ADD DATA
        self.add_image( int(rawdata.n_image), image_amp, append=False ) 
        
        # plot each single image
        title = "{} [run {}, skips {}, binning {}x{}, VCKDIRN {}]".format(
            rawdata._file_name.replace('_','\_'), rawdata.n_run, rawdata.nskips, rawdata.bin_row, rawdata.bin_col, rawdata.ampdir )

        self.output = "{}file_{}_run{}_".format(rawdata.me_img,rawdata._file_name.split(".")[0],rawdata.n_run)

        # vertical lines to highlight overscan region
        vlines = []
        for amp in rawdata.amplifier.keys():
            vlines.extend(list(rawdata.amplifier[amp]['cols']))

        self.execute_image_plot(image_amp, title, vlines)
         
    def execute_image_plot(self,image,title,vlines,**popts):
        self.palette = "jet"
        super().execute_plot( ptype='imshow', **{'title':title,'vlines':vlines} )

        return

    def  execute_plot(self, ptype='scatter', **popts):
        ### all images are plotted as MECCD, this function runs after the full run is reprocessed where is nothing left to do
        return
        
class MECCDMeanImage(MEAbstract):
    """CCD Image, equalized histogram
    """
    __name__= "MECCDMeanImage"
    __sequence_id__= 71

    def __init__(self):
        super().__init__()

        # only once per CCD
        self._per_amp = False

        # is not a Run or Image ME
        self.is_MEImage = False
        self.is_MERun = False

        # image can be defined by the user
        self.image   = "mean_compressed_pedestal_subtracted"

        self.title   = "Averaged image"
        self.caption = "Averaged skipper image and pedestal subtracted (in units of ADU)."
        self.xlabel  = "columns"
        self.ylabel  = "rows"
        self._ylabel = self.ylabel


        self.__units__.update({'image':1})

    def execute_process(self, rawdata, **kwargs):
        """
        """
        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        if itype is None:
            return
        image_amp = rawdata.get_image_by_attribute_name(self.image)
        image_amp = image_amp.data

        # ADD DATA
        self.add_image( int(rawdata.n_image), image_amp, append=False ) 

        # this is just a ME for the web, will not be included in the report file
           # plot each single image
        title = "{} [run {}, skips {}, binning {}x{}, VCKDIRN {}]".format(
            rawdata._file_name.replace('_','\_'), rawdata.n_run, rawdata.nskips, rawdata.bin_row, rawdata.bin_col, rawdata.ampdir )

        self.output = "{}file_{}_run{}_".format(rawdata.me_img,rawdata._file_name.split(".")[0],rawdata.n_run)

        # vertical lines to highlight overscan region
        vlines = []
        for amp in rawdata.amplifier.keys():
            vlines.extend(list(rawdata.amplifier[amp]['cols']))

        self.execute_image_plot(image_amp, title, vlines)
    
    def execute_image_plot(self,image,title,vlines,**popts):
        super().execute_plot( ptype='imshow', **{'title':title,'vlines':vlines} )

        return

    def  execute_plot(self, ptype='scatter', **popts):
        ### all images are plotted as MECCD, this function runs after the full run is reprocessed where is nothing left to do
        return
        
class MECCDStdImage(MEAbstract):
    """CCD Image, equalized histogram
    """
    __name__= "MECCDStdImage"
    __sequence_id__= 72

    def __init__(self):
        super().__init__()

        # only once per CCD
        self._per_amp = False

        # is not a Run or Image ME
        self.is_MEImage = False
        self.is_MERun = False

        # image can be defined by the user
        self.image = "std_compressed"

        self.title   = "Averaged image"
        self.caption = "Standard deviation skipper image (in units of ADU)"
        self.xlabel  = "columns"
        self.ylabel  = "rows"
        self._ylabel = self.ylabel


        self.__units__.update({'image':1})

    def execute_process(self, rawdata, **kwargs):
        """
        """
        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        if itype is None or rawdata.nskips==1:
            self.failed = True
            return 
        
        image_amp = rawdata.get_image_by_attribute_name(self.image)
        try:
            image_amp = image_amp.data
        except AttributeError:
            raise AttributeError("<MECCDStdImage>: CompressSkipperProcess does not have 'std' in its funciton list!")
        
        # ADD DATA
        self.add_image( int(rawdata.n_image), image_amp, append=False ) 

        # this is just a ME for the web, will not be included in the report file
       # this is just a ME for the web, will not be included in the report file
           # plot each single image
        title = "{} [run {}, skips {}, binning {}x{}, VCKDIRN {}]".format(
            rawdata._file_name.replace('_','\_'), rawdata.n_run, rawdata.nskips, rawdata.bin_row, rawdata.bin_col, rawdata.ampdir )

        self.output = "{}file_{}_run{}_".format(rawdata.me_img,rawdata._file_name.split(".")[0],rawdata.n_run)

        # vertical lines to highlight overscan region
        vlines = []
        for amp in rawdata.amplifier.keys():
            vlines.extend(list(rawdata.amplifier[amp]['cols']))

        self.execute_image_plot(image_amp, title, vlines)
    
    def execute_image_plot(self,image,title,vlines,**popts):
        super().execute_plot( ptype='imshow', **{'title':title,'vlines':vlines} )

        return

    def  execute_plot(self, ptype='scatter', **popts):
        ### all images are plotted as MECCD, this function runs after the full run is reprocessed where is nothing left to do
        return
            
class MECCDCalImage(MEAbstract):
    """CCD Image, equalized histogram
    """
    __name__= "MECCDCalImage"
    __sequence_id__= 73

    def __init__(self):
        super().__init__()

        # only once per CCD
        self._per_amp = False

        # is not a Run or Image ME
        self.is_MEImage = False
        self.is_MERun = False

        # image can be defined by the user
        self.image   = "mean_compressed_pedestal_subtracted_e"

        self.title   = "compressed, pedestal-subtracted calibrated image"
        self.caption = "Averaged, pedestal subtracted and calibrated image."
        self.xlabel  = "columns"
        self.ylabel  = "rows"
        self._ylabel = self.ylabel


        self.__units__.update({'image':1})

    def execute_process(self, rawdata, **kwargs):
        """
        """
        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        if itype is None:
            return
        image_amp = rawdata.get_image_by_attribute_name(self.image)
        image_amp = image_amp.data

        # ADD DATA
        self.add_image( int(rawdata.n_image), image_amp, append=False ) 

        # this is just a ME for the web, will not be included in the report file
           # plot each single image
        title = "{} [run {}, skips {}, binning {}x{}, VCKDIRN {}]".format(
            rawdata._file_name.replace('_','\_'), rawdata.n_run, rawdata.nskips, rawdata.bin_row, rawdata.bin_col, rawdata.ampdir )

        self.output = "{}file_{}_run{}_".format(rawdata.me_img,rawdata._file_name.split(".")[0],rawdata.n_run)

        # vertical lines to highlight overscan region
        vlines = []
        for amp in rawdata.amplifier.keys():
            vlines.extend(list(rawdata.amplifier[amp]['cols']))

        self.execute_image_plot(image_amp, title, vlines)
    
    def execute_image_plot(self,image,title,vlines,**popts):
        super().execute_plot( ptype='imshow', **{'title':title,'vlines':vlines, 'vrange':[-0.5,4.5]} )

        return

    def  execute_plot(self, ptype='scatter', **popts):
        ### all images are plotted as MECCD, this function runs after the full run is reprocessed where is nothing left to do
        return

class MESkewnessCoeff(MEAbstract):
    """CCD Image, equalized histogram
    """
    __name__= "MESkewnessCoeff"
    __sequence_id__= 80

    def __init__(self):
        super().__init__()

        # is not a Run or Image ME
        self.is_MEImage = False
        self.is_MERun = False

        # parameters of the process
        self.id_skip_reference = -1
        self.id_skip_start = 0
        self.id_skip_end = -1
        self.skip_step = 10
        self.kcl_threshold = 3.2
        self.kcl_n_sig = 8

        self.title   = "Skewness coefficient to look for non-conservative skip measurements (charge loss between skips)"
        #self.caption = "Evolution of the averaged pixel charge values in single-skip images as a function of the ordering of the NDCM (i.e. skip number)."
        self.caption = "Skewness coefficient is a parameter to evaluate the cahrge loss between skipper measurements. Values higher than 3 suggest possible charge lost."
        self.xlabel  = "NDCM ordinal"
        self.ylabel  = "$k_{CL}$"
        self._ylabel = "k(CL)"


        self.__units__.update({'id_skip_reference':1,'id_skip_start':1, 'id_skip_end':1, 'skip_step':1, 'kcl_threshold':1, 'kcl_n_sig':8})

    def execute_process(self, rawdata, **kwargs):
        """
        """
        if rawdata.nskips==1:
            self.failed = True
            return 
        
        pkcl = ChargeLossSkewnessProcess()
        # parameters of the process
        pkcl.id_skip_reference = self.id_skip_reference
        pkcl.id_skip_start = self.id_skip_start
        pkcl.id_skip_end = self.id_skip_end
        pkcl.skip_step = self.skip_step
        pkcl.kcl_threshold = self.kcl_threshold
        pkcl.kcl_n_sig = self.kcl_n_sig
        pkcl.display = False

        try:
            pkcl.execute_process(rawdata)
        except RuntimeError:
            # the algorithm has some issues with some images, and not always works .... 
            # Ignore it for now, but this algorithm should be optimized, understand why does not work, and solve it
            self.failed = True
            return
        
        # set output directory
        self.output = rawdata.me_output

        # ADD DATA
        self.add_data( pkcl.skips , pkcl.kcl_coeff, rawdata.execute_process_in_amp )                
        
        # add extra info parameter
        if not hasattr(self,'extrainfo'):
            self.extrainfo = []
        self.extrainfo.append( pkcl.kcl_extrainfo )

        # set ME coming from this process
        setattr(rawdata,'me_{}_kcl_mean_pcdd'.format(rawdata.execute_process_in_amp), pkcl.mean_pcdd )
        setattr(rawdata,'me_{}_kcl_skewness'.format(rawdata.execute_process_in_amp), pkcl.kcl_skew )
        setattr(rawdata,'me_{}_kcl_gauss_fit_mu0'.format(rawdata.execute_process_in_amp), pkcl.kcl_gauss_fit[1] )


class MEDefectsMask(MEAbstract):
    """CCD Image, equalized histogram
    """
    __name__= "MEDefectsMask"
    __sequence_id__= 100

    def __init__(self):
        super().__init__()
        # image can be defined by the user
        self.image = "mean_compressed_pedestal_subtracted"

        ## Qtest
        self.n_image = 0
        self.n_mad = 3

        self.hot_pixel_fraction = 0.4

        self.mask_hot_cols = True
        self.hot_cols_fraction = 0.2

        self.mask_hot_rows = True
        self.hot_rows_fraction = 0.2

        self.q_max_adu = 100.0

        # image can be defined by the user
        self.title = "Masked pixels"
        self.caption = "Mask of the image Hot Pixels"
        self.xlabel  = "columns"
        self.ylabel  = "rows"
        self._ylabel = self.ylabel

        self.__units__.update({'image':1,'n_mad':1,
            'hot_pixel_fraction':1,'hot_cols_fraction':1,'hot_rows_fraction':1,'q_max_adu':1})
    
    def execute_process(self,rawdata,**kwargs):
        """Masks the pixels above the threshold.
        """
        # do not append to the mongo DB
        self.failed = True

        # set MASK ARRAYS FOR THE PRE/OVER-SCAN AND ACTIVE REGIONS related to the proper amplifier
        itype = rawdata.get_image_attribute_name(self.image)
        image_amp = rawdata.get_image_by_attribute_name(self.image)

        # median and mad pixel charge in the overscan region
        ovs_med, ovs_mad = self.median_and_mad(np.ma.array(image_amp,mask=rawdata.mask_image_overscan_cols), axis=1, keepdims=1)
        
        ### mask all pixels above med + n_mad*mad
        #image_mad_mask = rawdata.image_mean_compressed > ovs_med + self.n_mad*ovs_mad
        image_mad_mask = image_amp > ovs_med + self.n_mad*ovs_mad
        
        ### ignore pixels above a given value --- belongs to a cluster
        image_mad_mask[image_amp>self.q_max_adu] = False

        if not hasattr(self, "image_mask"):
            self.image_mask = np.zeros_like(image_amp.data).astype(int)
        
        self.image_mask[~rawdata.mask_image_active_region] = self.image_mask[~rawdata.mask_image_active_region] + image_mad_mask.astype(int)[~rawdata.mask_image_active_region]
        
        self.n_image +=1

        if not hasattr(self,'image_header'):
            self.run    = rawdata.run
            self.output = rawdata.me_output
            self.image_header = rawdata.image_header.copy()
            self.n_amp = rawdata.n_amp
        return

    def execute_plot(self,**kwargs):
        """Calculates the fraction of images a pixel is above the threshold. Plots an image with the
        frecuency the pixel satisfies the condition. Plots an image with those pixels with
        frecuencies above the fraction threshold.
        """
        
        ### No overscan
        if self.n_image == 0:
            return
        # renormalize to the number of amplifiers
        self.n_image = self.n_image / float(self.n_amp)

        if self.hot_pixel_fraction is None:
            self.hot_pixel_fraction=0.5

        ### if a pixel has fired for more than self.hot_pixel_fraction per cent images, is
        #       considered as hot pixel
        self.freq = self.image_mask/float(self.n_image)

        #### hot pixels?
        #       hot pixel is defined as a pixel that appears as masked more than a given frequency
        self.mask_freq = (self.freq > self.hot_pixel_fraction).astype(int)

        #### detecting hot columns to be masked
        row_slice = slice(0, self.image_mask.shape[0])
        col_slice = slice(0, self.image_mask.shape[1])
        img_slice = [row_slice, col_slice]

        if self.mask_hot_rows:
            freq_rows = self.mask_freq.sum(axis=1)/self.image_mask.shape[0]
            hot_rows = np.where(freq_rows > self.hot_rows_fraction)
        print("Hot rows ", hot_rows)

        if self.mask_hot_cols:
            freq_cols = self.mask_freq.sum(axis=0)/self.image_mask.shape[1]
            hot_cols = np.where(freq_cols > self.hot_cols_fraction)[0]
        print("Hot cols:", hot_cols)
        
        for row in hot_rows:
            self.mask_freq[row,:] = 1
        for col in hot_cols:
            self.mask_freq[:,col] = 1

        #for ax_name in ["col","row"]:
        #    if getattr(self, "mask_hot_{}s".format(ax_name)):
        #        axis = self.__axis_id__[ax_name]
        #        axis_to_collapse = int(axis)#int(1-axis)
        #        hot_axis  = img_slice

        #        ### number of masked pixel normalized to the number of pixels along this axis
        #        freq_axis = self.mask_freq.sum(axis=axis_to_collapse)/self.image_mask.shape[axis_to_collapse]

        #        ### if a col or row has more than N masked pixels: mask the full row/col
        #        frac = getattr(self, "hot_{}s_fraction".format(ax_name))
        #        hot_axis[axis] = np.where( freq_axis > frac )[0]

        #        #### masking the full rows or columns
        #        self.mask_freq[hot_axis[0], hot_axis[1]] = 1
        #        print( "        MASKING THE FULL {} ".format(ax_name), hot_axis[axis] , " using frac ", frac)

        # the model is just a constanti
        setattr(self,'fit_model', np.poly1d([self.mask_freq.sum()]))
        setattr(self,'y',[self.mask_freq.sum()])
        setattr(self,'x',[self.run])

        print("         Total number of masked pixels: ", self.mask_freq.sum())

        #################################################################################################
        ##### copy image as fits file
        ###################### ###########################################################################
        self.image_header['run']    = (self.run, 'Run ID number')
        self.image_header['Nfiles'] = (self.n_image, 'Total number of images')
        self.image_header['Fpix']   = (self.hot_pixel_fraction, 'Max frac for a masked pixel')
        if self.mask_hot_rows:
            self.image_header['Fpprow'] = (self.hot_rows_fraction, 'Max frac of masked pixels in one col')
        if self.mask_hot_cols:
            self.image_header['Fppcol'] = (self.hot_cols_fraction, 'Max frac of masked pixels in one row')

        #### masked image not allowed for fits
        if np.ma.isarray( self.freq ):
            self.freq = self.freq.data
        if np.ma.isarray( self.mask_freq ):
            self.mask_freq = self.mask_freq.data

        #self.save_as_fits(
        #        self.output+"../mask/mask_run{:03}.fits".format(int(self.run)),
        #        [self.freq, self.mask_freq],
        #        self.image_header)

        ############# SAVE IMAGE WITH ALL ATRRIBUTES
        ######################################################################################################
        # add file name to the list of files to include variables on the fits file header
        #rawdata.fnames_to_update_header.append(outfile)
        self.SaveAsFits(
            self.output+"../mask/mask_run{:03}.fits".format(int(self.run)),
            [self.freq, self.mask_freq],
            self.image_header)

        #### ADD MASK TO u (SINGLETON Units) TO BE used for other ME
        setattr(u,"run_mask",self.mask_freq)

        return

    def execute_image_plot(self,image,title,**popts):
        super().execute_plot( ptype='imshow', **{'title':title} )

        return

    def SaveAsFits(self,output,hdus,header,dtype=None,overwrite=True):
        """Store a list of HDU data into a fits file
        """

        print("     INFO. Data have been save as fits file: {}\n".format(output))
	### header
        hdu_list = [fits.PrimaryHDU(data=hdus[0], header=header)]

        if len(hdus)>1:
            for data in hdus[1:]:
                if dtype is not None:
                    data = data.astype(dtype)
                hdu_list.append(fits.ImageHDU(data=data.data))
        
        hdul = fits.HDUList(hdu_list)
        hdul.writeto(output,overwrite=overwrite)        


