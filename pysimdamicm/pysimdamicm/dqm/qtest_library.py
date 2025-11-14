import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import datetime

DISPLAY = False 
DISPLAY_PLOTS = False

#DISPLAY = True#False 
#DISPLAY_PLOTS = True#False

###############################################################################################
#####       Qtests Methods 
#####           
###############################################################################################

#def qtest_compareValue(val, ref):
#    """ Compares if the value is greater than the expected value.
#        If it's the expected value returns True.
#    """
#    return val < ref


def qtest_chi2FromFit(x, y, ref, amp, n_std=3, cut_val=2, is_image=False):
    """ Only for ME obtained from a fit. Looks if the chi2 is above the threshold.
    """
    
    Result = y < cut_val
    if DISPLAY:
        print("Chi2 of the Fit is below {}?: ".format(cut_val), Result)
        print({'name':'QTChi2FromFit','result':Result, 'chi2':y, 'run_ref':None})
    return {'name':'QTChi2FromFit','result':Result, 'chi2':y, 'run_ref':None}


def qtest_chi2(x, y, ref, amp, n_std = 3, cut_val = 2, is_image = False):
    """ Uses the chi2 test to see if the distribution came
        from the reference model.
    """
    
    ### Loads the refence values
    #x_ref,y_ref = np.array(ref["x"]), np.array(ref["y"])
    x_ref,y_ref = np.array([np.array(x_i) for x_i in np.array(ref["x"])[np.array(ref["amplifier"])==amp]]), np.array([np.array(y_i) for y_i in np.array(ref["y"])[np.array(ref["amplifier"])==amp]])

    if is_image:
        ### Applies chi2-test to each row/col point compared to the distribution of points per row/col
        #XXX This only applies for ME that has only one value per image. If this is the case
        #It returns None
        if type(y) is list:
            ## Averages the columns or rows values of the images
            y_mean_ref = np.median(y_ref, axis = 0)
            x_mean_ref = np.median(x_ref, axis = 0)
            y_std_ref = np.std(y_ref, axis = 0)
            ### Masks those rows/cols above n_std and fits the averaged rows/cols
            mask_high_ref = np.abs(y_mean_ref-np.mean(y_mean_ref))<n_std*np.std(y_mean_ref)
            m_ref, y0_ref = np.polyfit(x_mean_ref[mask_high_ref],y_mean_ref[mask_high_ref],1)
            ### Calculates chi2
            chi2 = np.sum(np.power((y-m_ref*np.array(x)-y0_ref)/y_std_ref,2))
            if DISPLAY_PLOTS:
                m,y0 = np.polyfit(x,y,1)
                plt.plot(x,y, 'o',color = 'k')
                plt.plot(x,m_ref*np.array(x)+y0_ref,'--',color = "r")
                plt.plot(x,m*np.array(x)+y0,'--',color="b")
                plt.show()
            ### Calculates the reduced chi2/NdF
            chi2_ndf = chi2/(len(x)-2) 
        else:
            ### Dont apply to image ME
            return {'name':'QTChi2','result':np.nan, 'chi2_ndf':np.nan, 'run_ref':ref['run']} 
    else:
        ### Applies chi2-test to the whole RUN
        try:
            if x_ref.shape[1] > 0:
                ### Chi2test over rows/cols ME
                y_mean_ref = np.median(y_ref, axis = 0)
                x_mean_ref = np.median(x_ref, axis = 0)
                y_std_ref = np.std(y_ref, axis = 0)
                ### Calculates the reference model
                mask_high_ref = np.abs(y_mean_ref-np.mean(y_mean_ref))<n_std*np.std(y_mean_ref)
                m_ref, y0_ref = np.polyfit(x_mean_ref[mask_high_ref],y_mean_ref[mask_high_ref],1)
                ### Calculates chi2 and chi2/NdF
                chi2 = np.sum(np.power((np.mean(y, axis = 0) - m_ref*np.array(x)-y0_ref)/y_std_ref,2))
                chi2_ndf = chi2/(len(x)-2)
        except IndexError:
            ### Chi2test over image ME
            y_std_ref = np.std(y_ref)
            mask_high_ref = np.abs(y_ref-np.mean(y_ref))<n_std*y_std_ref
            m_ref, y0_ref = np.polyfit(x_ref[mask_high_ref],y_ref[mask_high_ref],1)
            chi2 = np.sum(np.power((y - m_ref*np.array(x)-y0_ref)/y_std_ref,2))
            chi2_ndf = chi2/(len(x)-2)
    ### Threshold to confirm that chi2 / ndf is bad 
    Result = (chi2_ndf < cut_val) #2)  

    if DISPLAY:
        print('Chi2 / ndf: ', chi2_ndf)
        print({'name':'QTChi2','result':Result, 'p_value':chi2_ndf, 'run_ref':ref['run']})
    return {'name':'QTChi2','result':Result, 'chi2_ndf':chi2_ndf, 'run_ref':ref['run']}


def qtest_ks(x, y, ref, amp, n_std = 3, cut_val = 5e-3, is_image = False):
    """ Uses the kolmogorvo-smirnov test to see if the distribution came
        from the reference.
    """
    
    ### Loads the refence values
    #x_ref,y_ref = np.array(ref["x"]), np.array(ref["y"])
    x_ref,y_ref = np.array(ref["x"])[np.array(ref["amplifier"])==amp], np.array(ref["y"])[np.array(ref["amplifier"])==amp]

    if is_image:
        ### Applies ks-test to each row/col point compared to the distribution of points per row/col
        #XXX This only applies for ME that has only one value per image. If this is the case
        #It returns None
        if type(y) is list:
            ### Compares each image of reference with the current image
            ks = []
            ## Calculates ks test of each reference image with the current image
            for y_i in y_ref:
                ks.append(stats.ks_2samp(y,y_i)[1])
            ## Averages the p-values of all the images
            p_value = np.mean(ks)
        else:
            ### Dont apply to image ME
            return {'name':'QTKstest','result':np.nan, 'p_value':np.nan, 'run_ref':ref['run']} 
    else:
        ### Applies ks-test to the whole RUN
        try:
            if x_ref.shape[1] > 0:
                ### Ttest over rows/cols ME
                ### Compares run images with the reference images and obtains the pvalue
                p_value = stats.ks_2samp(np.mean(y, axis = 0), np.mean(y_ref, axis = 0))[1]
        except IndexError:
            ### Kstest over image ME
            p_value = stats.ks_2samp(y, y_ref)[1]

    ### Threshold to confirm that p-value is 0 (95 confidence level)
    Result = p_value > cut_val #5e-3

    if DISPLAY:
        print('P-Value KS-est: ', p_value)
        print({'name':'QTKstest','result':Result, 'p_value':p_value, 'run_ref':ref['run']})
    return {'name':'QTKstest','result':Result, 'p_value':p_value, 'run_ref':ref['run']}


def qtest_tt(x, y, ref, amp, n_std = 3, cut_val = 5e-3, is_image = False):
    """ Uses the ttest test (Welchs version for different std) to see if the distribution came
        from the reference.
    """
    ### Loads the refence values
    #x_ref,y_ref = np.array(ref["x"]), np.array(ref["y"])
    x_ref,y_ref = np.array(ref["x"])[np.array(ref["amplifier"])==amp], np.array(ref["y"])[np.array(ref["amplifier"])==amp]
    
    if is_image:
        ### Applies t-test to each row/col point compared to the distribution of points per row/col
        #XXX This only applies for ME that has only one value per image. If this is the case
        #It returns None
        if type(y) is list:
            ### Compares each ref image with the current image
            #tt = []
            #for y_i in y_ref:
            #    tt.append(stats.ttest_ind(y,y_i, equal_var=False)[1])
 
            tt = stats.ttest_ind(y,y_ref, equal_var=False)[1]
            ### Averages the obtained row/col p-values
            p_value = np.mean(tt)
            print("Ttest p-value: ", p_value)
        else:
            ### Dont apply to image ME
            return {'name':'QTTtest','result':np.nan, 'p_value':np.nan, 'run_ref':ref['run']} 
    else:
        ### Applies t-test to the whole RUN
        try:
            if x_ref.shape[1] > 0:
                ### Ttest over rows/cols ME
                p_value = stats.ttest_ind(np.mean(y, axis = 0), np.mean(y_ref, axis = 0), equal_var=False)[1]
        except IndexError:
            ### Ttest over image ME
            p_value = stats.ttest_ind(y, y_ref, equal_var=False)[1]

    ### Threshold to confirm that p-value is 0
    Result = p_value > cut_val #5e-3

    if DISPLAY:
        print('P-Value TTest: ', p_value)
        print({'name':'QTTtest','result':Result, 'p_value':p_value, 'run_ref':ref['run']})
    return {'name':'QTTtest','result':Result, 'p_value':p_value, 'run_ref':ref['run']}


def qtest_residuals_per_image(x ,y, ref, amp, n_std = 3, cut_val = 5, is_image = False):
    """ Calculates the residuals for those ME that are plotted against the image number.
        The residuals are calculated by subtracting the linear fit to the data.
        Gets the number of points outside mean +- n_std*std and compares them against the expected value.
    """
    
    if type(x) is not list:
        x = [x]
        y = [y]

    # REFERECE DATA
    #x_ref,y_ref = np.array(ref["x"]), np.array(ref["y"])
    x_ref,y_ref = np.array(ref["x"])[np.array(ref["amplifier"])==amp], np.array(ref["y"])[np.array(ref["amplifier"])==amp]

    if datetime.datetime == type(ref["x"][0]):
        ### If x value is date from start 
        y_ref = np.array(y_ref.tolist())
        x_ref = np.array([i+1 for i in range(len(y_ref))])
    ## Masks pixels above 3sigma to fit the reference data to a linear model
    mask_high_ref = np.abs(y_ref-np.mean(y_ref))<=n_std*np.std(y_ref)
    m_ref, y0_ref = np.polyfit(x_ref[mask_high_ref],y_ref[mask_high_ref],1)
    ### Calculates the residuals reference and the residuals of the data with reference

    if is_image:
        ### There will be only one point so we cant subtract the model
        y_ref_res = np.array(y_ref) - (m_ref*x_ref)
        y_data_res = np.array(y)
    else:    
        ### Subtract the fitted model to data and ref 
        y_ref_res = np.array(y_ref) - (m_ref*x_ref)
        y_data_res = np.array(y)-(m_ref*np.array(x)+y0_ref)

    ## Calculates Standard deviation of reference
    std_ref = np.std(y_ref_res)
    
    #### NURIA ####
    #  - Mask pixels above 3sigma to fit the reference data to a linear model
    # FIXME we really need to remove the outliers?
    #mask_high_ref = np.abs(y_ref-np.mean(y_ref))<=3*np.std(y_ref)
    #_model_pars = np.polyfit(x_ref[mask_high_ref],y_ref[mask_high_ref],1)
    #model = np.poly1d(_model_pars)
    
    ### Calculates the residuals reference and the residuals of the data with reference 
    #y_ref_res = np.array(y_ref) - model(x_ref)
    
    # CALCULATE THE RESIDUALS AGAINTS THE REFERENCE MODEL
    # FIXME the tag is_image is misleading because not all ME for image has a single value, this
    # should be something like is_point
    #if is_image:
    #    ### There will be only one point so we cant do the fit
    #    y_data_res = np.array(y)
    #else:    
    ### Fit the data
    #mask_high = np.abs(y-np.median(y))<=3*np.std(y)
    #y_data_res = np.array(y) - model(x)

    ## Calculates Standard deviation of reference
    #std_ref = np.std(y_ref_res)
    ###############
    mean_ref = np.mean(y_ref_res)

    ## Contains the images that are above the threshold
    images_above = np.abs(y_data_res-mean_ref)>n_std*std_ref

    if DISPLAY:
        print('Check this images: ', np.array(x)[images_above])
        print('Number of sigmas above: ', np.round(y_data_res[images_above]/std_ref,0))
    
    ### x and y value of the failed images
    x_fail = np.array(x)[images_above].tolist()
    y_fail = np.array(y)[images_above].tolist()

    if DISPLAY:
        print("NUMBER OF FAILS :", len(y_fail))

    ### If all the y values passed the test stores the result as None
    if len(x_fail) == 0:
        x_fail = None
        y_fail = None

    ### Mean and std of the data
    mean_data, std_data = np.mean(y_data_res), np.std(y_data_res)

    ### Calculates the fraction of rows or columns above the expected value
    if is_image:
        ### If the image is bad np.sum > 0
        Result = np.sum(images_above) == 0
        if DISPLAY:
            print("Image is out of sigma region")
    else:
        ### Apply Qtest to all the RUN
        ### We accept only a certan percentage of bad images
        Result = np.sum(images_above)*100/len(x) <= cut_val #less than 5 per cent of row/cols bad

        #### Check if median is displaced.
        std_r = np.std(y_ref_res)
        if (mean_data-std_data > mean_ref+std_r) or (mean_data+std_data < mean_ref-std_r):
            ## If the median is displaced the given number of std the Qtest fails
            Result = False
            print("The mean is displaced!!!!")
        else:
            ## If not it stores the original result
            Result = Result

    ### Calculates the fraction of images above the expected value   
    if DISPLAY:
        print("Residuals Passed: ", Result)   
    if DISPLAY_PLOTS:
        ### Plot points
        plt.clf()
        plt.scatter(np.array(x_ref)-np.min(x_ref),y_ref_res,color="g", alpha=0.7, label="Ref")
        plt.scatter(np.array(x)-np.min(x),y_data_res,color="k", alpha=0.7, label="Data")
        plt.scatter((np.array(x)-np.min(x))[images_above], np.array(y_data_res)[images_above], color = "r", alpha = 0.7, label='Failed')
        plt.scatter(np.array(x_ref)-np.min(x_ref),y_ref_res,color="g", alpha=0.7, label="Ref")
        ### Plot limits
        x_max_lim = np.max(x_ref) if np.max(x)<=np.max(x_ref) else np.max(x)
        x_min_lim = np.min(x) if np.min(x)<=np.min(x_ref) else np.min(x_ref)
        plt.hlines(mean_ref+n_std*std_ref, 0, x_max_lim-x_min_lim, linestyle='dashed', color = "g", label="{}Sigma".format(n_std))
        plt.hlines(mean_ref-n_std*std_ref, 0, x_max_lim-x_min_lim, linestyle='dashed', color = "g")
        plt.hlines(mean_ref, 0, x_max_lim-x_min_lim, linestyle='solid', color = "g")
        ### Plot labels
        #plt.legend(loc="best")    
        plt.xlabel(ref["xlabel"])
        plt.ylabel(ref["ylabel"])
        plt.title(ref["title"])
        plt.show()
        #plt.clf()

    if DISPLAY:
        print( {'name':'QTresiduals','result':Result, 'x':x_fail, 'y':y_fail, 'mean_ref':mean_ref,'n_std':n_std, 'std_ref':std_ref, 'run_ref':ref['run']} )

    return {'name':'QTresiduals','result':Result, 'x':x_fail, 'y':y_fail, 'mean_ref':mean_ref,'n_std':n_std, 'std_ref':std_ref, 'run_ref':ref['run']}


def qtest_residuals_per_rowCol(x, y, ref, amp, n_std = 3, cut_val = 5, is_image = False):
    """ Calculates the residuals for those ME that are plotted against the row or column number.
        All the individual columns or rows are averaged for every the images.
        The residuals are calculated by subtracting a linear fit per image.
        Gets the number of points outside mean +- n_std*std and compares them against the expected value.
    """
    
    ### numpy arrays
    x = np.array(x)
    y = np.array(y)

    ### Loads the refence values
    #x_ref,y_ref = np.array(ref["x"]), np.array(ref["y"])
    x_ref,y_ref = np.array(ref["x"])[np.array(ref["amplifier"])==amp], np.array(ref["y"])[np.array(ref["amplifier"])==amp]

    ## Averages the columns or rows values of the images
    y_mean_ref = np.median(y_ref, axis = 0)
    x_mean_ref = np.median(x_ref, axis = 0)
    ### Masks those rows/cols above n_std and fits the averaged rows/cols
    mask_high_ref = np.abs(y_mean_ref-np.mean(y_mean_ref))<n_std*np.std(y_mean_ref)
    #_model_pars = np.polyfit(x_mean_ref[mask_high_ref],y_mean_ref[mask_high_ref],1)
    #model = np.poly1d(_model_pars)
    m_ref, y0_ref = np.polyfit(x_mean_ref[mask_high_ref],y_mean_ref[mask_high_ref],1) 
    
    ### Calculates the residuals of the reference   
    y_ref_full_res = np.array(y_ref) - (m_ref*x_ref+y0_ref)
    y_ref_res = np.mean( np.array(y_ref) - (m_ref*x_ref+y0_ref), axis = 0 )
    ## Computes the median colums/rows values of the images.
    if is_image:
        y_mean = y
        x_mean = x
    else:    
        y_mean = np.median(y, axis = 0)
        x_mean = np.median(x, axis = 0)
    
    ## Masks rows/cols above n_std
    mask_high = np.abs(y_mean-np.median(y_mean))<3*np.std(y_mean)
    ## Fit the distribution
    m_mean,y0_mean = np.polyfit(x_mean[mask_high],y_mean[mask_high],1)

    ### Obtains the std of each image and then calculates the average std
    std_ref = np.mean(np.std(y_ref_res,axis = 0))

    if is_image:
        ### Calculates each row/col std and averages all the row/col std
        std_ref = np.mean(np.std(y_ref_full_res,axis = 0))
    mean_ref = np.mean(np.mean(y_ref_res,axis = 0))

    ## Calculates the residuals of averaged data and fit
    if is_image:
        y_data_res =  np.array(y)-(m_ref*np.array(x)+y0_ref)
    else:
        ### Calculates each image row/col residuals and averages to one value per RUN
        y_data_res = np.mean( np.array(y)-(m_ref*np.array(x)+y0_ref), axis = 0 )
    
    ### Calculates the rows/cols that are above the threshold in the RUN
    images_above = np.abs(y_data_res-mean_ref)>n_std*std_ref

    if DISPLAY:
        print('Check Rows/Cols (Average value high): ', np.array(x_mean)[images_above])
        print('Number of sigmas above: ', np.round(y_data_res[images_above]/std_ref,0))

    ### Failed positions
    x_fail = np.array(x_mean)[images_above].tolist()
    y_fail = np.array(y_mean)[images_above].tolist()

    if len(x_fail) == 0:
        ### If all the cols/rows passed the test results are stored as None
        x_fail, y_fail = None, None

    ### Median RUN and Ref
    mean_data, std_data = np.mean(y_data_res), np.std(y_data_res)
    ### Calculates the fraction of rows or columns above the expected value
    ## Only allows to have a certain percentage of bad images
    if is_image:
        Result = np.sum(images_above)*100/len(x) <= cut_val#5 #less than 5 per cent of row/cols bad
        if DISPLAY:
            print("Number of images: ", len(x))
        print("Percentage of Bad images: ", np.sum(images_above)*100/len(x))
    else:
        Result = np.sum(images_above)*100/len(x) <= cut_val#5 #less than 5 per cent of row/cols bad

    #### Check if mean is displaced.
    std_r = np.std(y_ref_res)
    if (mean_data-std_data > mean_ref+std_r) or (mean_data+std_data < mean_ref-std_r):
        ## If the mean is displaced the Qtest fails
        Result = False
        print("The mean is displaced!!!!")
    else:
        ## If not stores the original value
        Result = Result

    if DISPLAY:
        print("Residuals Passed: ", Result)   

    if DISPLAY_PLOTS:
        print('Run Mean above threshold: ',np.abs(np.mean(y_data_res)-np.mean(y_ref_res))>n_std*np.std(y_ref_res) )
        plt.clf()
        if not is_image:
            #plt.plot(x_mean,y_mean,'o',color="b",alpha=0.7,label="Data")
            #plt.plot(x_mean,m_mean*np.array(x)+y0_mean,'-',color="b",alpha=0.7,label="Fit")
            plt.scatter(x_ref[0],y_ref_res,color="g", alpha=0.7, label="Ref")
            plt.scatter(x_mean,y_data_res,color="k", alpha=0.7, label="Data")
            plt.scatter(np.array(x_mean)[images_above], np.array(y_data_res)[images_above], color = "r", alpha = 0.7, label='Failed')
            plt.hlines(np.mean(y_ref_res)+n_std*std_ref,x[0][0], x[0][-1], linestyle='dashed', color = "g", label="{}Sigma".format(n_std))
            plt.hlines(np.mean(y_ref_res)-n_std*std_ref,x[0][0], x[0][-1], linestyle='dashed', color = "g")
            plt.hlines(np.mean(y_data_res),x[0][0], x[0][0], linestyle='dashed', color = "k")
            plt.hlines(np.mean(y_ref_res),x[0][0], x[0][-1], color = "g")
        else:
            plt.scatter(x_mean,y_data_res,color="k", alpha=0.7, label="Data")
            plt.scatter(np.array(x_mean)[images_above], np.array(y_data_res)[images_above], color = "r", alpha = 0.7, label='Failed')
            plt.scatter(x_ref[0],y_ref_res,color="g", alpha=0.7, label="Ref")
            plt.hlines(np.mean(y_ref_res)+n_std*std_ref,x[0], x[-1], linestyle='dashed', color = "g", label="{}Sigma".format(n_std))
            plt.hlines(np.mean(y_ref_res)-n_std*std_ref,x[0], x[-1], linestyle='dashed', color = "g")
            plt.hlines(np.mean(y_data_res),x[0], x[0], linestyle='dashed', color = "k")
            plt.hlines(np.mean(y_ref_res),x[0], x[-1], color = "g")

        plt.legend(loc="best")    
        plt.xlabel(ref["xlabel"])
        plt.ylabel(ref["ylabel"])
        plt.title(ref["title"])
        plt.show()
    if DISPLAY:
        print({'Result':Result, 'x_fail':x_fail, 'y_fail':y_fail, 'mean_ref':mean_ref, 'n_std':n_std, 'std_ref':std_ref})

    return {'name':'QTresiduals','result':Result, 'x':x_fail, 'y':y_fail, 'mean_ref':mean_ref,'n_std':n_std, 'std_ref':std_ref, 'run_ref':ref['run']}
