
from pysimdamicm.utils.units import Units
u=Units()
from pysimdamicm.dqm.dqm_manager import DQMManager
from pysimdamicm.scripts import anaSKImg

from matplotlib import use,get_backend
use(get_backend())

import time
from astropy.io import fits
from matplotlib import pyplot as plt
from glob import glob 
import os

def raise_error(msm,EType):
    raise EType("\x1b[35m {}\n\x1b[m".format(msm))

def main(run,infile_list,config,
        me_reference=None,
        substring_skip_list=['_full_1.fits'],
        image_HR='skip_1.fits',
        dqm_data_dir=None,dqm_out_dir=None,display=False,all_me=False,
        dopdf=True, image_id=None, ccd_id=None, run_tag=None,nopdf=True):
    
    # If called by panaSKImg this two parameters can be None
    if substring_skip_list is None:
        substring_skip_list=[]

    print("DQM starting ... ")

    #DQM mode. In this mode several things are required
    #1. DQM_OUT_DIR: output where to create AN STRUCTURE OF DIRECTORIES FOR THE OUTPUTS
    #####################################################################################
    dqm_data_dir_env = os.environ.get('DQM_DATA_DIR')
    if dqm_data_dir_env is not None:
        dqm_data_dir = dqm_data_dir_env
    if (dqm_data_dir is None):
        dqm_data_dir = os.environ.get('DQM_DATA_DIR')
        if dqm_data_dir is None:
            msm ="Input data directory is mandatory: use --dqm-data-dir or "
            msm+="define the environment variable 'DQM_DATA_DIR' pointing to your data"
            raise_error(msm,IOError)
    # Is directory?
    dqm_data_dir = os.path.abspath(dqm_data_dir)
    if not os.path.isdir(dqm_data_dir):
        raise_error("{} is not a directory".format(dqm_data_dir),ValueError)
       
    #2. DQM_DATA_DIR:   where data is and a patter file name
    #####################################################################################
    if dqm_out_dir is None:
        dqm_out_dir = os.environ.get('DQM_OUT_DIR')
        if dqm_out_dir is None:
            msm ="Directory for the output is mandatory for the DQM running mode: use --outdir or "
            msm+="define the environment variable 'DQM_OUT_DIR' pointing where reprocessed"
            msm+=" data should be stored"
            raise_error(msm,IOError)
    # Is directory?
    dqm_out_dir = os.path.abspath(dqm_out_dir)
    if not os.path.isdir(dqm_out_dir):
        raise_error("{} is not a directory".format(dqm_out_dir),ValueError)
    
    #3. DQM_RUN:  number of the run to be processed
    #####################################################################################
    if run is None:
        run = os.environ.get('DQM_RUN')
        if run is None:
            msm = "Run number ID is mandatory for during DQM running mode: use --run or "
            msm+= "define the environment variable 'DQM_RUN' pointing to the run number"
            raise_error(msm,IOError)
       

    # CREATE DIRECTORY STRUCTURE FOR THE OUTPUTS
    #####################################################################################
    output = {
            ### root directory for the output structure
            'out_path': dqm_out_dir,
            ### run number
            'run':"{:3}".format(int(run)),
            ###
            'run_path':"{}/run{:03}".format(dqm_out_dir,int(run))
            }
    # create dictionary with all sub-directories
    output.update({
        'avg':"{}/avgimg".format(output['run_path']),
        'me':"{}/me".format(output['run_path']),
        'me_pcd':"{}/me/pcd".format(output['run_path']),
        'me_ect':"{}/me/ect".format(output['run_path']),
        'me_img':"{}/me/img".format(output['run_path']),
        'me_dcfit':"{}/me/dcfit".format(output['run_path']),
        'me_sed':"{}/me/sed".format(output['run_path']),
        'me_recon':"{}/recon".format(output['run_path']),
        'me_report':"{}/me/report".format(output['run_path']),
        'mask':"{}/mask".format(output['run_path']),
        'logs':"{}/logs".format(output['run_path']),
        'others':"{}/others".format(output['run_path'])
        })
    # create structure of directories under the root directory `out_path`
    for sub_dir in ['run_path','avg','me','logs','others','me_pcd','me_ect','me_img','me_dcfit','me_sed','me_recon',
            'mask','me_report']:
        subfolder = output[sub_dir]
        if not os.path.exists(subfolder):
            os.mkdir(subfolder)
        
    # IMAGES FOUND IN THE DQM_DATA_DIR
    run_images = sorted(glob("{}/{}".format(dqm_data_dir,infile_list[0])))
    for pattern in infile_list[1:]:
        run_images.extend(glob("{}/{}".format(dqm_data_dir,pattern)))
    # XXX Assuming patterns lile /.../.../string_some_NNN_XX.fits
    #   order the list with string XX, assuming is a number and the image number
    #run_images = sorted(run_images, key=lambda f: int(f.split("/")[-1].split(".")[0].split("_")[-1]) )

    # filtering images if skip has been used
    if substring_skip_list is not None:
        ### looks for a file type input and appends internal files to substring list for removal
        for s in substring_skip_list:
            if s.endswith('.txt'):
                cuts = open(s, 'r')
                for i in cuts.readlines():
                    substring_skip_list.append(i.strip())

        selected_images = []
        ignored_images = []
        for image in run_images:
            if any(map(image.__contains__, substring_skip_list)):
                ignored_images.append(image)
            else:
                selected_images.append(image)
        if len(ignored_images)>0:
            print("INFO. {} have been ignored: \n".format(len(ignored_images),ignored_images))
    else:
        selected_images = run_images

    # FIND High RESOLUTION IMAGE, IF exists
    #found_image_HR = False
    #if image_HR is not None:
    #    image_HR_in_run = next((s for s in selected_images if image_HR in s),None)
    #    if image_HR_in_run is not None:
    #        selected_images.remove(image_HR_in_run)
    #        ### Move image_HR to be the first to process: 
    #        #       the gain is needed for the rest of the images
    #        selected_images.insert(0,image_HR_in_run)
    #        found_image_HR = True
        
    #################################################################################
    ### dqm: DATA QUALITY MONITOR
    #################################################################################
    # ADD ALL ME BY DEFAULT
    dqm_man = DQMManager(run)
    #for attr in dir(dqm_man):
    #    if attr.count('active')>0:
    #        setattr(dqm_man,attr,True)
    dqm_man.set_configuration(config[1])
    # ADD REFERENCE IF ANY
    if me_reference is not None:
        dqm_man.me_summary_ref = me_reference

    # RUNNING SEQUENTIALLY OVER ALL IMAGES
    ##########################################################################################
    t_init = time.time()
    Nimgs = len(selected_images)
    if Nimgs == 0:
        raise RuntimeError("<DQM>: Number of images is 0: DQM killed!")

    print("DQM INFO: Processing Run number {} with {} images".format(run,Nimgs))
    for n_id,image in enumerate(selected_images):
        if not os.path.isfile(image):
            Nimgs -=1
            continue

        rawdata = anaSKImg.main([image],config[0], outdir=output, display=display, run_id=run,
                image_id=image_id, ccd_id=ccd_id, run_tag=run_tag)
        plt.close('all')
       
        # if image number is 0, change it to be the ID from the list
        if rawdata.n_image==0:
            rawdata.n_image=int(n_id)
        
        # execute all ME
        print("DQM INFO: Running all ME for this image")
        dqm_man.execute_process(rawdata)
        
        # save all MEs as pickle file
        pkl_name = dqm_man.save_MEs_as_dataframe(output['me_report'],rawdata._file_name.split(".")[0])
        
        ### WRITE KEYWORDS IN THE FITS FILE HEADER (IF ANY)
        # panaSKImg generates a multiextension fits file, add it to update header 
        if hasattr(rawdata,'_multiext_waders_image'):
            rawdata.fnames_to_update_header.append(rawdata._multiext_waders_image)
        for fname in rawdata.fnames_to_update_header:
            for pkey in u._tohdr.keys():
                val,txt = u._tohdr[pkey]
                fits.setval(fname, pkey, value=val, ext=rawdata.extension, comment=txt)
            
    #### Abstract for the DQM report 
    dqm_report_abstract = [
            r"\noindent","\n",
            "Data directory:","\n",
            r"\verb|{}|".format(dqm_data_dir),"\n\n",
            "Output directory:","\n",
            r"\verb|{}|".format(dqm_out_dir),"\n\n",
            "Reference used:","\n",
            r"\verb|{}|".format(me_reference),"\n\n",
            "Total images: {} ".format(Nimgs),
            ]

    
    #### Build summary ME pkl object and REPORT
    #dqm_man.dump_run_summary(run,output['me'])

    print("DQM INFO: Create DQM report")
    if nopdf:
        print("DQM INFO: No PDF report is being created")
    else:
        ### if MEDCFit was booked, get the gain
        print("DQM INFO: Creating ME plots: ")
        dqm_man.execute_plot()

        if dopdf:
            dqm_man.build_report(run,output['me_report'],abstract=dqm_report_abstract,dopdf=dopdf)
        else:
            dqm_man.build_report(run,output['me_report'],abstract=dqm_report_abstract,dopdf=True,file_name=rawdata._file_name.replace('.fits',''))
    
    # save mongoDB document
    # dqm_man.close_DB_document(run,output['me_report'])

    t_end = time.time()
    print("DQM INFO: Total number of reprocessed images ", Nimgs)
    print("DQM INFO: Total execution time: {} s".format(t_end-t_init))

    print("DQM INFO: Completed Successfully!")
    print("DQM PROCESS STATUS::OK, RUN::{}, PKL::{}".format(run,pkl_name.split('/')[-1]))
 
    return

