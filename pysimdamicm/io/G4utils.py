import numpy as np
import re

class G4Volume(object):
    """ 

    multiplicity    Set if the volume mass has some multiplicity, i.e. if the volume is repeated
                    several times, so the total mass, is the mass given for all the reapeated
                    voluems


    return
        
        mass        in kg
        density     in gr/cm^2
        volume      in cm^3
    
    """
    def __init__(self, runtree, PVname ):
        

        runtree.GetEntry(0)
        self.PVname = PVname.replace('\x00','')
        #concatedVolumeNames = str(runtree.concatedVolumeNames).split(';')
        concatedVolumeNames=np.array(runtree.concatedVolumeNames.split(';'))

        #### initializing DM
        self.name=[]
        self.density=[]
        self.mass=[]
        self.volume=[]
        self.surface=[]

        #index = np.where( concatedVolumeNames == self.PVname )[0]
        index = np.where(concatedVolumeNames == self.PVname)[0]

        if len(index) == 1:
            self.index = [index]
        else:
            self.index = []
            ##### assuming the following patter for the collection of simualted volumes
            ##      i.e. KConccd_1, KConccd_2, ... where KConccd is PVname
            patter_vol_name=r'^'+self.PVname+'_'+'[0-9]+'
            for index, vol_name in enumerate(concatedVolumeNames):
                if re.match(patter_vol_name,vol_name):
                    self.index.append( index )
                elif vol_name == self.PVname:
                    self.index.append( index )

            self.index=np.array(self.index) 
        
        if len(self.index)>0:
            print(" -- Found volume ",PVname)
            for ind in self.index:
                ind=int(ind)
                self.name.append(concatedVolumeNames[ind])
                self.density.append(runtree.volumeDensity[ind])
                self.mass.append(runtree.volumeMass[ind])
                self.volume.append(runtree.volumeVolume[ind])
                self.surface.append(runtree.volumeSurface[ind])
        else:
            print(" -- Not found volume ",PVname)
            self.name.append(PVname)
            self.density.append(0.0)
            self.mass.append(0.0)
            self.volume.append(0.0)
            self.surface.append(0.0)

        for attr in ['name','density','mass','volume','surface']:
            setattr(self,attr,np.array(getattr(self,attr)))

    def SetActivity(self, isotope, Texp=61.0, Tcool=183.0, Trun=365.):

        self.isotope = isotope
        self.isRadiogenic = isotope in radiogenic_isotopes_list

        if self.isRadiogenic:
            self.SetRadiogenicActivity(isotope, Texp, Tcool, Trun)
        else:
            self.activity_unit = "decays/Kg/day"
            self.Texp = 0
            self.Trun = 0
            self.Tcool = 0

            try:
                self.activity = isotope_activity_parameters[self.isotope][self.PVname][0]
                self.activity_err = isotope_activity_parameters[self.isotope][self.PVname][1]
            except KeyError:
                mssg = "\t WARNING. Activity for isotope {0} is set to 1.".format(self.isotope)
                mssg += "No values are set for this isotope in volume {0}.".format(self.PVname)
                print( mssg )

                self.activity = 1.0
                self.activity_err = 1.0

        print( "{0} decays/kg/day activity for isotope {1}".format(self.activity, isotope) )



    def SetRadiogenicActivity(self, isotope, Texp=61.0, Tcool=0.0, Trun=0.0):
        """
        isotope ::  name of the isotope as follows, i.e. 60Co will be 60a27z

        Texp    ::  equivalent activation time in units of days
        Tcool   ::  the cooldown time for the simulated volume, in days
        Trun    ::  time since the detector is taking data, in days

        """

        self.Tcool = Tcool
        self.Trun = Trun
        self.Texp = Texp

        umBqkg2decaysPerKgPerday =  0.000001*60*60*24
        S_umBq_kg = isotope_activity_parameters[isotope]["sat"]
        S_decays_perKg_perday = umBqkg2decaysPerKgPerday * S_umBq_kg

        t_half_life = isotope_activity_parameters[isotope]["t_half_life"]

        lambda_iso = np.log(2.0)/t_half_life

        # due to exposure time
        act_Texp  = 1.0 - np.exp( -lambda_iso * Texp )
        # due to cooldown
        act_Tcool = np.exp( -lambda_iso * Tcool )

        self.activity = S_decays_perKg_perday * act_Texp * act_Tcool

        if Trun>0.0:
            # due to running time
            act_Trun  = 1.0/(lambda_iso*Trun) * (1.0 - np.exp(-lambda_iso * Trun) )
            self.activity = act_Trun

        self.activity_unit = "decays/Kg/day"

        return


