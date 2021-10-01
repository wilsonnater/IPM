import numpy as np
class LJpotetnial:
    """Does stuff
    """
    def __init__(self,hyperparameters=None, StructureContainer=[]):
        """
        Attributes
        -----------
        hyperparameters : dictionary(hyperparamters)
            dicitonary with relevant hyperparameters for model
        StructureContainer : StructureContainer2Body or StructureContainer3Body(Not yet implemented)
            Structure Container with structures to be featurized
        """
        
        self._hyperparameters = hyperparameters
        
        if hyperparameters != None:
            print("LJ potential does not use any hypermateres")
            
        self._StructureContainer=StructureContainer
            
        self._all_features=[]
        self._all_dfeatures=[]
        self._all_sfeatures=[]
        
    @property
    def Features(self,Return_Energy=True,Return_Forces=True,Return_Stress=True):
        """ ClusterSpace : copy of the cluster space the structure
        container is based on"""
        
        AllFeatures=[]
        
        EnergyFeatures=[]
        ForceFeatures=[]
        StressFeatures=[]
        for structure_number in range(len(self._all_features)):
            if Return_Energy:
                tempENEfeat=np.array([0.0,0.0])
                
            if Return_Forces:
                tempForcefeat=np.zeros([len(self._all_dfeatures[structure_number]),3,2])
                    
            if Return_Stress:
                tempStressfeat=np.zeros([6,2])
                
            for site in range(len(self._all_features[structure_number])):
                if Return_Energy:
                    tempENEfeat+=self._all_features[structure_number][site].sum(0)
                    
                if Return_Forces:
                    tempForcefeat[site] = self._all_dfeatures[structure_number][site].sum(1)
                
                if Return_Stress:
                    tempStressfeat+=self._all_sfeatures[structure_number][site].sum(1)
                    
            EnergyFeatures.append(0.5*np.array(tempENEfeat))
            ForceFeatures.append(np.array(tempForcefeat))
            StressFeatures.append(np.array(tempStressfeat))

        if Return_Energy:
            AllFeatures.append(EnergyFeatures)

        if Return_Forces:
            AllFeatures.append(ForceFeatures)

        if Return_Stress:
            AllFeatures.append(StressFeatures)


        return AllFeatures

        
    def StructureProcesseror(self,StructureContainer=None):
        if StructureContainer==None:
            StructureContainer = self._StructureContainer


        for structure in StructureContainer:
            energy,forces,stress=self._feature_generator(structure)
            self._all_features.append(energy)
            self._all_dfeatures.append(forces)
            self._all_sfeatures.append(stress)

    def _feature_generator(self,Structure_NeighborInfo):
        feature=[]
        dfeature=[]
        sfeatures=[]
        for site in Structure_NeighborInfo:
            feature.append(np.array([site['radius']**(-12),-site['radius']**(-6)]).T)

            tempdfeat = (-1)*np.array([-12*site['dradius']*((site['radius'].reshape(-1,1))**(-13)),
                           6*site['dradius']*((site['radius'].reshape(-1,1))**(-7))])

            #tempdfeat.sum(1).T for right feature style
            dfeature.append(tempdfeat.T)
            #.sum(1) for six stress features, xx,yy,zz,xy,yz,xz order
            sfeatures.append(np.vstack([(tempdfeat*-site['relativepos']).T,
                                        (tempdfeat*-np.roll(site['relativepos'], 1,1)).T]))
            
        return feature,dfeature,sfeatures
