import numpy as np
class GeneralizedLJpotetnial:
    """Does stuff
    """
    def __init__(self, StructureContainer=[],A=None,B=None):
        """
        Attributes
        -----------
        hyperparameters : dictionary(hyperparamters)
            dicitonary with relevant hyperparameters for model
        StructureContainer : StructureContainer2Body or StructureContainer3Body(Not yet implemented)
            Structure Container with structures to be featurized
        """
            
        self._StructureContainer=StructureContainer
            
        self._all_features=[]
        self._all_dfeatures=[]
        self._all_sfeatures=[]
        
        Anynone=False
        
        if A ==None:
            self._A = 12
            print("Using A=12 from Lennard Jones")
            Anynone=True
        else:
            self._A = A
        
        if B ==None:
            self._B = 6
            print("Using B=6 from Lennard Jones")
            Anynone=True
        else:
            self._B = B
            
        if Anynone:
            print("CITATION")
            print("See NEED TO PUT IN SOURCE FOR EQUATION for equation details")
        
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
            feature.append(np.array([site['radius']**(-self._A),-site['radius']**(-self._B)]).T)

            tempdfeat = (-1)*np.array([-self._A*site['dradius']*((site['radius'].reshape(-1,1))**((-self._A)-1)),
                           self._B*site['dradius']*((site['radius'].reshape(-1,1))**((-self._B)-1))])

            #tempdfeat.sum(1).T for right feature style
            dfeature.append(tempdfeat.T)
            #.sum(1) for six stress features, xx,yy,zz,xy,yz,xz order
            sfeatures.append(np.vstack([(tempdfeat*-site['relativepos']).T,
                                        (tempdfeat*-np.roll(site['relativepos'], 1,1)).T]))
            
        return feature,dfeature,sfeatures    
    
    
class LJpotetnial(GeneralizedLJpotetnial):
    def __init__(self,StructureContainer=[]):
        super().__init__(StructureContainer,12,6)
        
        
class GeneralizedLJpotetnialWithCutoff:
    """Does stuff
    """
    def __init__(self,cutoff, StructureContainer=[],A=None,B=None):
        """
        Attributes
        -----------
        cutoff : cutoff function that provides a value for each structure
        StructureContainer : StructureContainer2Body or StructureContainer3Body(Not yet implemented)
            Structure Container with structures to be featurized
        """
        self._cutoff=cutoff
        
        self._StructureContainer=StructureContainer
            
        self._all_features=[]
        self._all_dfeatures=[]
        self._all_sfeatures=[]
        
        if A ==None:
            self._A = 12
            print("Using A=12 from Lennard Jones")
            Anynone=True
        else:
            self._A = A
        
        if B ==None:
            self._B = 6
            print("Using B=6 from Lennard Jones")
            Anynone=True
        else:
            self._B = B
            
        if Anynone:
            print("CITATION")
            print("See NEED TO PUT IN SOURCE FOR EQUATION for equation details")
            
            
        
        
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


        for count,structure in enumerate(StructureContainer):
            energy,forces,stress=self._feature_generator(structure,count)
            self._all_features.append(energy)
            self._all_dfeatures.append(forces)
            self._all_sfeatures.append(stress)

    def _feature_generator(self,Structure_NeighborInfo,StructNumber):
        feature=[]
        dfeature=[]
        sfeatures=[]
        for sitenumber,site in enumerate(Structure_NeighborInfo):
            
            cutoff_feats=self._cutoff._all_features[StructNumber][sitenumber]
            cutoff_dfeats=self._cutoff._all_dfeatures[StructNumber][sitenumber]
            
            feature.append(np.array([cutoff_feats.flatten()*((site['radius'])**(-self._A)),cutoff_feats.flatten()*(-(site['radius']**(-self._B)))]).T)

            tempdfeat = (-1)*np.array([(((cutoff_feats*(-self._A)).reshape(-1,1))*site['dradius']*((site['radius'].reshape(-1,1))**((-self._A)-1)))+
                                       cutoff_dfeats*((site['radius'].reshape(-1,1))**(-self._A)),
                                       (cutoff_feats*self._B*site['dradius']*((site['radius'].reshape(-1,1))**((-self._B)-1)))+
                                        cutoff_dfeats*((site['radius'].reshape(-1,1))**(-self._B))])

            #tempdfeat.sum(1).T for right feature style
            dfeature.append(tempdfeat.T)
            #.sum(1) for six stress features, xx,yy,zz,xy,yz,xz order
            sfeatures.append(np.vstack([(tempdfeat*-site['relativepos']).T,
                                        (tempdfeat*-np.roll(site['relativepos'], 1,1)).T]))
            
        return feature,dfeature,sfeatures