import numpy as np

class Structure2BodyContainer:
    """
    This class serves as a container for structures and the
    2 body parameter (radius and local distance)

    Parameters
    -----------
    cutoff_radius : float or int
        cutoff radius for checking for neighbors

    structurelist : list(pymatgen strucutres)
         structures to be added to the container
    """
    def __init__(self,cutoff_radius, structure_list=None):
        """
        Attributes
        -----------
        _structure_list : list(pymatgen strucutres)
            structures to add to container
        _cutoff_radius : float
            cutoff distance it will look for nearest neighbors
        """

        self._cutoff_radius = cutoff_radius
        
        # Add atoms from atoms_list
        self._featurized_structure_list = []
        if structure_list is not None:
            self.Add_Structure_list(structure_list)
                
                
        '''
        I think this should check if structure is structure type
        also should be processing structures
        if _structure_list is not None:
            for fit_structure in _featurized_structure_list:
                if not isinstance(fit_structure, FitStructure):
                    raise TypeError('Can only add FitStructures')
                self._structure_list.append(fit_structure)
        '''
    def __len__(self):
        return len(self._featurized_structure_list)
    
    def __getitem__(self, ind):
        return self._featurized_structure_list[ind]
    
    @property
    def featurized_structure_list(self):
        """ _featurized_structure_list: copy of the cluster space the structure
        container is based on"""
        return self._featurized_structure_list

    @property
    def cutoff_radius(self):
        """ _featurized_structure_list: copy of the cluster space the structure
        container is based on"""
        return self._cutoff_radius

    
    def GetAllNeigbhors(self, structure):
        '''structure : pymatgen.core.structure.Structure or pymatgen.core.structure.IStructure
        '''
        return structure.get_all_neighbors(self._cutoff_radius)
    
    def Add_Structure(self, structure):
        """structure : pymatgen.core.structure.Structure or pymatgen.core.structure.IStructure
        """
        self._featurized_structure_list.append(self.FitStructure(structure))
        
    def Add_Structure_list(self, structure_list):
        """_structure_list : list(pymatgen.core.structure.Structure or pymatgen.core.structure.IStructure)
        """
        for structure in structure_list:
            self.Add_Structure(structure)
        
        
    def FitStructure(self, structure, **meta_data):
        """Fits a structure.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure or pymatgen.core.structure.IStructure
            the structure to be added
        meta_data : dict
            dict with meta_data about the atoms
        """

        #atoms_copy = atoms.copy()
        neigh_list=self.GetAllNeigbhors(structure)


        '''
        # check if an identical atoms object already exists in the container
        for i, structure in enumerate(self._structure_list):
            if are_configurations_equal(atoms_copy, structure.atoms):
                raise ValueError('Atoms is identical to structure {}'.format(i))
        '''
        
        #logger.debug('Adding structure')
        #M = self._compute_fit_matrix(atoms)
        #structure = FitStructure(atoms_copy, M, **meta_data)
        allneighinfo=[]
        for site,neigh in enumerate(neigh_list):
            if len(neigh)==0:
                print("Site {} has no neighbors within cutoff radius".format(site))
            site_cords =structure[site].coords
            relativeposlist=[]
            rlist=[]
            drlist=[]
            indexlist=[]
            atomtypelist=[]
            for i in neigh:
                relativepos=i.coords-site_cords
                r = np.linalg.norm(relativepos)
                drdx = relativepos[0]/r
                drdy = relativepos[1]/r
                drdz = relativepos[2]/r

                relativeposlist.append(relativepos)
                rlist.append(r)
                drlist.append([drdx,drdy,drdz])
                indexlist.append(i.index)
                atomtypelist.append(i.specie.number)

            allneighinfo.append({'Site':site,'relativepos':np.array(relativeposlist),
                                 'radius':np.array(rlist),'dradius':np.array(drlist),
                                 'neighborSiteNumber':np.array(indexlist),'NeighborAtomType':np.array(atomtypelist)})
        return allneighinfo
        

    def delete_all_structures(self):
        """ Remove all current structures in StructureContainer. """
        self._featurized_structure_list = []



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
                    
            EnergyFeatures.append(np.array(tempENEfeat))
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

            tempdfeat = np.array([-12*site['dradius']*((site['radius'].reshape(-1,1))**(-13)),
                           6*site['dradius']*((site['radius'].reshape(-1,1))**(-7))])

            #tempdfeat.sum(1).T for right feature style
            dfeature.append(tempdfeat.T)
            #.sum(1) for six stress features, xx,yy,zz,xy,yz,xz order
            sfeatures.append(np.vstack([(tempdfeat*-site['relativepos']).T,
                                        (tempdfeat*-np.roll(site['relativepos'], 1,1)).T]))
            
            return feature,dfeature,sfeatures
            

class Cosine_cutoff:
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
        
        self._StructureContainer=StructureContainer
        
        #Need to figure out how to do hyperparameters better
        if hyperparameters == None:
            print("Cosine cutoff is using the Structure container cutoff as cutoff value")
            self._cutoffradius = self._StructureContainer.cutoff_radius
        else:
            self._cutoffradius = hyperparameters['cutoffradius']
            
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
                tempForcefeat=np.zeros([len(self._all_dfeatures[structure_number]),2,3])
                    
            if Return_Stress:
                tempStressfeat=np.zeros([6,2])
                
            for site in range(len(self._all_features[structure_number])):
                if Return_Energy:
                    tempENEfeat+=self._all_features[structure_number][site].sum(0)
                    
                if Return_Forces:
                    tempForcefeat[site] = self._all_dfeatures[structure_number][site].sum(1)
                
                if Return_Stress:
                    tempStressfeat+=self._all_sfeatures[structure_number][site].sum(1)
                    
            EnergyFeatures.append(np.array(tempENEfeat))
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
            feature.append((np.array([np.cos(site['radius']*np.pi/(2*self._cutoffradius))])*(site['radius']<self._cutoffradius)).T)

            tempdfeat = ((site['radius']<self._cutoffradius).reshape(-1,1))*site['dradius']*np.pi*(-np.sin(((site['radius']).reshape(-1,1))*np.pi/(2*self._cutoffradius)))/(2*self._cutoffradius)

            #tempdfeat.sum(1).T for right feature style
            dfeature.append(tempdfeat)
            #.sum(1) for six stress features, xx,yy,zz,xy,yz,zx order (first is force direction, second is displacement direciton)
            sfeatures.append(np.vstack([(tempdfeat*-site['relativepos']).T,
                                        (tempdfeat*-np.roll(site['relativepos'], 1,1)).T])[...,np.newaxis])
                
                #np.vstack([tempdfeat.T[...,np.newaxis],tempdfeat.T[...,np.newaxis]]))
            
            return feature,dfeature,sfeatures
            
