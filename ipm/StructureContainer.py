
import numpy as np
#import pymatgen as mg

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
                                 'neighborSiteNumber':np.array(indexlist),'NeighborAtomType':np.array(atomtypelist),
                                 'volume':structure.volume})
        return allneighinfo
        

    def delete_all_structures(self):
        """ Remove all current structures in StructureContainer. """
        self._featurized_structure_list = []


'''        
Class 2BodyContainer:
    def __init__(self,cutoff_radius, structure):
        
        
        self._cutoff_radius = cutoff_radius
        self._3body=False
        
        self.allneighinfo=[]
    
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
    
    
    
        # check if an identical atoms object already exists in the container
        #for i, structure in enumerate(self._structure_list):
        #    if are_configurations_equal(atoms_copy, structure.atoms):
        #        raise ValueError('Atoms is identical to structure {}'.format(i))

        
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
    
            #allneighinfo.append({'Site':site,'relativepos':np.array(relativeposlist),
            #                     'radius':np.array(rlist),'dradius':np.array(drlist),
            #                     'neighborSiteNumber':np.array(indexlist),'NeighborAtomType':np.array(atomtypelist)})
            self.allneighinfo.append(site,np.array(relativeposlist),np.array(rlist),np.array(drlist),np.array(indexlist),np.array(atomtypelist))
            '''
            
            
            
            