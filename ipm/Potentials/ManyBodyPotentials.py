import numpy as np

class GeneralizedSuttonChenEAM:
    def __init__(self,StructureContainer=[],A=None,B=None,C=None):
        """
        Attributes
        -----------
        hyperparameters : dictionary(hyperparamters)
            dicitonary with relevant hyperparameters for model
        StructureContainer : StructureContainer2Body or StructureContainer3Body(Not yet implemented)
            Structure Container with structures to be featurized
        """
        
        
        Anynone=False
        
        if A ==None:
            self._A = 9
            print("Using A=9 from Sutton Chen Cite")
            Anynone=True
        else:
            self._A = A
        
        if B ==None:
            self._B = 6
            print("Using B=6 from Sutton Chen Cite")
            Anynone=True
        else:
            self._B = B
            
        if C ==None:
            self._C = 0.5
            print("Using C=0.5 from Sutton Chen Cite")
            Anynone=True
        else:
            self._C = C
            
        if Anynone:
            print("CITATION")
            print("See NEED TO PUT IN SOURCE for equation details")
        
            
        self._StructureContainer=StructureContainer
            
        self._all_features=[]
        self._all_dfeatures=[]
        self._all_sfeatures=[]
        
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

        sitedenistyterm=[]

        for site in Structure_NeighborInfo:
            #-6 is sutton chen value for the potential from POET
            sitedenistyterm.append((site['radius']**(-self._B)).sum())

        for site in Structure_NeighborInfo:
            #-9 and 0.5 is sutton chen value for the potential from POET
            feature.append(np.array([(site['radius']**(-self._A)).sum(),-sitedenistyterm[site['Site']]**(self._C)]).T)


            tempsitedensity=[]
            for nieghs in site['neighborSiteNumber']:
                tempsitedensity.append(sitedenistyterm[nieghs])

            #-0.5 is sutton chen value after derivative for the potential from poet
            tempsitedensity=((np.array(tempsitedensity)**(-self._C))+((sitedenistyterm[site['Site']])**(-self._C))).reshape(-1,1)
            print(tempsitedensity)
            #-9 and -6 is sutton chen value for the potential from POET
            tempdfeat = np.array([((-self._A)*site['dradius']*((site['radius'].reshape(-1,1))**((-self._A)-1))),
                           ((self._C)*tempsitedensity*(self._B)*site['dradius']*((site['radius'].reshape(-1,1))**((-self._B)-1)))])

            #tempdfeat.sum(1) for right feature style
            dfeature.append(tempdfeat.T)
            #.sum(1) for six stress features, xx,yy,zz,xy,yz,xz order
            sfeatures.append(np.vstack([(tempdfeat*-site['relativepos']).T,
                                        (tempdfeat*-np.roll(site['relativepos'], 1,1)).T]))

        return feature,dfeature,sfeatures

    
class SuttonChenEAM(GeneralizedSuttonChenEAM):
    def __init__(self,StructureContainer=[]):
        super().__init__(StructureContainer,9,6,0.5)
        
        
