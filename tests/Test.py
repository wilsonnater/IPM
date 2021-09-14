#struc=mg.Structure([[10,0,0],[0,10,0],[0,0,10]],['Cu','Cu'],[[1,1,1],[2.0,1,1]],coords_are_cartesian=True)
#struc=mg.Structure([[20,0,0],[0,20,0],[0,0,20]],['Cu','Cu','Cu'],[[2,2,2],[0.92,2,2],[3.08,2,2]],coords_are_cartesian=True)
#struc=mg.Structure([[10,0,0],[0,10,0],[0,0,10]],['Cu','Cu','Cu'],[[1,1,1],[1.5,1,1],[7.5,5,1]],coords_are_cartesian=True)

#struc=mg.Structure([[10,0,0],[0,10,0],[0,0,10]],['Cu','Cu','Cu','Cu'],[[1,1,1],[2.0,1,1],[1,6,6],[2.0,6,6]],coords_are_cartesian=True)

#struc=mg.Structure([[10,0,0],[0,10,0],[0,0,10]],['Cu','Cu','Cu'],[[1,1,1],[0.566987298,1.5,1],[1.4330127,1.5,1]],coords_are_cartesian=True)

struc=mg.Structure([[1.9,0,0.1],[0,2,0],[0,0,1.9]],['Cu'],[[1,1,1]],coords_are_cartesian=True)
struc2=mg.Structure([[2.1,0,0],[0,2.1,0],[0,0,2.1]],['Cu'],[[1,1,1]],coords_are_cartesian=True)
struc3=mg.Structure([[4,0,0],[0,2,0],[0,0,2]],['Cu','Cu'],[[1,1,1],[3,1,1]],coords_are_cartesian=True)

#istruc=mg.IStructure([[10,0,0],[0,10,0],[0,0,10]],['Cu','Cu'],[[1,1,1],[1.5,1,1]],coords_are_cartesian=True)



sc=Structure2BodyContainer(2.5,[struc,struc2])
sc.Add_Structure(struc3)
sc.delete_all_structures()
sc.Add_Structure(struc2)


LJ=LJpotetnial(StructureContainer=sc)
LJ.StructureProcesseror()

cut=Cosine_cutoff(StructureContainer=sc)
cut.StructureProcesseror()
