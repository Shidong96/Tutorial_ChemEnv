from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy, MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments

from pymatgen.analysis.chemenv.connectivity.structure_connectivity import StructureConnectivity
from pymatgen.analysis.chemenv.connectivity.connectivity_finder import ConnectivityFinder
#from pymatgen.analysis.chemenv.connectivity.connected_components 

from pymatgen.core.structure import Structure
import logging

logging.basicConfig(filename='test.log',
                    format='%(levelname)s:%(module)s:%(funcName)s:%(message)s',
                    level=logging.DEBUG)

#add structure
s000=Structure.from_file("CONTCAR_000")
S600=Structure.from_file("CONTCAR_600")

# Setup the local geometry finder
lgf = LocalGeometryFinder()
lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True)
lgf.setup_structure(structure=s000)

# Use Voronoi analysis to get the StructureEnvironments
se = lgf.compute_structure_environments(
                                        maximum_distance_factor=1.41,
                                        minimum_angle_factor=0.29,
                                        excluded_atoms=['O', 'K']) #se is a StructureEnvironments object, explore other parameters!

# Analyze the information from Voronoi analysis 
strategy = SimplestChemenvStrategy(distance_cutoff=1.4, angle_cutoff=0.3) #unweighted method
lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)
print(lse.coordination_environments[:])

strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters() #using some kind of defalut setting here, should explore more