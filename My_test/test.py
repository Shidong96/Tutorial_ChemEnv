from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy, MultiWeightsChemenvStrategy, NormalizedAngleDistanceNbSetWeight
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments

from pymatgen.analysis.chemenv.connectivity.structure_connectivity import StructureConnectivity
from pymatgen.analysis.chemenv.connectivity.connectivity_finder import ConnectivityFinder
#from pymatgen.analysis.chemenv.connectivity.connected_components 

from pymatgen.core.structure import Structure
import logging

from sympy import E

logging.basicConfig(filename='test.log',
                    format='%(levelname)s:%(module)s:%(funcName)s:%(message)s',
                    level=logging.DEBUG)

#add structure
s000=Structure.from_file("CONTCAR_000")
s600=Structure.from_file("CONTCAR_600")

# Setup the local geometry finder
lgf = LocalGeometryFinder()
lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True)
lgf.setup_structure(structure=s600)

# Use Voronoi analysis to get the StructureEnvironments
se = lgf.compute_structure_environments(
                                        maximum_distance_factor=1.41,
                                        minimum_angle_factor=0.29,
                                        excluded_atoms=['O', 'K']) #se is a StructureEnvironments object, explore other parameters

# Choose strategy to analyze the information from Voronoi analysis 
#strategy = SimplestChemenvStrategy(distance_cutoff=1.4, angle_cutoff=0.3) #unweighted method
strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters() #using some kind of defalut setting here, should explore more

# Now it's applying the strategy on the se object to get a lse object? the name is confusing, get a se from a se.
lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)
#print(lse.strategy, lse.coordination_environments[:]) # this output only has strategy and ce, input is strategy and se

sc = ConnectivityFinder().get_structure_connectivity(lse)
sc.setup_environment_subgraph(["O","O","T","T","T","T","T","T"])
print (sc.environment_subgraph)
sc.print_links()
connectcomp = sc.get_connected_components
print (connectcomp)

############
# The connectivity package is not fully written, so far I can't see how this can work, especially on my problem. shit.
############

