from collections import OrderedDict
import json

# PROVIDE PATH TO REACTION SMARTS CONSTANTS
reaction_constants_path = "/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/D68EV3CPROA/RXN_SMARTS_CONSTANTS.json"
with open(reaction_constants_path) as f:
    reaction_smarts = json.load(f)

def fromReactionFullNameToReactantName(reactionName, reactantNum):
    return ("react%d_" + reactionName.replace(" ", "_"))% reactantNum

def fromReactantNameToReactionFullName(reactantName):
    return reactantName.replace("react1", "").replace("react2", "")

REACTANT1_PREFIX = fromReactionFullNameToReactantName("", 1)
REACTANT2_PREFIX = fromReactionFullNameToReactantName("", 2)

REACTIONS_NAMES = list(reaction_smarts.keys())

############### REACTION SMARTS GRAVEYARD (not used) ################
# # IMPORTANT: DO NOT CHANGE THE ORDER, OTHERWISE YOU WILL MESS UP REACTIONSMARTS.PY SCRIPTS
# # THESE ARE VERY SPECIFIC REACTION SMARTS
# reaction_smarts = OrderedDict({
#     # DONE Check amine is not an amide, check that it is not a sulfonamide
#     "Amidation": "[#6:1](=[#8:2])-[#8;H1].[#7&X3;H2,H1;!$(NC=*);!$(NS):3]>>[#6:1](=[#8:2])-[#7:3]",
#     # Need to specify
#     "Amide_schotten-baumann": "[#7;H2,H1:3].[#6:1](=[#8:2])-[#17]>>[#6:1](=[#8:2])-[#7:3]",
#     # Need to specify
#     "Reductive_amination": "[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]",
#     # Need to specify
#     "N-nucleophilic_aromatic_substitution": "[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]",
#     # Need to specify (NEED TO TEST)
#     "Sp2-sp2_Suzuki_coupling": "[c:1]-[F,Cl,Br,I].[#6:2]-[B]>>[c:1]-[#6:2]",
#     # DONE Formation of urea from two primary or secondary amines, specifically not amides
#     "Formation_of_urea_from_two_amines": "[N&X3;H2,H1;!$(NC=*):3].[N&X3;H2,H1;!$(NC=*):4]>>[#7;H1:3]-[#6](=[#8])-[#7;H1:4]",
#     # DONE, high specficity for sulfone and primary or secondary amine
#     "Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)": "[#16X4:1](=[OX1:2])(=[OX1:3])([#6:4])[#17,#9].[#7&X3;H2,H1;!$(NC=*);!$(NS):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6:4])-[#7:5]",
#     # DONE, only4 looking for OTf group (not halide) and primary or secondary amine
#     #"Buchwald-Hartwig_amination": "[#16X4:1](=[OX1:2])(=[OX1:3])([OX2H0:4][ar:5])[#6X4:6](-F)(-F)(-F).[#7&X3;H2,H1;!$(NC=*);!$(NS);!$(NN):7]>>[ar:5]-[N:7]"
#     # Using halide configuration
#     "Buchwald-Hartwig_amination": "[ar:1]-[Cl,Br,I].[#7&X3;H2,H1;!$(NC=*);!$(NS);!$(NN):7]>>[ar:1]-[N:7]"
# })

# REACTION_ATOM_NUMS = OrderedDict({
#     "Amidation": [3, 1],
#     "Amide_schotten-baumann": [1, 3],
#     "Reductive_amination": [3, 1],
#     "N-nucleophilic_aromatic_substitution": None,
#     "Sp2-sp2_Suzuki_coupling": [2, 2],
#     "Formation_of_urea_from_two_amines": [1],
#     "Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)": [5, 1]
# })