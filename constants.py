from collections import OrderedDict


def fromReactionFullNameToReactantName(reactionName, reactantNum):
    return ("react%d_" + reactionName.replace(" ", "_"))% reactantNum

def fromReactantNameToReactionFullName(reactantName):
    return reactantName.replace("react1", "").replace("react2", "")

REACTANT1_PREFIX = fromReactionFullNameToReactantName("", 1)
REACTANT2_PREFIX = fromReactionFullNameToReactantName("", 2)

reaction_smarts = OrderedDict({
    "Amidation": "[#6:1](=[#8:2])-[#8].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]",
    "Amide_schotten-baumann": "[#7;H2,H1:3].[#6:1](=[#8:2])-[#17]>>[#6:1](=[#8:2])-[#7:3]",
    "Reductive_amination": "[#6:2](=[#8])(-[#6:1]).[#7;H3,H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]",
    "N-nucleophilic_aromatic_substitution": "[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]",
    "Sp2-sp2_Suzuki_coupling": "[c:1]-[F,Cl,Br,I].[#6:2]-[B]>>[c:1]-[#6:2]"
})

REACTIONS_NAMES = list(reaction_smarts.keys())
