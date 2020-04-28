from rdkit.Chem import Fragments
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors

from scipy.spatial import distance


#####################
# FUNCTIONAL GROUPS #
#####################

def get_mol(smiles):
    return Chem.MolFromSmiles(smiles)


def fp(mol):
    return list(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=512))


def logp(mol):
    return round(Descriptors.MolLogP(mol), 4)

def molwt(mol):
    return Descriptors.ExactMolWt(mol)


def NH2(mol):
    return Fragments.fr_NH2(mol)


def NR2(mol):
    patt1 = Chem.MolFromSmarts('c(N(C)(C))')
    patt2 = Chem.MolFromSmarts('C(N(C)(C))')
    total = 0
    for p in [patt1, patt2]:
        total += len(mol.GetSubstructMatches(p))
    return total


def OH(mol):
    total = Fragments.fr_Ar_OH(mol)
    total += Fragments.fr_Al_OH(mol)
    return total


def OR(mol):
    patt1 = Chem.MolFromSmarts('c(OC)')
    patt2 = Chem.MolFromSmarts('C(OC)')
    total = 0
    for p in [patt1, patt2]:
        total += len(mol.GetSubstructMatches(p))
    return total


def SH(mol):
    return Fragments.fr_SH(mol)


def SR(mol):
    return Fragments.fr_sulfide(mol)


def NO2(mol):
    return Fragments.fr_nitro(mol)


def CN(mol):
    patt = Chem.MolFromSmarts('C#N')
    return len(mol.GetSubstructMatches(patt))


def SO3H(mol):
    patt = Chem.MolFromSmarts('S(=O)(=O)(O)')
    return len(mol.GetSubstructMatches(patt))


def CF3(mol):
    patt = Chem.MolFromSmarts('FC(F)(F)')
    return len(mol.GetSubstructMatches(patt))


def COOH(mol):
    return Fragments.fr_Ar_COO(mol)


def F(mol):
    patt1 = Chem.MolFromSmarts('c(F)')
    patt2 = Chem.MolFromSmarts('C(F)')
    total = 0
    for p in [patt1, patt2]:
        total += len(mol.GetSubstructMatches(p))

    cf3 = Chem.MolFromSmarts('FC(F)(F)')
    total -= 3 * len(mol.GetSubstructMatches(cf3))
    return total

def Br(mol):
    patt = Chem.MolFromSmarts('c(Br)')
    total = len(mol.GetSubstructMatches(patt))
    return total


####################
# REACTION CLASSES #
####################

def suzuki(mol):
    grps = ['[#6;H0;D3:1]B([OH])[OH]', '[#6;H0;D3:2][Br,I,Cl]']

    total = 0
    for grp in grps:
        patt = Chem.MolFromSmarts(grp)
        total += len(mol.GetSubstructMatches(patt))
    return total


def buchwald_hartwig(mol):
    grps = [
        '[Cl,Br,I][c;$(c1:[c,n]:[c,n]:[c,n]:[c,n]:[c,n]:1):1]',
        '[N;$(NC)&!$(N=*)&!$([N-])&!$(N#*)&!$([ND3])&!$([ND4])&!$'
        '(N[c,O])&!$(N[C,S]=[S,O,N]),H2&$(Nc1:[c,n]:[c,n]:[c,n]:[c,n]'
        ':[c,n]:1):2]'
    ]

    total = 0
    for grp in grps:
        patt = Chem.MolFromSmarts(grp)
        total += len(mol.GetSubstructMatches(patt))
    return total


def schotten_baumann_amide(mol):
    grps = [
        '[C;$(C=O):1][OH1]',
        '[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]'
    ]

    total = 0
    for grp in grps:
        patt = Chem.MolFromSmarts(grp)
        total += len(mol.GetSubstructMatches(patt))
    return total


def reductive_amination(mol):
    grps = [
        '[#6:4]-[C;H1,$([CH0](-[#6])[#6]):1]=[OD1]',
        '[N;H2,$([NH1;D2](C)C);!$(N-[#6]=[*]):3]-[C:5]'
    ]

    total = 0
    for grp in grps:
        patt = Chem.MolFromSmarts(grp)
        total += len(mol.GetSubstructMatches(patt))
    return total


def mida_deprotection(mol):
    grps = ['[#6:1]B12OC(=O)C[N+](C)1CC(=O)O2']

    total = 0
    for grp in grps:
        patt = Chem.MolFromSmarts(grp)
        total += len(mol.GetSubstructMatches(patt))
    return total


####################
# SIMILARITY CALCS #
####################

def calculate_pairwise_similarities(df):

    def calc_tanimoto_similarity(fp1, fp2):
        return DataStructs.TanimotoSimilarity(fp1, fp2)

    n_molecules = len(df.smiles)

    similarities = []
    for i in range(n_molecules):
        for j in range(n_molecules - i):
            similarity = 1 - distance.jaccard(df.fp.iloc[i], df.fp.iloc[j])
            similarities.append(similarity)

    return similarities


#################
# PREPROCESSING #
#################

def preprocess(df_from_upload):

    df = df_from_upload
    df.columns = ['smiles']

    preprocess_functions_fgroups = {
        'fp': fp,
        'logp': logp,
        'molwt': molwt,
        'NH2': NH2,
        'NR2': NR2,
        'OH': OH,
        'OR': OR,
        'SH': SH,
        'SR': SR,
        'NO2': NO2,
        'CN': CN,
        'SO3H': SO3H,
        'CF3': CF3,
        'COOH': COOH,
        'F': F,
        'Br': Br,
    }

    preprocess_functions_rxns = {
        'suzuki-miyaura': suzuki,
        'buchwald-hartwig': buchwald_hartwig,
        'schotten-baumann-amide': schotten_baumann_amide,
        'reductive-amination': reductive_amination,
        'mida-deprotection': mida_deprotection,
    }

    df['mol'] = df.smiles.apply(get_mol)

    for name, function in preprocess_functions_fgroups.items():
        df[name] = df.mol.apply(function)

    for name, function in preprocess_functions_rxns.items():
        df[name] = df.mol.apply(function)

    return df
