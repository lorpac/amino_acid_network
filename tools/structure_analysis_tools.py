import sys
import os
from Bio.PDB import PDBList
import networkx as nx
import biographs as bg
import csv
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import tools.network_visualization
import pickle
import networkx.algorithms.clique as clique
import tools.sidechain_length as scl
import tools.amino_acids_conversion as aaconv
import glob
from biopandas.pdb import PandasPdb
from shutil import copyfile


def download_pdb(pdb_id, pdbs_path):
    """Downloads a pdb file
    
    Parameters
    ----------
    pdb_id : str
        pdb id
    pdbs_path: str, optional
        path of folder containing the pdbs (default is "pdbs")
    """
    pdbl = PDBList(obsolete_pdb=True)
    pdbl.download_pdb_files(pdb_id, file_format="pdb", pdir=pdbs_path)


def get_pdb_path(pdb_id, pdbs_path="pdbs"):
    """Gets the location of a pdb file
    
    Parameters
    ----------
    pdb_id : str
        pdb id
    pdbs_path: str, optional
        path of folder containing the pdbs (default is "pdbs")
    
    Returns
    -------
    str
        pdb file path
    str
        True if the pdb has been downloaded
    """

    downloaded = False

    if not pdbs_path:
        pdbs_path = "pdbs"

    package_path = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]
    abs_path = os.path.join(package_path, "data", pdbs_path)
    abs_file_path = os.path.join(abs_path, pdb_id + ".*")

    if len(glob.glob(abs_file_path)) == 0:
        abs_file_path = os.path.join(abs_path, "pdb" + pdb_id + ".*")
        downloaded = True

        if len(glob.glob(abs_file_path)) == 0:
            os.makedirs(abs_path, exist_ok=True)
            download_pdb([pdb_id], abs_path)

        else:
            pdb_id = "pdb" + pdb_id
            abs_file_path = os.path.join(abs_path, pdb_id + ".*")
    pdb_path = glob.glob(abs_file_path)[0]

    return pdb_path, downloaded


def remove_hydrogens(pdb_file):
    """Remove hydrogens from a pdb file (saves copy of original file)
    
    Parameters
    ----------
    pdb_file : str
        pdb file path
    """

    pdb_id = pdb_file.rsplit(".", 1)[0]
    copyfile(pdb_file, pdb_file.rsplit(".", 1)[0] + "_ORIG.pdb")

    with open(pdb_file, "r") as f:
        pdb = [line.rsplit("\n")[0] for line in f]
    pdb_new = []
    for line in pdb:
        if line[0:4] == "ATOM" and (line[13] == "H" or line[12] == "H"):
            pass
        else:
            pdb_new.append(line)
    with open(pdb_file, "w") as f:
        f.writelines([line + "\n" for line in pdb_new])


def list_relations(dim="all"):
    """Returns the list of relations to analyze
    
    Parameters
    ----------
    dim : str, optional
        "1D", "2D", "1-2D", "3-4D", "all" or "" (default is "all")
        "" is interpreted as "all"
    
    Returns
    -------
    list
        list of relations to analyze
    """

    if dim == "2D":
        rel_list = ["2D2", "2D3", "2D4"]
    elif dim == "1D":
        rel_list = ["1D"]
    elif dim == "1-2D":
        rel_list = ["1D", "2D2", "2D3", "2D4"]
    elif dim == "3-4D":
        rel_list = ["3D", "4D"]
    elif dim == "intramolecular":
        rel_list = ["1D", "2D2", "2D3", "2D4", "3D"]
    elif dim == "all" or dim == "" :
        rel_list = ["1D", "2D2", "2D3", "2D4", "3D", "4D"]
    else:
        rel_list = [dim]

    return rel_list


def assign_secondary_structure(pdb):
    """Returns the secondary structure elements of a pdb
    
    Parameters
    ----------
    pdb : str
        pdb file path
    
    Returns
    -------
    dict
        dictionary of secondary structure elements
    """

    ppdb = PandasPdb().read_pdb(pdb)

    secondary_structure = {}

    helices_from_pdb = ppdb.df["OTHERS"][ppdb.df["OTHERS"]["record_name"] == "HELIX"][
        "entry"
    ]
    for helix in helices_from_pdb:
        identifier_h = helix[5:8].strip()
        initial_chain_h = helix[13].strip()
        initial_pos_h = helix[16:19].strip()
        final_pos_h = helix[28:31].strip()
        for i in range(int(initial_pos_h), int(final_pos_h) + 1):
            secondary_structure[initial_chain_h + str(i)] = (
                "helix" + identifier_h + "-" + initial_chain_h
            )

    sheets_from_pdb = ppdb.df["OTHERS"][ppdb.df["OTHERS"]["record_name"] == "SHEET"][
        "entry"
    ]
    for sheet in sheets_from_pdb:
        identifier_s = sheet[6:8].strip()
        initial_chain_s = sheet[15].strip()
        initial_pos_s = sheet[17:20].strip()
        final_pos_s = sheet[28:31].strip()
        for i in range(int(initial_pos_s), int(final_pos_s) + 1):
            secondary_structure[initial_chain_s + str(i)] = (
                "sheet" + identifier_s + "-" + initial_chain_s
            )

    mol = bg.Pmolecule(pdb)
    net = mol.network()

    residues_type = {}
    for residue in mol.model.get_residues():
        res_type = residue.resname
        res_pos = residue.parent.id + str(residue.id[1])
        residues_type[res_pos] = res_type

    residues = list(net.nodes)  # assume they are ordered
    last_structure = None
    last_chain = None
    i = 0
    for residue in residues:
        chain = residue[0]
        try:
            structure = secondary_structure[residue]
            if structure != last_structure:
                i += 1
        except KeyError:
            if chain != last_chain:
                i += 1
            structure = "loop" + str(i)
            secondary_structure[residue] = structure
        last_structure = structure
        last_chain = chain

    return secondary_structure


def get_neighbor_structure_relation(secondary_structure, u, v):
    """Returns the relation (1D, 2D, 3D, 4D) between to neighboring nodes
    
    Parameters
    ----------
    secondary_structure : dict
        dictionary of secondary structure elements
    u: str
        node label
    v: str
        node label
    
    Returns
    -------
    str
        relation (1D, 2D, 3D, 4D)
    """

    chain_u = u[0]
    chain_v = v[0]
    pos_u = u[1::]
    pos_v = v[1::]
    struct_u = secondary_structure[u]
    struct_v = secondary_structure[v]

    if chain_u == chain_v:
        dist = np.abs(int(pos_u) - int(pos_v))
        if dist == 1:
            relation = "1D"
        elif struct_u == struct_v:
            if dist < 5:
                relation = "2D" + str(dist)
            else:
                relation = "3D"
        else:
            relation = "3D"
    else:
        relation = "4D"

    return relation


def import_mutations(folder_path):
    """Imports a list (or lists) of single mutations
    
    Parameters
    ----------
    folder_path : str
        The path of the colder containing the mutation list files: 
        "single_mutations.csv" and eventually "single_mutations_non_pathogenic.csv".
    
    Returns
    -------
    dict
        dictionary of type of mutation (pathogenic or non pathogenic)
    list
        list of positions having a pathogenic mutation
    list
        list of positions having a non pathogenic mutation
    list
        list of positions having both a pathogenic and a non pathogenic mutation
    list
        list of pathogenic mutations
    list
        list of non pathogenic mutations
    """

    path = os.path.split(folder_path)[0]

    mutations_type = {}
    mutations = []
    pathogenic = []
    pathogenic_names = []
    non_pathogenic = []
    non_pathogenic_names = []
    both = []
    if "single_mutations.csv" in os.listdir(
        path
    ) and "single_mutations_non_pathogenic.csv" in os.listdir(path):

        with open(os.path.join(path, "single_mutations.csv"), "r") as f:
            reader = csv.reader(f)
            for line in reader:
                pathogenic.append(line[1])
                name = "".join(line).replace(" ", "")
                pathogenic_names.append(name)
                mutations_type[name] = "pathogenic"
                mutations.append(line)

        with open(os.path.join(path, "single_mutations_non_pathogenic.csv"), "r") as f:
            reader = csv.reader(f)
            for line in reader:
                non_pathogenic.append(line[1])
                name = "".join(line).replace(" ", "")
                non_pathogenic_names.append(name)
                mutations_type[name] = "non_pathogenic"
                mutations.append(line)

        def remove_duplicates(values):
            output = []
            seen = set()
            for value in values:
                # If value has not been encountered yet,
                # ... add it to both list and set.
                if value not in seen:
                    output.append(value)
                    seen.add(value)
            return output

        pathogenic = remove_duplicates(pathogenic)
        non_pathogenic = remove_duplicates(non_pathogenic)

        for p in pathogenic:
            if p in non_pathogenic:
                pathogenic.pop(pathogenic.index(p))
                non_pathogenic.pop(non_pathogenic.index(p))
                both.append(p)

        print(len(pathogenic_names), " pathogenic mutations loaded")
        print(len(non_pathogenic_names), " non pathogenic mutations loaded")
        print(len(pathogenic), " positions involved in pathogenic mutations")
        print(len(non_pathogenic), " positions involved in non-pathogenic mutations")
        print(
            len(both),
            " positions involved in both pathogenic and non-pathogenic mutations",
        )

    else:
        print("No mutations loaded")

    return (
        mutations_type,
        pathogenic,
        non_pathogenic,
        both,
        pathogenic_names,
        non_pathogenic_names,
    )


def list_aa(sort_by_size=False, sort_by_sidechain=False):
    """Imports the list of amino acids.
        Requires "aminoacids_size.txt" to be present in the directory.
    
    Parameters
    ----------
    sort_by_size : boolean, optional
        if True, sorts the amino acids by number of atoms (default is False)
    sort_by_sidechain : boolean, optional
        if True, sorts the amino acids by side chain length (default is False)
    
    Returns
    -------
    list
        list of amino acids
    """

    package_path = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]

    if sort_by_size:
        aa_path = os.path.join(package_path, "data", "aminoacids_size.txt")
    elif sort_by_sidechain:
        aa_path = os.path.join(package_path, "data", "aminoacids_sidechain.txt")
    else:
        aa_path = os.path.join(package_path, "data", "aminoacids.txt")

    with open(aa_path, "r") as f:
        amino_acids = [aa.rsplit("\n")[0] for aa in f]

    return amino_acids


def dict_aa_schl(amino_acids):
    """Creates a dictionary of amino acids side chain lengths.
    
    Parameters
    ----------
    amino_acids : list
        list of amino acids
    
    Returns
    -------
    dict
        dictionary of amino acids side chain lengths
    """

    amino_acids1l = [aaconv.three2one(aa) for aa in amino_acids]
    schl_dict = scl.dict_classif
    amino_acids_schl = [schl_dict[aa] for aa in amino_acids1l]

    return amino_acids_schl

def get_chains(mol):
    """Gets the chain of a Pmolecule object.
        
    Parameters
    ----------
    mol : Pmolecule
        Biographs Pmolecule object
    
    Returns
    -------
    list
        list of chains id's
    """

    chains = mol.model.child_list
    chains_id = [c.get_id() for c in chains]

    return chains_id

def create_aa_network(
    pdb_id,
    rel_list,
    folder_path,
    selected_positions=None,
    cutoff=5,
    kw_reprod=False,
    k_w=None,
    db_1_name=None,
    db_2_name=None,
    separate_jaccard=False,
    separate_weights=False,
    pdbs_path="pdbs",
    save_csv=True,
    remove_hydrogen_atoms=False
):
    """Creates the amino acid network from a pdb id.
        
    Parameters
    ----------
    pdb_id : str
        pdb id of the protein
    rel_list: list
        list of relation (1D, 2D, 3D, 4D) to consider.
    folder_path: str
        path of the output folder
    selected_positions: None or list, optional
        list of sequence positions to consider. If None, all positions are considered (default is None)
    cutoff: int, optional
        cutoff threshold for the connection of nodes in the amino acid network (dafault is 5).
    kw_reprod: boolean, optional
        if True, adds a column that checks if (k, w) of a node exist in a database of nodes of robust proteins (default is False).
    k_w: dict or None, optional
        dictionary of weight values associated to each k value in the database of nodes of robust proteins (default is None).
    db_1_name: str, optional
        name of the database of type 1. If None and save_csv is True, saves the database as database_pos_1.csv (default is None).
    db_2_name: str, optional
        name of the database of type 2. If None and save_csv is True, saves the database as database_pos_1.csv (default is None).
    separate_jaccard: boolean, optional
        if True, separates the Jaccard vector in the database of type 2 based on the size of he neighbors (small, medium, large) (default is False).
    separate_weights: boolean, optional
        if True, separates the weight in the database of type 2 based on the size of he neighbors (small, medium, large) (default is False).
    pdbs_path: str, optional
        path of the pdb files folder (default is "pdbs")
    save_csv: boolean, optional
        if True, saves the database as a csv file in the directory specified by folder_path (default is True).
    remove_hydrogen_atoms: boolean, optional
        if True, saves removes the hydrogen atoms from the pdb file (default is True).
    Returns
    -------
    Graph
        NetworkX Graph object
    DataFrame
        pandas DataFrame object
    DataFrame
        pandas DataFrame object
    Pmolecule
        Biograph Pmolecule object
    boolean
        True if the pdb file was downloaded
    """

    amino_acids = list_aa()

    if separate_jaccard:
        amino_acids_schl = dict_aa_schl(amino_acids)

    if not separate_jaccard:  # DO I WANT TO LEAVE IT LIKE THIS?
        separate_weights = False

    pdb, downloaded = get_pdb_path(pdb_id, pdbs_path)

    if remove_hydrogen_atoms:
        remove_hydrogens(pdb)

    # initialize databases to report
    database_1 = []
    database_2 = []

    mol = bg.Pmolecule(pdb)
    net = mol.network(cutoff=cutoff)

    # take only selected positions:
    if selected_positions:
        for node in list(net.nodes):
            pos = int(node[1::])
            if pos not in selected_positions:
                net.remove_node(node)
    else:
        positions = [int(node[1::]) for node in list(net.nodes)]
        pos_min = min(positions)
        pos_max = max(positions)
        selected_positions = range(pos_min, pos_max + 1)

    secondary_structure = assign_secondary_structure(pdb)

    residues_dict = {}
    for residue in mol.model.get_residues():
        res_type = residue.resname.strip()
        if len(res_type) < 3:
            res_type = aaconv.one2three(res_type)
        res_pos = residue.parent.id + str(residue.id[1])
        residues_dict[res_pos] = res_type

    for residue in mol.model.get_residues():
        adj_vector = [0] * 20
        weight_vector = [0] * 20
        node_name = residue.parent.id + str(residue.id[1])
        deg = nx.degree(net, residue.parent.id + str(residue.id[1]))
        if deg == 0:
            net.remove_node(residue.parent.id + str(residue.id[1]))
        else:
            weight = nx.degree(
                net, residue.parent.id + str(residue.id[1]), weight="weight"
            )
            restype = residue.resname
            resname = (
                os.path.split(pdb)[1][:-4] + residue.parent.id + str(residue.id[1])
            )
            size = len(residue)
            seqpos = residue.id[1]
            if seqpos not in selected_positions:
                continue
            structure = secondary_structure[node_name]

            # check how many other aas can have the same k and w in the database
            if kw_reprod:
                n_others = 0
                for aa in amino_acids:
                    if aa != restype:
                        try:
                            [w_min, w_max] = k_w[aa][deg]
                            if w_min <= weight and w_max >= weight:
                                n_others += 1
                        except KeyError:
                            pass

            if separate_weights:
                w_separated = {"s": 0, "m": 0, "l": 0}

            line_neigh = []
            for neighbor in list(nx.neighbors(net, node_name)):
                neighbor_type = residues_dict[neighbor]
                edge_weight = nx.edges(net)[(node_name, neighbor)]["weight"]
                aa_num = amino_acids.index(neighbor_type)
                relation = get_neighbor_structure_relation(
                    secondary_structure, node_name, neighbor
                )
                # select only edges of desired relations
                if relation in rel_list:
                    adj_vector[aa_num] += 1
                    weight_vector[aa_num] += edge_weight
                    line_neigh.append(neighbor)
                    line_neigh.append(neighbor_type)
                    line_neigh.append(edge_weight)
                    line_neigh.append(relation)
                else:
                    net.remove_edge(neighbor, node_name)

                # separate weights
                if separate_weights:
                    neigh_size = scl.dict_classif[aaconv.three2one(neighbor_type)]
                    w_separated[neigh_size] += edge_weight

            # check if the residue became of degree zero:
            deg = nx.degree(net, residue.parent.id + str(residue.id[1]))

            if deg == 0:
                net.remove_node(residue.parent.id + str(residue.id[1]))
            else:
                weight = nx.degree(
                    net, residue.parent.id + str(residue.id[1]), weight="weight"
                )

                line = [
                    resname,
                    node_name,
                    seqpos,
                    restype,
                    deg,
                    weight,
                    weight / deg,
                    size,
                    structure,
                ]
                line_2 = [
                    resname,
                    node_name,
                    seqpos,
                    restype,
                    deg,
                    weight,
                    weight / deg,
                    size,
                    structure,
                ]

                if kw_reprod:
                    line.append(n_others)
                    line_2.append(n_others)

                line += line_neigh

                database_1.append(line)
                line_2 += adj_vector

                if separate_jaccard:
                    sep_jaccard_dict = {"s": 0, "m": 0, "l": 0}
                    for k in range(len(amino_acids)):
                        num = adj_vector[k]
                        size = amino_acids_schl[k]
                        sep_jaccard_dict[size] += num
                    sep_jacc = [
                        sep_jaccard_dict["s"],
                        sep_jaccard_dict["m"],
                        sep_jaccard_dict["l"],
                    ]
                    line_2 += sep_jacc

                    if separate_weights:
                        w_separated = list(w_separated.values())
                        line_2 += w_separated

                database_2.append(line_2)

    sortedlist_pos = sorted(database_1, key=lambda row: row[2])

    sortedlist_pos_2 = sorted(database_2, key=lambda row: row[2])

    if not db_1_name:
        db_1_name = "database_pos_1.csv"
    if not db_2_name:
        db_2_name = "database_pos_2.csv"

    if kw_reprod:
        columns1 = [
            "Residue name",
            "Position",
            "Sequence position",
            "Type of residue",
            "Degree",
            "Weight",
            "Weight/Degree",
            "Atomic number",
            "Secondary structure",
            "N. others",
            "Neighbor position",
            "Neighbor type",
            "Pairwise weight",
            "Relation",
        ]
        columns2 = [
            "Residue name",
            "Position",
            "Sequence position",
            "Type of residue",
            "Degree",
            "Weight",
            "Weight/Degree",
            "Atomic number",
            "Secondary structure",
            "N. others",
        ] + amino_acids
    
    else:
        columns1 = [
            "Residue name",
            "Position",
            "Sequence position",
            "Type of residue",
            "Degree",
            "Weight",
            "Weight/Degree",
            "Atomic number",
            "Secondary structure",
            "Neighbor position",
            "Neighbor type",
            "Pairwise weight",
            "Relation",
        ]
        columns2 = [
            "Residue name",
            "Position",
            "Sequence position",
            "Type of residue",
            "Degree",
            "Weight",
            "Weight/Degree",
            "Atomic number",
            "Secondary structure"
        ] + amino_acids

    if separate_jaccard:
        columns2 += ["small neighbors", "medium neighbors", "large neighbors"]

        if separate_weights:
            columns2 += ["w_small", "w_medium", "w_large"]

    if save_csv:
        db_1_path = os.path.join(folder_path, db_1_name)
        with open(db_1_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(columns1)
            writer.writerows(sortedlist_pos)

        db_2_path = os.path.join(folder_path, db_2_name)
        with open(db_2_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(columns2)
            writer.writerows(sortedlist_pos_2)

    # database 1 has to have all rows of the same lenght to be read as a dataframe
    lengths = [len(row) for row in database_1]
    max_length = max(lengths)

    db_1 = []
    missing_header = max_length - len(columns1)
    for i in range(int(missing_header / 4)):
        columns1.append("Neighbor position")
        columns1.append("Neighbor type")
        columns1.append("Pairwise weight")
        columns1.append("Relation")

    for row in database_1:
        missing = max_length - len(row)
        for i in range(int(missing)):
            row.append("-")
        db_1.append(row)

    db_1 = pd.DataFrame(db_1, columns=columns1)

    db_2 = pd.DataFrame(database_2, columns=columns2)

    return net, db_1, db_2, mol, downloaded


def create_databases(
    pdbs,
    rel_list,
    folder_path,
    pdbs_path,
    cutoff=5,
    db_1_name=None,
    db_2_name=None,
    separate_jaccard=False,
    separate_weights=False,
    save_csv=True,
):
    """Creates a databases of node properties from a list of pdbs.
    
    Parameters
    ----------
    pdbs : list
        list of pdb id's
    rel_list: list
        list of relation (1D, 2D, 3D, 4D) to consider.
    folder_path: str
        path of the output folder
    pdbs_path: str
        path of the pdb files folder
    cutoff: int, optional
        cutoff threshold for the connection of nodes in the amino acid network (dafault is 5).
    db_1_name: str, optional
        name of the database of type 1. If None and save_csv is True, saves the database as database_pos_1.csv (default is None).
    db_2_name: str, optional
        name of the database of type 2. If None and save_csv is True, saves the database as database_pos_1.csv (default is None).
    separate_jaccard: boolean, optional
        if True, separates the Jaccard vector in the database of type 2 based on the size of he neighbors (small, medium, large) (default is False).
    separate_weights: boolean, optional
        if True, separates the weight in the database of type 2 based on the size of he neighbors (small, medium, large) (default is False).
    save_csv: boolean, optional
        if True, saves the database as a csv file in the directory specified by folder_path (default is True).
    Returns
    -------
    DataFrame
        pandas DataFrame object
    DataFrame
        pandas DataFrame object
    """

    os.makedirs(folder_path, exist_ok=True)

    amino_acids = list_aa()

    if separate_jaccard:
        amino_acids_schl = dict_aa_schl(amino_acids)

    if not separate_jaccard:  # DO I WANT TO LEAVE IT LIKE THIS?
        separate_weights = False

    # initialize databases to report
    databases_1 = []
    databases_2 = []

    for pdb_id in pdbs:

        net, db_1, db_2, mol, downloaded = create_aa_network(
            pdb_id,
            rel_list,
            folder_path,
            cutoff=cutoff,
            separate_jaccard=separate_jaccard,
            separate_weights=separate_weights,
            pdbs_path=pdbs_path,
            save_csv=False,
        )

        databases_1.append(db_1)
        databases_2.append(db_2)

    # databases of type 1 need to have the same column length
    lengths = [len(database.columns) for database in databases_1]
    max_length = max(lengths)

    databases_1_new = []
    for db in databases_1:
        columns = list(db.columns)
        missing_header = max_length - len(columns)
        for i in range(int(missing_header / 4)):
            columns.append("Neighbor position")
            columns.append("Neighbor type")
            columns.append("Pairwise weight")
            columns.append("Relation")

        db_rows = list(db.values)
        db_rows_new = []
        for row in db_rows:
            missing = max_length - len(row)
            for i in range(int(missing)):
                row = np.append(row, '-')
            db_rows_new.append(row)
        db = pd.DataFrame(db_rows_new, columns=columns)
        databases_1_new.append(db)
    databases_1 = databases_1_new

    database_1 = pd.concat(databases_1)
    database_2 = pd.concat(databases_2)
    
    database_1 = database_1.sort_values(by="Type of residue")
    database_2 = database_2.sort_values(by="Type of residue")
    
    db1_path = os.path.join(folder_path, "database1.csv")
    db2_path = os.path.join(folder_path, "database2.csv")

    database_1.to_csv(db1_path, header=True, index=False)
    database_2.to_csv(db2_path, header=True, index=False)

    return database_1, database_2


def create_database_nodes(
    pdbs, folder_path, pdbs_path, cutoff=5, save_csv=True, db_name=None
):
    """Creates a database of node properties from a list of pdbs.
    
    Parameters
    ----------
    pdbs : list
        list of pdb id's
    folder_path: str
        path of the output folder
    pdbs_path: str
        path of the pdb files folder
    cutoff: int, optional
        cutoff threshold for the connection of nodes in the amino acid network (dafault is 5).
    save_csv: boolean, optional
        if True, saves the database as a csv file in the directory specified by folder_path (default is True).
    db_name: str, optional
        name of the database. If None and save_csv is True, saves the database as database_nodes.csv (default is None).
    
    Returns
    -------
    DataFrame
        pandas DataFrame object
    """

    amino_acids = list_aa()
    # initialize database to report
    database = []
    for pdb_id in pdbs:

        pdb, downloaded = get_pdb_path(pdb_id, pdbs_path)

        mol = bg.Pmolecule(pdb)
        net = mol.network(cutoff=cutoff)

        secondary_structure = assign_secondary_structure(pdb)

        residues_dict = {}
        for residue in mol.model.get_residues():
            res_type = residue.resname.strip()
            if len(res_type) < 3:
                res_type = aaconv.one2three(res_type)
            res_pos = residue.parent.id + str(residue.id[1])
            residues_dict[res_pos] = res_type

        for residue in mol.model.get_residues():
            node_name = residue.parent.id + str(residue.id[1])
            deg = nx.degree(net, residue.parent.id + str(residue.id[1]))
            if deg == 0:
                net.remove_node(residue.parent.id + str(residue.id[1]))
            else:
                weight = nx.degree(
                    net, residue.parent.id + str(residue.id[1]), weight="weight"
                )
                restype = residue.resname
                resname = (
                    os.path.split(pdb)[1][:-4] + residue.parent.id + str(residue.id[1])
                )
                size = len(residue)
                structure = secondary_structure[node_name]
                if structure[0] == "h":
                    structure = "helix"
                elif structure[0] == "s":
                    structure = "sheet"
                else:
                    structure = "loop"

                line = [resname, restype, deg, weight, weight / deg, size, structure]

                # divide 1D from other neighbors
                w1D = 0
                k1D = 0  # can be 1 or 2
                wOTH = 0
                kOTH = 0

                for neighbor in list(nx.neighbors(net, node_name)):
                    edge_weight = nx.edges(net)[(node_name, neighbor)]["weight"]
                    relation = get_neighbor_structure_relation(
                        secondary_structure, node_name, neighbor
                    )

                    if relation == "1D":
                        w1D += edge_weight
                        k1D += 1
                    else:
                        wOTH += edge_weight
                        kOTH += 1

                if k1D == 0 or kOTH == 0:
                    # print(resname)
                    continue

                nw1D = w1D / k1D
                nwOTH = wOTH / kOTH

                not_terminal = k1D == 2

                line += [nw1D, nwOTH, not_terminal]
                database.append(line)

    if not db_name:
        db_name = "database_nodes.csv"
    elif db_name[-3::] != ".csv":
        db_name += ".csv"

    columns = [
        "Residue name",
        "Type of residue",
        "Degree",
        "Weight",
        "Weight/Degree",
        "Atomic number",
        "Secondary structure",
        "NW1D",
        "NWothers",
        "Not terminal",
    ]

    if save_csv:
        db_path = os.path.join(folder_path, db_name)
        with open(db_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(columns)
            writer.writerows(database)

    db = pd.DataFrame(database, columns=columns)

    return db


def divided_nw_lists(db_2):
    """Returns the 3 lists of observed nw values of amino acids belonging to
    helices, sheets and loops, respectively.
        
    Parameters
    ----------
    db_2 : DataFrame
        Pandas DataFrame containing node properties
    Returns
    -------
    list
        list of Nw in alpha helices
    list
        list of Nw in beta sheets
    list
        list of Nw in loops
    """

    nw_helix = []
    nw_sheet = []
    nw_loop = []

    for _, row in db_2.iterrows():
        nw = row["Weight/Degree"]
        secstr = row["Secondary structure"]
        if secstr[0] == "h":
            nw_helix.append(nw)
        elif secstr[0] == "s":
            nw_sheet.append(nw)
        else:
            nw_loop.append(nw)

    return nw_helix, nw_sheet, nw_loop


def herfindhal_index(x):
    """Returns Herfindhal index of a list.
        
    Parameters
    ----------
    x : list
        list of values
    Returns
    -------
    float
        Herfindhal index
    """
    num = sum([v ** 2 for v in x])
    den = (sum(x)) ** 2
    H = num / den
    return H


def plot_nw_hist(
    nw_helix,
    nw_sheet,
    nw_loop,
    folder_path,
    bins=None,
    inset=None,
    file_name="neighborhood_watch_distribution.pdf",
):
    """Plots the histogram of the Nw values based on the secondary structure of the node.
        
    Parameters
    ----------
    nw_helix : list
        list of Nw in alpha helices
    nw_sheet : list
        list of Nw in beta sheets
    nw_loop : list
        list of Nw in loops
    folder_path: str
        output folder path
    bins: int or list or None, optional
        number of list of bins. In None, bins of size 1 from zero to the max value of Nw are employed (default is None).
    inset: "herfindhal" of "boxplot" or None, optional
        inset to include (default is None)
    file_name: str, optional
        output file name (default is "neighborhood_watch_distribution.pdf")
    """

    if inset == "herfindhal":
        herfindhal = True
        bplot = False
    elif inset == "boxplot":
        bplot = True
        herfindhal = False
    else:
        bplot = False
        herfindhal = False

    neigh_watch = nw_helix + nw_sheet + nw_loop
    max_value = max(neigh_watch)

    if not bins:
        bins = range(0, int(max_value + 1) + 1)

    fig, ax1 = plt.subplots(figsize=(10, 10))
    plt.hist(
        nw_helix + nw_sheet + nw_loop,
        bins=bins,
        label="total",
        alpha=0.3,
        histtype="bar",
        rwidth=0.9,
    )
    plt.hist(
        [nw_helix, nw_sheet, nw_loop],
        bins=bins,
        label=["helix", "sheet", "loop"],
        alpha=1,
        histtype="bar",
        rwidth=0.9,
    )
    plt.legend()
    plt.title("Node's average link weight distribution")
    plt.xlabel("w/k")
    plt.ylabel("N(w/k)")

    if bplot:
        left, bottom, width, height = [0.15, 0.65, 0.35, 0.1]
        ax2 = fig.add_axes([left, bottom, width, height])
        ax2.boxplot(neigh_watch, vert=False)
        plt.yticks([])
        plt.xticks(fontsize=18)
        plt.xlabel("w/k", fontsize=18)
        plt.xlim(-1, int(max_value + 1) + 1)

    elif herfindhal:
        #        left, bottom, width, height = [0.15, 0.65, 0.3, 0.13]
        #        ax2 = fig.add_axes([left, bottom, width, height])
        H = herfindhal_index(neigh_watch)
        N = len(neigh_watch)
        t = "$\\newline$".join(
            [
                "",
                "Herfindhal index:",
                "H = %.3f" % (H),
                "1/H = %.1f" % (1 / H),
                "N = %s" % (N),
                "N. outliers = %.3f * N" % (1 - 1 / (H * N)),
            ]
        )
        plt.text(
            0.05,
            0.7,
            t,
            verticalalignment="center",
            horizontalalignment="left",
            fontsize=18,
            transform=ax1.transAxes,
        )

    plt.savefig(os.path.join(folder_path, file_name))
    plt.close()


def plot_nw_sequence(db_2, folder_path, file_name="nw_sequence.pdf", onlyA=False):
    """Plots the histogram of the Nw values along the seuqence.
        
    Parameters
    ----------
    db_2 : DataFrame
        Pandas DataFrame containing node properties
    folder_path: str
        output folder path
    bins: int or list or None, optional
        number of list of bins. In None, bins of size 1 from zero to the max value of Nw are employed (default is None).
    inset: "herfindhal" of "boxplot" or None, optional
        inset to include (default is None)
    file_name: str, optional
        output file name (default is "nw_sequence.pdf")
    onlyA: boolean, optional
        if True, considers only chain A (default is False)
    """

    nw_chain = []
    seq_positions = []
    # select only chain A
    for _, row in db_2.iterrows():
        if onlyA and row["Position"][0] == "B":
            break
        else:
            nw_chain.append(row["Weight/Degree"])
            seq_positions.append(row["Sequence position"])

    colors = ["b" if nw >= 10 and nw <= 16 else "r" for nw in nw_chain]

    plt.figure(figsize=(35, 3.5))
    plt.scatter(seq_positions, nw_chain, c=colors)
    plt.ylabel("w/k")
    plt.xlabel("sequence")
    plt.xticks(seq_positions)
    # plt.yticks([min(nw_chain), 10, 16, max(nw_chain)])
    plt.xlim(min(seq_positions) - 0.5, max(seq_positions) + 0.5)
    # plt.ylim(6, 25)
    plt.grid()
    plt.savefig(os.path.join(folder_path, file_name))
    plt.close()


def plot_wij_distribution(
    net,
    folder_path,
    bins=None,
    normed=True,
    yscale="linear",
    boxplot=False,
    title="Pairwise weight distribution in the amino acids network",
    file_name="pairwise_distrib.pdf"
):
    """Plots the pairwise weights distribution.
        
    Parameters
    ----------
    net : Graph
        NetworkX Graph object
    folder_path: str
        output folder path
    bins: int or list or None, optional
        number of list of bins. In None, bins of size 1 from zero to the max value of Nw are employed (default is None).
    normed: boolean, optional
        if True, normalized the distribution (default is True)
    yscale: str, optional
        type of yscale (linear or logarithmic) (default is linear)
    boxplot: boolean, optional
        if True, add boxplot (default is False)
    tile: str, optional
        title of the plot (default is "Pairwise weight distribution in the amino acids network")
    file_name: str, optional
        output file name (default is "pairwise_distrib.pdf")
    """

    wij_list = []
    for u, v in net.edges:
        wij = net.get_edge_data(u, v)["weight"]
        wij_list.append(wij)

    if not bins:
        bins = range(1, max(wij_list) + 1)

    hist, bin_edges = np.histogram(wij_list, bins=bins, normed=normed)
    xticks = [x + 0.5 for x in bin_edges[0:-1]]

    fig, ax1 = plt.subplots()
    plt.hist(wij_list, bins=range(1, max(wij_list) + 1, 3), normed=normed)
    # plt.scatter(xticks, hist)
    plt.title(title)
    plt.xlabel("w(i,j)")
    if normed:
        plt.ylabel("P(w(i, j))")
    else:
        plt.ylabel("N(w(i, j))")
    # plt.ylim(min(hist), max(hist))
    plt.yscale(yscale)
    if boxplot:
        left, bottom, width, height = [0.38, 0.75, 0.49, 0.1]
        ax2 = fig.add_axes([left, bottom, width, height])
        ax2.boxplot(wij_list, vert=False)
        plt.yticks([])
        plt.xlim(-1, max(wij_list) + 1)
        plt.xlabel('w(i,j)')
    plt.savefig(os.path.join(folder_path, file_name))
    plt.close()


def subnet_highnw(net):
    """Retuns subnetwork of high Nw
        
    Parameters
    ----------
    net : Graph
        NetworkX Graph object
    
    Returns
    -------
    Graph
        NetworkX Graph object
    """

    subnet = net.copy()
    to_remove = []
    for index, node in enumerate(subnet.nodes):
        k = nx.degree(subnet, node)
        w = nx.degree(subnet, node, weight="weight")
        nw = w / k
        if nw < 13:
            to_remove.append(node)
    subnet.remove_nodes_from(to_remove)

    return subnet


def subnet_highwij(net):
    """Returns the subnetwork of high wij
        
    Parameters
    ----------
    net : Graph
        NetworkX Graph object
    
    Returns
    -------
    Graph
        NetworkX Graph object
    """

    subnet = net.copy()
    to_remove = []

    for index, (u, v) in enumerate(subnet.edges):
        wij = subnet.get_edge_data(u, v)["weight"]
        if wij < 30:
            to_remove.append((u, v))
    subnet.remove_edges_from(to_remove)

    # remove nodes of degree zero
    degrees_zero = []
    for node in subnet.nodes:
        k = nx.degree(subnet, node)
        if k == 0:
            degrees_zero.append(node)
    subnet.remove_nodes_from(degrees_zero)

    return subnet


def lcc(net):
    """Returns the largest connected component of a graph
        
    Parameters
    ----------
    net : Graph
        NetworkX Graph object
    
    Returns
    -------
    Graph
        NetworkX Graph object
    """
    lcc = max(nx.connected_component_subgraphs(net), key=len)

    return lcc


def get_edge_data(net):
    """Gets edge data of a graph
        
    Parameters
    ----------
    net : Graph
        NetworkX Graph object
    
    Returns
    -------
    list
        list of weights
    list
        list of colors, based on the weight values (from blue to red increasing weight)
    """

    wij_list = []
    edge_colors_net = []
    for u, v in net.edges:
        wij = net.get_edge_data(u, v)["weight"]
        wij_list.append(wij)
        if wij < 10:
            edge_colors_net.append("blue")
        elif wij < 20:
            edge_colors_net.append("cyan")
        elif wij < 30:
            edge_colors_net.append("greenyellow")
        elif wij < 40:
            edge_colors_net.append("yellow")
        elif wij < 50:
            edge_colors_net.append("orange")
        else:
            edge_colors_net.append("red")

    return wij_list, edge_colors_net


def get_size_from_nw(net):
    """Gets node size from Nw value of nodes of a graph
        
    Parameters
    ----------
    net : Graph
        NetworkX Graph object
    
    Returns
    -------
    list
        list of sizes
    """

    nw_net_size = []
    for index, node in enumerate(net.nodes):
        k = nx.degree(net, node)
        w = nx.degree(net, node, weight="weight")
        nw = w / k
        if nw < 5:
            nw_net_size.append(1 * 500)
        elif nw < 10:
            nw_net_size.append(2 * 500)
        elif nw < 15:
            nw_net_size.append(3 * 500)
        elif nw < 20:
            nw_net_size.append(4 * 500)
        elif nw < 25:
            nw_net_size.append(5 * 500)
        else:
            nw_net_size.append(6 * 500)

    return nw_net_size


def cliques(net):
    """Returns the cliques of a network
        
    Parameters
    ----------
    net : Graph
        NetworkX Graph object
    
    Returns
    -------
    dict
        dictionary of cliques of given size
    """
    cliques_list = list(clique.enumerate_all_cliques(net))

    cliques_dict = dict()
    for c in cliques_list:
        l = len(c)
        try:
            cliques_dict[l].append(c)
        except KeyError:
            cliques_dict[l] = [c]

    return cliques_dict


def get_edge_color(edge, edge_color_dict):
    """Returns the color of an edge
        
    Parameters
    ----------
    edge : tuple
        edge (u, v)
    edge_color_dict: dict
        dictionary of edge colors
    
    Returns
    -------
    str
        edge color
    """

    u, v = edge
    try:
        color = edge_color_dict[(u, v)]
    except KeyError:
        color = edge_color_dict[(v, u)]

    return color


def triangles_ipn(ipn, edge_color_dict):
    """Returns the triangles of an induced perturbation network
        
    Parameters
    ----------
    ipn : Graph
        NetworkX Graph object, an induced perturbation network
    edge_color_dict: dict
        dictionary of edge colors
    
    Returns
    -------
    dict
        dictionary of triangles of given type (rrr, ggg, rrg, ggr)
    """

    cliques_dict = cliques(ipn)
    try:
        triangles = cliques_dict[3]
    except KeyError:
        triangles = []

    triangles_dict = {"ggg": [], "rrr": [], "ggr": [], "rrg": []}

    for (u, v, w) in triangles:

        edge1 = (u, v)
        edge2 = (v, w)
        edge3 = (u, w)

        color1 = get_edge_color(edge1, edge_color_dict)[0]
        color2 = get_edge_color(edge2, edge_color_dict)[0]
        color3 = get_edge_color(edge3, edge_color_dict)[0]

        colors = color1 + color2 + color3

        if colors == "rgr" or colors == "grr":
            colors = "rrg"
        elif colors == "grg" or colors == "rgg":
            colors = "ggr"

        triangles_dict[colors] += [[u, v, w]]

    return triangles_dict
