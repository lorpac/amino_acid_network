import sys
import os
import tools.structure_analysis_tools as sa
import tools.network_visualization as nv
import tools.helpers as helpers
import networkx as nx
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import csv
from shutil import copyfile

plt.rc("text", usetex=True)

data_folder = "data"

# Read configuration file
config_file_path = "analysis_config.json"
parameters = helpers.get_configuration_parameters(config_file_path)

name = parameters["name"]
pdb_id = parameters["pdb_id"]
cutoff = parameters["cutoff"]
dim = parameters["dim"]
if dim == "all":
    dim = ""
rel_list = sa.list_relations(dim)
pdbs_path = parameters["pdbs_path"]
select_positions = parameters["select_positions"]
if select_positions:
    start = parameters["start"]
    stop = parameters["stop"]
    selected_positions = range(start, stop)
else:
    selected_positions = None
variant = parameters["variant"]
if variant:
    reference_folder = parameters["reference_folder"]
    threshold = parameters["threshold"]
    single_mutation = parameters["single_mutation"]
    if single_mutation:
        source = parameters["source"]
subnet_highnw = parameters["subnet_highnw"]
subnet_highwij = parameters["subnet_highwij"]
draw_neighborhoods = parameters["draw_neighborhoods"]
calculate_cliques = parameters["calculate_cliques"]
if calculate_cliques:
    draw_large_cliques = parameters["draw_large_cliques"]
remove_hydrogen_atoms = parameters["remove_hydrogen_atoms"]

results_to_report = []

folder_path = os.path.join("results", name, dim + "analysis")
os.makedirs(folder_path, exist_ok=True)

helpers.save_configuration_parameters(config_file_path, folder_path)

for k, v in parameters.items():
    print(k, "\t", v)

# import mutations
(
    mutations_type,
    pathogenic,
    non_pathogenic,
    both,
    pathogenic_names,
    non_pathogenic_names,
) = sa.import_mutations(folder_path)

# load k_w dictionary (from database)
k_w = pickle.load(open(os.path.join(data_folder, "k_w.p"), "rb"))

# create aa network and save database of its aa (2 formats)
net, db_1, db_2, mol, downloaded = sa.create_aa_network(
    pdb_id,
    rel_list,
    folder_path,
    selected_positions,
    cutoff=cutoff,
    kw_reprod=True,
    k_w=k_w,
    pdbs_path=pdbs_path,
    remove_hydrogen_atoms=remove_hydrogen_atoms,
    db_1_name=pdb_id + 'nodes_1.csv',
    db_2_name=pdb_id + 'nodes_2.csv'
)

print("AAN created")
print("N. nodes: ", len(net.nodes))
print("N. edges: ", len(net.edges))

contact_list = list(net.edges)
file_contacts = os.path.join(folder_path, pdb_id + "contact_list.csv")
with open(file_contacts, "w", newline="") as f:
    writer = csv.writer(f)
    for line in contact_list:
        writer.writerow(line)



results_to_report.append(["AAN"])
results_to_report.append(["N. nodes: ", len(net.nodes)])
results_to_report.append(["N. edges: ", len(net.edges)])

chains = sa.get_chains(mol)

# plot nw distribution
neigh_watch = db_2["Weight/Degree"].values

nw_helix, nw_sheet, nw_loop = sa.divided_nw_lists(db_2)

sa.plot_nw_hist(
    nw_helix,
    nw_sheet,
    nw_loop,
    folder_path,
    inset="herfindhal",
    file_name="neighborhood_watch_distribution1.pdf",
)

sa.plot_nw_hist(
    nw_helix,
    nw_sheet,
    nw_loop,
    folder_path,
    inset="boxplot",
    file_name="neighborhood_watch_distribution2.pdf",
)

# plot nw along the sequence

if name[0:3] == "ttr":
    onlyA = True
else:
    onlyA = False

sa.plot_nw_sequence(db_2, folder_path, onlyA=onlyA)

# plot wij distribution

sa.plot_wij_distribution(net, folder_path)

sa.plot_wij_distribution(
    net, folder_path, file_name="pairwise_distrib_log.pdf", yscale="log",
    boxplot=True, normed=False
)


# DRAW PROTEIN NETWORK
pict_folder = os.path.join(folder_path, "network_pictures")
os.makedirs(pict_folder, exist_ok=True)

if downloaded:
    pdb_id = "pdb" + pdb_id

node_labels, sizes, color_map, _ = nv.network_visualization_parameters(
    pdb_id,
    database=db_2,
    net=net,
    rel_list=rel_list,
    kw_reprod=True,
    k_w=k_w,
    colors="neighborhood_watch_sharp",
    sizes=None,
)
node_labels2, sizes2, color_map2, _ = nv.network_visualization_parameters(
    pdb_id,
    database=db_2,
    net=net,
    rel_list=rel_list,
    kw_reprod=True,
    k_w=k_w,
    colors=None,
    sizes=None,
)
node_labels3, sizes3, color_map3, edge_color_map = nv.network_visualization_parameters(
    pdb_id,
    database=db_2,
    net=net,
    rel_list=rel_list,
    kw_reprod=True,
    k_w=k_w,
    colors="pairwise_sharp",
    sizes="neighborhood_watch_sharp",
)


pos_spring = nx.spring_layout(net)
pos_circular = nx.circular_layout(net)

name1 = os.path.join(pict_folder, pdb_id + "_color")
name2 = os.path.join(pict_folder, pdb_id + "_size")
name3 = os.path.join(pict_folder, pdb_id + "_colorsize")
name3b = os.path.join(pict_folder, pdb_id + "_colornosize")
name4 = os.path.join(pict_folder, pdb_id + "_nodes")
name5 = os.path.join(pict_folder, pdb_id + "_links")
name6 = os.path.join(pict_folder, pdb_id + "_links_color")

nv.draw_network(
    net,
    node_labels,
    color_map=color_map,
    draw_edges=False,
    draw_labels=False,
    name=name1,
    figsize=(50, 50),
    pos=pos_spring,
)
nv.draw_network(
    net,
    node_labels,
    sizes=sizes,
    draw_labels=False,
    name=name2,
    figsize=(50, 50),
    pos=pos_spring,
)
nv.draw_network(
    net,
    node_labels,
    sizes=sizes,
    color_map=color_map,
    draw_labels=False,
    name=name3,
    figsize=(50, 50),
    pos=pos_spring
)
nv.draw_network(
    net,
    node_labels,
    color_map=color_map,
    draw_labels=False,
    name=name3b,
    figsize=(50, 50),
    pos=pos_spring,
    edge_thickness=False
)

nv.draw_network(
    net,
    node_labels2,
    color_map="deepskyblue",
    draw_edges=False,
    draw_labels=False,
    name=name4,
    figsize=(50, 50),
    pos=pos_spring
)
nv.draw_network(
    net,
    node_labels2,
    sizes=sizes2,
    color_map="deepskyblue",
    draw_labels=False,
    name=name5,
    figsize=(50, 50),
    pos=pos_spring,
    legend_thickness=True
)
nv.draw_network(
    net,
    node_labels3,
    sizes=sizes3,
    edge_color_map=edge_color_map,
    draw_labels=False,
    name=name6,
    figsize=(50, 50),
    pos=pos_spring,
    facecolor="silver",
)

net_path = os.path.join(folder_path, dim)

pickle.dump(net, open(net_path + "net.p", "wb"))

nv.plot_legends(folder_path=folder_path)

# sizes_dict = {node: size_map[i] for i, node in enumerate(net.nodes)}
edge_color_dict = {edge: edge_color_map[i] for i, edge in enumerate(net.edges)}
color_dict = {node: color_map[i] for i, node in enumerate(net.nodes)}

# CLIQUES
if calculate_cliques:
    cliques = sa.cliques(net)
    print("Cliques of AAN created")
    results_to_report.append(["Cliques of AAN"])

    if draw_large_cliques:
        cliques_folder = os.path.join(pict_folder, "large_cliques")
        os.makedirs(cliques_folder, exist_ok=True)

    for l in cliques:
        c_list = cliques[l]
        n = len(c_list)
        print(n, " cliques of size ", l)
        results_to_report.append(["%s cliques of size %s" % (n, l)])
        if draw_large_cliques and l > 4:
            for j, c in enumerate(c_list):
                c_net = nx.Graph()
                c_net.add_nodes_from(c)
                for i, v in enumerate(c[:-1]):
                    for u in c[i + 1 : :]:
                        wuv = net.get_edge_data(v, u)["weight"]
                        c_net.add_edge(u, v, weight=wuv)

                color_map_c = []
                edge_color_map_c = []
                node_labels_c = {}
                for v in c_net.nodes:
                    color = color_dict[v]
                    label = node_labels[v]
                    node_labels_c[v] = label
                    color_map_c.append(color)
                for u, v in c_net.edges:
                    edge_color = edge_color_dict[(u, v)]
                    edge_color_map_c.append(edge_color)

                file_name = os.path.join(
                    cliques_folder, "%sclique%s-%s" % (pdb_id, l, j)
                )
                nv.draw_network(
                    c_net,
                    node_labels_c,
                    edge_color_map=edge_color_map_c,
                    name=file_name,
                    figsize=(10, 10),
                    node_edge_size=1,
                    node_edge_colors="gray",
                )


# DRAW NEIGHBORHOODS
if draw_neighborhoods:
    neigh_path = os.path.join(pict_folder, "neighborhoods")
    os.makedirs(neigh_path, exist_ok=True)

    for node in net.nodes:
        if node[0] != "A":
            break
        file_name = os.path.join(neigh_path, pdb_id + node)
        nv.draw_neighborhood(
            net,
            node,
            node_labels,
            pathogenic,
            non_pathogenic,
            both,
            save_fig=True,
            file_name=file_name,
            threshold_edge=20,
        )

if variant:
    results_perturbation = []
    # DRAW PERTURBATION NETWORK
    #    path_ref = os.path.join(data_dir, reference_folder)
    net_ref_path = os.path.join(
        reference_folder, dim + "analysis", dim + "net.p"
    )
    net_ref = pickle.load(open(net_ref_path, "rb"))
    db_ref_path = os.path.join(reference_folder, "analysis", "database_pos_2.csv")
    db_ref = pd.DataFrame(pd.read_csv(db_ref_path))

    (
        net_p,
        node_labels_p,
        size_map_p,
        color_map_p,
        edge_color_dict_p,
        node_borders_p,
        edges_1only,
        edges_2only,
        common_edges,
    ) = nv.create_perturbation_network(
        net_ref, net, db_ref, db_2, threshold=threshold
    )

    print("Perturbation network")
    print("Threshold: ", threshold)
    print("N. nodes in PN: ", len(net_p.node))
    print("N. edges in PN: ", len(net_p.edges))
    n_components = nx.algorithms.components.number_connected_components(net_p)
    print("N. connected components in PN: ", n_components)

    results_perturbation.append(["Perturbation network created"])
    results_perturbation.append(["Threshold: ", threshold])
    results_perturbation.append(["N. nodes in PN: ", len(net_p.node)])
    results_perturbation.append(["N. edges in PN: ", len(net_p.edges)])
    results_perturbation.append(["N. connected components in PN: ", n_components])

    edge_color_map_p = []
    for u, v in net_p.edges:
        try:
            edge_color_map_p.append(edge_color_dict_p[u, v])
        except KeyError:
            edge_color_map_p.append(edge_color_dict_p[v, u])

    nodes_pos = nv.nodes_position_perturbation(
        net_p, dim, chains, scale=4, max_x=30
    )
    nodes_pos2 = nv.nodes_position_perturbation(
        net_p, dim, chains, columns=True
    )

    file_name = os.path.join(
        pict_folder, pdb_id + "perturbation_net_links_color" + str(threshold)
    )
    file_name2 = os.path.join(
        pict_folder,
        pdb_id + "perturbation_net_links_color" + str(threshold) + "_columns",
    )

    plt.rc("text", usetex=False)

    nv.draw_network(
        net_p,
        node_labels_p,
        sizes=size_map_p,
        edge_color_map=edge_color_map_p,
        color_map=color_map_p,
        draw_labels=True,
        name=file_name,
        figsize=(20, 15),
        node_edge_colors=node_borders_p,
        pos=nodes_pos,
    )

    nv.draw_network(
        net_p,
        node_labels_p,
        sizes=size_map_p,
        edge_color_map=edge_color_map_p,
        color_map=color_map_p,
        draw_labels=True,
        name=file_name2,
        figsize=(20, 30),
        node_edge_colors=node_borders_p,
        pos=nodes_pos2,
    )

    if draw_neighborhoods:
        neigh_pert_path = os.path.join(pict_folder, "neighborhoods_perturbation")
        os.makedirs(neigh_pert_path, exist_ok=True)
        for node in net_p.nodes:
            if node[0] != "A":
                continue
            file_name = os.path.join(neigh_pert_path, pdb_id + node)
            nv.draw_neighborhood_perturbation(
                net_p,
                node,
                node_labels_p,
                size_map_p,
                color_map_p,
                edge_color_map_p,
                node_borders_p,
                file_name=file_name,
            )

    sizes_dict = {node: size_map_p[i] for i, node in enumerate(net_p.nodes)}
    edge_color_dict = {edge: edge_color_map_p[i] for i, edge in enumerate(net_p.edges)}
    color_dict = {node: color_map_p[i] for i, node in enumerate(net_p.nodes)}
    node_borders_dict = {node: node_borders_p[i] for i, node in enumerate(net_p.nodes)}
    weights_p = nx.get_edge_attributes(net_p, "weight")

    net_p_file = os.path.join(folder_path, "net_p%s.csv" % (threshold))
    with open(net_p_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["u", "v", "weight", "color"])
        for i, (u, v) in enumerate(net_p.edges):
            w = net_p.get_edge_data(u, v)["weight"]
            color = edge_color_map_p[i]
            writer.writerow([u, v, w, color])

    # INDUCED PERTURBATION NETWORK

    if single_mutation:

        ipn = nv.induced_perturbation_network(net_p, source)

        if ipn:
            node_labels_ipn = {node: node_labels_p[node] for node in ipn.nodes}

            sizes_ipn = [sizes_dict[node] for node in ipn.nodes]

            edge_color_ipn = []
            weights_ipn = {}
            for u, v in ipn.edges:
                try:
                    edge_color_ipn.append(edge_color_dict[(u, v)])
                except KeyError:
                    edge_color_ipn.append(edge_color_dict[(v, u)])
                try:
                    weights_ipn[(u, v)] = weights_p[(u, v)]
                except KeyError:
                    weights_ipn[(u, v)] = weights_p[(v, u)]

            nx.set_edge_attributes(ipn, weights_ipn, name="weight")

            print("Induced perturbation network created")
            print("N. nodes in IPN: ", len(ipn.node))
            print("N. edges in IPN: ", len(ipn.edges))

            results_perturbation.append(["Induced perturbation network"])
            results_perturbation.append(["N. nodes in IPN: ", len(ipn.node)])
            results_perturbation.append(["N. edges in IPN: ", len(ipn.edges)])

            color_ipn = [color_dict[node] for node in ipn.nodes]
            node_borders_ipn = [node_borders_dict[node] for node in ipn.nodes]

            pickle.dump(
                [
                    ipn,
                    node_labels_ipn,
                    sizes,
                    edge_color_ipn,
                    color_ipn,
                    node_borders_ipn,
                ],
                open(os.path.join(folder_path, "ipn.p"), "wb"),
            )

            pos_ipn = nv.nodes_position_perturbation(
                ipn, dim, chains
            )

            pos_ipn2 = nv.nodes_position_perturbation(
                ipn, dim, chains, columns=True
            )

            node_borders_ipn = ["r"]
            node_borders_ipn = node_borders_ipn + [
                "gray" for node in list(ipn.nodes)[1::]
            ]

            color_ipn = ["y"]
            color_ipn = color_ipn + ["lightgray" for node in list(ipn.nodes)[1::]]
            color_ipn_dict = {node: color_ipn[i] for i, node in enumerate(ipn.node)}

            name_ipn = os.path.join(
                pict_folder,
                pdb_id + "perturbation_ipn_" + str(threshold) + source + "_" + dim,
            )
            nv.draw_network(
                ipn,
                node_labels_ipn,
                sizes=600,
                edge_color_map=edge_color_ipn,
                color_map=color_ipn,
                draw_labels=True,
                name=name_ipn,
                figsize=(30, 30),
                pos=pos_ipn,
                labels_size=18,
                node_edge_size=1,
            )

            name_ipn2 = os.path.join(
                pict_folder,
                pdb_id
                + "perturbation_ipn_"
                + str(threshold)
                + source
                + "_"
                + dim
                + "columns",
            )
            nv.draw_network(
                ipn,
                node_labels_ipn,
                sizes=600,
                edge_color_map=edge_color_ipn,
                color_map=color_ipn,
                draw_labels=True,
                name=name_ipn2,
                figsize=(30, 30),
                pos=pos_ipn2,
                labels_size=18,
                node_edge_size=1,
            )

            ipn_file = os.path.join(folder_path, "ipn%s.csv" % (threshold))
            with open(ipn_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["u", "v", "weight", "color"])
                for i, (u, v) in enumerate(ipn.edges):
                    w = ipn.get_edge_data(u, v)["weight"]
                    color = edge_color_ipn[i]
                    writer.writerow([u, v, w, color])

            triangles_ipn = sa.triangles_ipn(ipn, edge_color_dict)

            green_tr = triangles_ipn["ggg"]
            red_tr = triangles_ipn["rrr"]
            mostlygreen_tr = triangles_ipn["ggr"]
            mostlyred_tr = triangles_ipn["rrg"]

            print("N. green triangles: ", len(green_tr))
            print("N. mostly green triangles: ", len(mostlygreen_tr))
            print("N. red triangles: ", len(red_tr))
            print("N. mostly red triangles: ", len(mostlyred_tr))

            results_perturbation.append(["N. green triangles: ", len(green_tr)])
            results_perturbation.append(
                ["N. mostly green triangles: ", len(mostlygreen_tr)]
            )
            results_perturbation.append(["N. red triangles: ", len(red_tr)])
            results_perturbation.append(
                ["N. mostly red triangles: ", len(mostlyred_tr)]
            )

            # Green and red IPN
            green_ipn = nx.Graph()
            red_ipn = nx.Graph()
            for i, (u, v) in enumerate(ipn.edges):
                color = edge_color_ipn[i]
                wuv = ipn.get_edge_data(u, v)["weight"]
                if color == "green":
                    green_ipn.add_edge(u, v, weight=wuv)
                else:
                    red_ipn.add_edge(u, v, weight=wuv)

            color_green_ipn = [color_ipn_dict[node] for node in green_ipn.nodes]
            color_red_ipn = [color_ipn_dict[node] for node in red_ipn.nodes]

            name_green = os.path.join(
                pict_folder,
                pdb_id
                + "perturbation_green_ipn_"
                + str(threshold)
                + source
                + "_"
                + dim,
            )

            name_red = os.path.join(
                pict_folder,
                pdb_id + "perturbation_red_ipn_" + str(threshold) + source + "_" + dim,
            )

            nv.draw_network(
                green_ipn,
                node_labels_ipn,
                sizes=600,
                edge_color_map="green",
                color_map=color_green_ipn,
                draw_labels=True,
                name=name_green,
                figsize=(10, 10),
                pos=pos_ipn,
                labels_size=18,
                node_edge_size=1,
            )

            nv.draw_network(
                red_ipn,
                node_labels_ipn,
                sizes=600,
                edge_color_map="red",
                color_map=color_red_ipn,
                draw_labels=True,
                name=name_red,
                figsize=(10, 10),
                pos=pos_ipn,
                labels_size=18,
                node_edge_size=1,
            )

            print("Green IPN created")
            print("N. nodes in green IPN: ", len(green_ipn))
            print("N. edges in green IPN: ", len(green_ipn.edges))
            n_components = nx.algorithms.components.number_connected_components(
                green_ipn
            )
            print("N. connected components in green IPN: ", n_components)

            results_perturbation.append(["Green IPN"])
            results_perturbation.append(["N. nodes in green IPN: ", len(green_ipn)])
            results_perturbation.append(
                ["N.edges in green IPN: ", len(green_ipn.edges)]
            )
            results_perturbation.append(
                ["N. connected components in green IPN: ", n_components]
            )

            # CLIQUES
            green_cliques = sa.cliques(green_ipn)
            pickle.dump(
                green_cliques, open(net_path + "green_cliques%s.p" % (threshold), "wb")
            )
            print("Cliques of green IPN created")
            results_perturbation.append(["Cliques of green IPN created"])

            for l in green_cliques:
                c_list = green_cliques[l]
                n = len(c_list)
                print(n, " cliques of size ", l)
                results_perturbation.append(["%s cliques of size %s" % (n, l)])
                if l > 5:
                    green_cliques_folder = os.path.join(
                        pict_folder, "large_cliques_greenIPN"
                    )
                    os.makedirs(green_cliques_folder, exist_ok=True)
                    for j, c in enumerate(c_list):
                        c_net = nx.Graph()
                        c_net.add_nodes_from(c)
                        for i, v in enumerate(c[:-1]):
                            for u in c[i + 1 : :]:
                                wuv = net.get_edge_data(v, u)["weight"]
                                c_net.add_edge(u, v, weight=wuv)

                        color_map_c = []
                        edge_color_map_c = []
                        node_labels_c = {}
                        for v in c_net.nodes:
                            color = color_dict[v]
                            label = node_labels[v]
                            node_labels_c[v] = label
                            color_map_c.append(color)
                        for u, v in c_net.edges:
                            edge_color = edge_color_dict[(u, v)]
                            edge_color_map_c.append(edge_color)

                        file_name = os.path.join(
                            cliques_folder, "%sclique%s-%s" % (pdb_id, l, j)
                        )
                        nv.draw_network(
                            c_net,
                            node_labels_c,
                            edge_color_map=edge_color_map_c,
                            name=file_name,
                            figsize=(10, 10),
                            node_edge_size=1,
                            node_edge_colors="gray",
                        )

            print("Red IPN created")
            print("N. nodes in red IPN: ", len(red_ipn))
            print("N.edges in red IPN: ", len(red_ipn.edges))
            n_components = nx.algorithms.components.number_connected_components(red_ipn)
            print("N. connected components in red IPN: ", n_components)

            results_perturbation.append(["Red IPN"])
            results_perturbation.append(["N. nodes in red IPN: ", len(red_ipn)])
            results_perturbation.append(["N.edges in red IPN: ", len(red_ipn.edges)])
            results_perturbation.append(
                ["N. connected components in red IPN: ", n_components]
            )

            # CLIQUES
            red_cliques = sa.cliques(red_ipn)
            pickle.dump(
                red_cliques, open(net_path + "red_cliques%s.p" % (threshold), "wb")
            )
            print("Cliques of red IPN created")
            results_perturbation.append(["Cliques of red IPN created"])

            for l in red_cliques:
                c_list = red_cliques[l]
                n = len(c_list)
                print(n, " cliques of size ", l)
                results_perturbation.append(["%s cliques of size %s" % (n, l)])
                if l > 5:
                    red_cliques_folder = os.path.join(
                        pict_folder, "large_cliques_redIPN"
                    )
                    os.makedirs(red_cliques_folder, exist_ok=True)
                    for j, c in enumerate(c_list):
                        c_net = nx.Graph()
                        c_net.add_nodes_from(c)
                        for i, v in enumerate(c[:-1]):
                            for u in c[i + 1 : :]:
                                wuv = net.get_edge_data(v, u)["weight"]
                                c_net.add_edge(u, v, weight=wuv)

                        color_map_c = []
                        edge_color_map_c = []
                        node_labels_c = {}
                        for v in c_net.nodes:
                            color = color_dict[v]
                            label = node_labels[v]
                            node_labels_c[v] = label
                            color_map_c.append(color)
                        for u, v in c_net.edges:
                            edge_color = edge_color_dict[(u, v)]
                            edge_color_map_c.append(edge_color)

                        file_name = os.path.join(
                            cliques_folder, "%sclique%s-%s" % (pdb_id, l, j)
                        )
                        nv.draw_network(
                            c_net,
                            node_labels_c,
                            edge_color_map=edge_color_map_c,
                            name=file_name,
                            figsize=(10, 10),
                            node_edge_size=1,
                            node_edge_colors="gray",
                        )

            file_results_perturbation = os.path.join(
                folder_path, "results_perturbation%s.csv" % (threshold)
            )
            with open(file_results_perturbation, "w", newline="") as f:
                writer = csv.writer(f)
                for line in results_perturbation:
                    writer.writerow(line)


if subnet_highnw:
    # SUBNETWORK OF HIGH NW
    subnet = sa.subnet_highnw(net)

    print("Subnetwork of high nw created")
    N_subnet = len(subnet)
    print("N. nodes in subnet (before removing degrees zero): ", N_subnet)
    print("N. edges in subnet: ", len(subnet.edges))
    n_components = nx.algorithms.components.number_connected_components(subnet)
    print("N. connected components in subnet: ", n_components)

    results_to_report.append(["Subnetwork of high nw created"])
    results_to_report.append(
        ["N. nodes in subnet (before removing degrees zero): ", len(subnet)]
    )
    results_to_report.append(["N. edges in subnet: ", len(subnet.edges)])
    results_to_report.append(["N. connected components in subnet: ", n_components])

    # remove nodes of degree zero
    degrees_zero = []
    for node in subnet.nodes:
        k = nx.degree(subnet, node)
        if k == 0:
            degrees_zero.append(node)
    subnet.remove_nodes_from(degrees_zero)

    print("N. nodes in subnet (after removing degrees zero): ", len(subnet))
    print("N. components of size 1: ", N_subnet - len(subnet))

    results_to_report.append(
        ["N. nodes in subnet (after removing degrees zero): ", len(subnet)]
    )
    results_to_report.append(["N. components of size 1: ", N_subnet - len(subnet)])

    # plot P(wij) subnet
    wij_list_sn, edge_colors_sn = sa.get_edge_data(subnet)

    file_name = os.path.join(folder_path, "pairwise_distrib_highnw.pdf")
    title = "Pairwise weight distribution in the sub-network with $w/k\geq13$"
    sa.plot_wij_distribution(
        subnet, folder_path, title=title, file_name=file_name, yscale="log"
    )

    # plot subnet
    node_labels_subnet = {node: node_labels3[node] for node in subnet.nodes}
    nw_subnet_size = sa.get_size_from_nw(subnet)

    nv.draw_network(
        subnet,
        node_labels_subnet,
        sizes=nw_subnet_size,
        edge_color_map=edge_colors_sn,
        draw_labels=False,
        name=pict_folder + pdb_id + "_links_color" + "_highnw",
        figsize=(50, 50),
        pos=pos_spring,
        facecolor="silver",
    )

    # largest connected component
    lcc = sa.lcc(subnet)

    print("LCC of subnetwork of high nw created")
    print("N. nodes in LCC of subnet: ", len(lcc))

    results_to_report.append(["LCC of subnetwork of high nw"])
    results_to_report.append(["N. nodes in LCC of subnet: ", len(lcc)])

    # plot P(wij) of lcc of subnet

    wij_list_lcc, edge_colors_lcc = sa.get_edge_data(lcc)

    file_name = os.path.join(folder_path, "pairwise_distrib_highnwLCC.pdf")
    title = (
        "Pairwise weight distribution in the LCC of the sub-network with $w/k\geq13$"
    )
    sa.plot_wij_distribution(
        lcc, folder_path, title=title, file_name=file_name, yscale="log"
    )

    # plot lcc
    node_labels_lcc = {node: node_labels3[node] for node in lcc.nodes}
    nw_lcc_size = sa.get_size_from_nw(lcc)

    name = os.path.join(pict_folder, pdb_id + "_links_color" + "_highnwLCC")

    nv.draw_network(
        lcc,
        node_labels_lcc,
        sizes=nw_lcc_size,
        edge_color_map=edge_colors_lcc,
        draw_labels=False,
        name=name,
        figsize=(50, 50),
        pos=pos_spring,
        facecolor="silver",
    )

if subnet_highwij:
    # SUBNETWORK OF HIGH wij
    subnet2 = sa.subnet_highwij(net)

    print("Subnetwork of high wij created")
    N_subnet = len(subnet2)
    print("N. nodes in subnet: ", N_subnet)
    print("N. edges in subnet: ", len(subnet2.edges))
    n_components = nx.algorithms.components.number_connected_components(subnet2)
    print("N. connected components in subnet: ", n_components)

    results_to_report.append(["Subnetwork of high wij created"])
    results_to_report.append(
        ["N. nodes in subnet (before removing degrees zero): ", len(subnet2)]
    )
    results_to_report.append(["N. edges in subnet: ", len(subnet2.edges)])
    results_to_report.append(["N. connected components in subnet: ", n_components])

    wij_list_sn2, edge_colors_sn2 = sa.get_edge_data(subnet2)

    # plot subnet
    node_labels_subnet2 = {node: node_labels3[node] for node in subnet2.nodes}
    nw_subnet_size2 = sa.get_size_from_nw(subnet2)

    name = os.path.join(pict_folder, pdb_id + "_links_color" + "_highwij")

    # nv.draw_network(
    #     subnet2,
    #     node_labels_subnet2,
    #     sizes=nw_subnet_size2,
    #     edge_color_map=edge_colors_sn2,
    #     draw_labels=False,
    #     name=name,
    #     figsize=(50, 50),
    #     pos=pos_spring,
    #     facecolor="silver",
    # )

    nv.draw_network(
        subnet2,
        node_labels_subnet2,
        draw_labels=False,
        name=name,
        figsize=(50, 50),
        pos=pos_spring,
        facecolor="silver",
    )

    # largest connected component
    lcc2 = sa.lcc(subnet2)

    print("LCC of subnetwork of high wij created")
    print("N. nodes in LCC of subnet: ", len(lcc2))

    results_to_report.append(["LCC of subnetwork of high wij created"])
    results_to_report.append(["N. nodes in LCC of subnet: ", len(lcc2)])


file_results = os.path.join(folder_path, "results.csv")
with open(file_results, "w", newline="") as f:
    writer = csv.writer(f)
    for line in results_to_report:
        writer.writerow(line)

