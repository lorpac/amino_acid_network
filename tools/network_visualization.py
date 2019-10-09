import sys
import os
import numpy as np
import networkx as nx
import biographs as bg
import matplotlib.pyplot as plt
import tools.amino_acids_conversion as aac
import pandas as pd
from collections import deque
import matplotlib.patches as patches
from copy import deepcopy


def network_visualization_parameters(
    pdb_id,
    database=None,
    net=None,
    rel_list=["1D", "2D", "3D", "4D"],
    selected_positions=None,
    cutoff=5,
    kw_reprod=False,
    k_w=None,
    pdbs_path=None,
    database_path=None,
    db_1_name=None,
    db_2_name=None,
    separate_jaccard=False,
    separate_weights=False,
    save_csv=True,
    pathogenic=[],
    non_pathogenic=[],
    both=[],
    colors="mutations",
    sizes="neighborhood_watch_sharp",
    downloaded = False
):

    # TO TEST: create database and/or net if not given
    # TODO: create object to define database configuration
    if not isinstance(database, pd.DataFrame) or not net:
        net, _, db_2, _, downloaded = sa.create_aa_network(
            pdb_id,
            rel_list,
            folder_path=database_path,
            selected_positions=selected_positions,
            cutoff=cutoff,
            kw_reprod=kw_reprod,
            k_w=k_w,
            db_1_name=db_1_name,
            db_2_name=db_2_name,
            separate_jaccard=separate_jaccard,
            separate_weights=separate_weights,
            pdbs_path=pdbs_path,
            save_csv=save_csv,
        )

    node_labels = {}
    for node in net.nodes:
        if downloaded:
            info = database[database["Residue name"] == "pdb" + pdb_id + node]
        else:
            info = database[database["Residue name"] == pdb_id + node]
        if len(info) > 1:
            info = info.iloc[0]  # check why more than one
        if type(info) == pd.core.series.Series:
            info = pd.DataFrame(info).transpose()  # it happens when multiple aa
        # are reported for one position
        # --CHECK WHY
        type_aa = aac.three2one(info["Type of residue"].item())
        label = type_aa + node[1::] + ":" + node[0]
        node_labels[node] = label

    mutation_type = []
    neighborhood_watch_sharp = []
    neighborhood_watch_smooth = []
    degree = []
    pairwise_sharp = []

    for node in net.nodes:
        if colors == "mutations" or sizes == "mutations":
            seq_pos = node[1::]
            if seq_pos in pathogenic:
                mutation_type.append("tomato")
            elif seq_pos in non_pathogenic:
                mutation_type.append("limegreen")
            elif seq_pos in both:
                mutation_type.append("gold")
            else:
                mutation_type.append("lightgrey")

        elif colors or sizes:
            k = nx.degree(net, node)
            degree.append(k)
            weight = nx.degree(net, node, weight="weight")
            if colors == "neighborhood_watch_sharp":
                if weight / k < 5:
                    neighborhood_watch_sharp.append("blue")
                elif weight / k < 10:
                    neighborhood_watch_sharp.append("cyan")
                elif weight / k < 15:
                    neighborhood_watch_sharp.append("greenyellow")
                elif weight / k < 20:
                    neighborhood_watch_sharp.append("yellow")
                elif weight / k < 25:
                    neighborhood_watch_sharp.append("orange")
                else:
                    neighborhood_watch_sharp.append("red")

            elif sizes == "neighborhood_watch_sharp":
                if weight / k < 5:
                    neighborhood_watch_sharp.append(1 * 500)
                elif weight / k < 10:
                    neighborhood_watch_sharp.append(2 * 500)
                elif weight / k < 15:
                    neighborhood_watch_sharp.append(3 * 500)
                elif weight / k < 20:
                    neighborhood_watch_sharp.append(4 * 500)
                elif weight / k < 25:
                    neighborhood_watch_sharp.append(5 * 500)
                else:
                    neighborhood_watch_sharp.append(6 * 500)

            elif (
                colors == "neighborhood_watch_smooth"
                or sizes == "neighborhood_watch_smooth"
            ):
                neighborhood_watch_smooth.append(weight / k)

    if colors == "pairwise_sharp" or sizes == "pairwise_sharp":
        for u, v in net.edges:
            wij = net.get_edge_data(u, v)["weight"]
            if wij < 10:
                pairwise_sharp.append("blue")
            elif wij < 20:
                pairwise_sharp.append("cyan")
            elif wij < 30:
                pairwise_sharp.append("greenyellow")
            elif wij < 40:
                pairwise_sharp.append("yellow")
            elif wij < 50:
                pairwise_sharp.append("orange")
            else:
                pairwise_sharp.append("red")

    color_map = []
    edge_color_map = []
    if colors == "mutations":
        color_map = mutation_type
    elif colors == "degree":
        color_map = degree
    elif colors == "neighborhood_watch_sharp":
        color_map = neighborhood_watch_sharp
    elif colors == "neighborhood_watch_smooth":
        color_map = neighborhood_watch_smooth
    elif colors == "pairwise_sharp":
        edge_color_map = pairwise_sharp

    size_map = []
    #    edge_size_map = []
    if sizes == "mutations":
        size_map = mutation_type
    elif sizes == "degree":
        size_map = degree
    elif sizes == "neighborhood_watch_sharp":
        size_map = neighborhood_watch_sharp
    elif sizes == "neighborhood_watch_smooth":
        size_map = neighborhood_watch_smooth

    return node_labels, size_map, color_map, edge_color_map


def draw_network(
    net,
    node_labels,
    sizes=None,
    color_map=None,
    cmap=None,
    edge_color_map=None,
    name="",
    draw_edges=True,
    edge_thickness=True,
    draw_labels=True,
    pos=None,
    layout="spring",
    default_size=2000,
    figsize=(130, 100),
    facecolor="w",
    node_edge_colors=None,
    labels_size=30,
    node_edge_size=4,
    pos_label=None,
    legend_thickness=False
):
    net_copy = deepcopy(net)
    pos_copy = deepcopy(pos)
    if edge_color_map:
        edge_color_map_copy = deepcopy(edge_color_map)
    
    fig = plt.figure(figsize=figsize, facecolor=facecolor)

    if not pos_copy:
        if layout == "circular":
            pos_copy = nx.circular_layout(net)
        else:  # TODO: put all cases
            pos_copy = nx.spring_layout(net)

    
    if legend_thickness:
        samples = range(10, 70, 10)
        maxposx = max(np.array(list(pos_copy.values()))[:, 0])
        maxposy = max(np.array(list(pos_copy.values()))[:, 1])
        for i, s in enumerate(samples):
            net_copy.add_node('sample_%sL' %s)
            net_copy.add_node('sample_%sR' %s)
            pos_copy['sample_%sL' %s] = [maxposx +.1 , maxposy - .1 * i]
            pos_copy['sample_%sR' %s] = [maxposx + .3, maxposy - .1 * i]
            net_copy.add_edge('sample_%sL' %s, 'sample_%sR' %s, weight=s)
            if edge_color_map:
                if s < 10:
                    edge_color_map_copy.append("blue")
                elif s < 20:
                    edge_color_map_copy.append("cyan")
                elif s < 30:
                    edge_color_map_copy.append("greenyellow")
                elif s < 40:
                    edge_color_map_copy.append("yellow")
                elif s < 50:
                    edge_color_map_copy.append("orange")
                else:
                    edge_color_map_copy.append("red")
            # TODO: fix for cases when colormap, color of edges, sizes of nodes etc are not None


    width = nx.get_edge_attributes(net_copy, "weight")

    if not sizes:
        sizes = default_size
    if not color_map:
        color_map = "w"
    if not edge_color_map:
        edge_color_map_copy = "grey"
    else:
        edge_color_map_copy = edge_color_map

    if not node_edge_colors:
        node_edge_colors = ["k" for node in net_copy.nodes]

    nx.draw_networkx_nodes(
        net_copy,
        pos_copy,
        node_size=sizes,
        node_color=color_map,
        linewidths=node_edge_size,
        cmap=cmap,
        edgecolors=node_edge_colors,
    )
    if draw_labels:
        if not pos_label:
            pos_label = pos_copy
        nx.draw_networkx_labels(
            net_copy,
            pos_label,
            labels=node_labels,
            font_size=labels_size,
            font_weight="bold",
        )
    if draw_edges:
        if edge_thickness:
            nx.draw_networkx_edges(
                net_copy, pos_copy, width=[w / 2 for w in width.values()], edge_color=edge_color_map_copy
            )
        else:
            nx.draw_networkx_edges(
                net_copy, pos_copy, width=[5 for w in width.values()], edge_color=edge_color_map_copy
            )

    
    if legend_thickness:
        for i, s in enumerate(samples):
            plt.text(maxposx + .15 , maxposy - .1 * i + .02, "$w_{ij} = %s$" %s, fontsize=50)

    plt.axis("off")
    plt.axis("equal")
    plt.tight_layout()
    fig.savefig(name + ".pdf", facecolor=facecolor, bbox_inches="tight")
    # plt.show()
    plt.close()



def link_relation(node, neighbor, criterium="sequence_distance"):

    if criterium == "sequence_distance":
        node_pos = int(node[1::])
        node_chain = node[0]
        neighbor_pos = int(neighbor[1::])
        neighbor_chain = neighbor[0]
        dist = neighbor_pos - node_pos
        abs_dist = abs(dist)
        if node_chain != neighbor_chain:
            relation = "4D"
            dist = 99999
        elif abs_dist == 1:
            relation = "1D"
        elif abs_dist <= 4:
            relation = "2D"
        else:
            relation = "3D"
        return (structure, dist)

    elif criterium == "secondary_structure_elements":
        # TODO
        return


def create_ego_network(
    net,
    central_node,
    all_labels,
    color_map,
    node_borders_color,
    node_borders_secstruct,
    n_cases,
    sort,
    pathogenic=[],
    non_pathogenic=[],
    both=[],
    threshold_edge=20,
):
    ego = nx.ego_graph(net, central_node)

    if not color_map:
        color_list = []
    else:
        color_list = [color_map[n] for n in ego.nodes]
    sizes = []
    labels = {}

    for node in ego.nodes:
        if not color_map:
            seq_pos = node[1::]
            if seq_pos in pathogenic:
                color_list.append("tomato")
            elif seq_pos in non_pathogenic:
                color_list.append("limegreen")
            elif seq_pos in both:
                color_list.append("gold")
            else:
                color_list.append("lightgrey")

        k = nx.degree(net, node)
        weight = nx.degree(net, node, weight="weight")
        if weight / k < 10:
            sizes.append(200)
        elif weight / k > 16:
            sizes.append(1000)
        else:
            sizes.append(3500)
        labels[node] = all_labels[node]
    neighbors = ego.copy()
    neighbors.remove_node(central_node)

    neighbors_structure = {}

    if node_borders_color:
        if type(node_borders_color) == dict:
            node_borders = node_borders_color
        else:
            node_borders = {n: node_borders_color for n in neighbors}

    if n_cases and len(n_cases) > 0:
        check_colors = True
        edge_colors = {}
    else:
        check_colors = False
        edge_colors = "grey"

    if not node_borders_color or check_colors:

        for neighbor in neighbors:

            if not node_borders_color and node_borders_secstruct:
                struct, dist = link_relation(central_node, neighbor)  # to modify
                neighbors_structure[neighbor] = (struct, dist)
                if struct == "1D":
                    node_borders[neighbor] = "turquoise"
                elif struct == "2D":
                    node_borders[neighbor] = "turquoise"
                elif struct == "3D":
                    node_borders[neighbor] = "cadetblue"
                else:
                    node_borders[neighbor] = "cadetblue"

            else:
                node_borders[neighbor] = "k"

            if check_colors:
                num = n_cases[neighbor]
                if num < threshold_edge:
                    edge_colors[(central_node, neighbor)] = "red"
                else:
                    edge_colors[(central_node, neighbor)] = "green"

    if sort == True:
        sorted_neighbors = deque(
            sorted(neighbors_structure, key=lambda x: neighbors_structure[x][1])
        )
        initial = central_node[0] + str(int(central_node[1::]) + 1)
        if initial in sorted_neighbors:
            sorted_neighbors.rotate(-sorted_neighbors.index(initial))
        else:
            initial = central_node[0] + str(int(central_node[1::]) - 1)
            sorted_neighbors.rotate(-sorted_neighbors.index(initial) + 1)

    else:
        sorted_neighbors = list(neighbors)

    pos_original = nx.circular_layout(neighbors)
    pos = {}

    if len(sorted_neighbors) == 1:
        pos[sorted_neighbors[0]] = np.array([np.sqrt(2) / 2, np.sqrt(2) / 2])
    else:
        for i in range(len(sorted_neighbors)):
            pos[sorted_neighbors[i]] = list(pos_original.values())[i]
    pos[central_node] = np.array([0, 0])
    node_borders[central_node] = "k"
    width = nx.get_edge_attributes(ego, "weight")

    if check_colors:
        edges = ego.edges()
        edge_colors_list = []
        for u, v in edges:
            try:
                edge_colors_list.append(edge_colors[(u, v)])
            except KeyError:
                try:
                    edge_colors_list.append(edge_colors[(v, u)])
                except KeyError:
                    edge_colors_list.append("grey")
        edge_colors = edge_colors_list

    return ego, labels, pos, sizes, width, color_list, node_borders, edge_colors


def draw_neighborhood(
    net,
    pos,
    all_labels,
    pathogenic=[],
    non_pathogenic=[],
    both=[],
    save_fig=True,
    file_name="figure",
    n_cases=None,
    threshold_edge=20,
    sort=True,
    color_map=None,
    node_borders_color=None,
    node_borders_secstruct=True,
):
    plt.figure(figsize=(6, 6))
    try:
        int(pos[0])
        central_node = "A" + pos
    except ValueError:
        central_node = pos
    if central_node in net.nodes:
        ego, labels, pos, sizes, width, color_map, node_borders, edge_colors = create_ego_network(
            net,
            central_node,
            all_labels,
            color_map,
            node_borders_color,
            node_borders_secstruct,
            n_cases,
            pathogenic,
            non_pathogenic,
            both,
            threshold_edge=threshold_edge,
        )
        nx.draw_networkx_nodes(
            ego,
            pos,
            node_size=sizes,
            node_color=color_map,
            edgecolors=[node_borders[n] for n in ego.nodes],
            linewidths=4,
        )
        nx.draw_networkx_labels(
            ego, pos, labels=labels, font_size=12, font_weight="bold"
        )
        nx.draw_networkx_edges(
            ego, pos, width=[w / 5 for w in width.values()], edge_color=edge_colors
        )
        #        nx.draw_networkx_edge_labels(ego, pos, node_size=sizes, node_color=color_map, edge_labels=width)
        plt.axis("off")
        plt.tight_layout()
        if save_fig:
            plt.savefig(file_name + ".pdf")
        plt.close()
        return ego, labels, pos, sizes, width, color_map, node_borders
    else:
        print(central_node, "not in the network")
        return None, None, None, None, None, None, None


def create_perturbation_network(net1, net2, db1, db2, threshold=4):
    # check which edges are only in net1 or only in net2 and which are in common:
    edges1 = set(net1.edges)
    edges2 = set(net2.edges)
    common_edges = edges1.intersection(edges2)
    edges_1only = edges1.difference(edges2)
    edges_2only = edges2.difference(edges1)

    net = nx.Graph()
    edge_colors = {}
    size_map = []
    node_borders = []
    color_map = []
    node_labels = {}
    added_edges = []

    for u, v in common_edges:
        wij1 = net1.get_edge_data(u, v)["weight"]
        wij2 = net2.get_edge_data(u, v)["weight"]
        deltawij = wij2 - wij1
        if deltawij > threshold:
            color = "green"
            net.add_edge(u, v, weight=abs(deltawij))
            edge_colors[(u, v)] = color
            added_edges.append((u, v))
        elif deltawij < -threshold:
            color = "red"
            net.add_edge(u, v, weight=abs(deltawij))
            edge_colors[(u, v)] = color
            added_edges.append((u, v))

    for u, v in edges_1only:
        wij = net1.get_edge_data(u, v)["weight"]
        color = "red"
        if wij > threshold:
            net.add_edge(u, v, weight=wij)
            edge_colors[(u, v)] = color
            added_edges.append((u, v))

    for u, v in edges_2only:
        wij = net2.get_edge_data(u, v)["weight"]
        color = "green"
        if wij > threshold:
            net.add_edge(u, v, weight=wij)
            edge_colors[(u, v)] = color
            added_edges.append((u, v))

    for node in net.nodes:
        k1 = net1.degree(node)
        k2 = net2.degree(node)
        w1 = net1.degree(node, weight="weight")
        w2 = net2.degree(node, weight="weight")
        if isinstance(k1, int) and isinstance(k2, int):
            nw1 = w1 / k1
            nw2 = w2 / k2

            info1 = db1[db1["Position"] == node]
            if len(info1) > 1:
                info1 = info1.iloc[0]  # check why more than one
            type_aa_1 = aac.three2one(
                info1["Type of residue"].item()
            )

            info2 = db2[db2["Position"] == node]
            if len(info2) > 1:
                info2 = info2.iloc[0]  # check why more than one
            type_aa_2 = aac.three2one(
                info2["Type of residue"].item()
            )

            if type_aa_1 == type_aa_2:
                color_map.append("gray")
                label = type_aa_1 + node[1::] + ":" + node[0]
                node_labels[node] = label

            else:
                color_map.append("blue")
                label = type_aa_1 + "-" + type_aa_2 + node[1::] + ":" + node[0]
                node_labels[node] = label

        elif isinstance(k2, int):
            nw1 = 0
            nw2 = w2 / k2
            color_map.append("green")
            info2 = db2[db2["Position"] == node]
            if len(info2) > 1:
                info2 = info2.iloc[0]  # check why more than one
            type_aa_2 = aac.three2one(
                info2["Type of residue"].item()
            )
            label = type_aa_2 + node[1::] + ":" + node[0]
            node_labels[node] = label

        else:
            nw1 = w1 / k1
            nw2 = 0
            color_map.append("red")
            info1 = db1[db1["Position"] == node]
            if len(info1) > 1:
                info1 = info1.iloc[0]  # check why more than one
            type_aa_1 = aac.three2one(
                info1["Type of residue"].item()
            )
            label = type_aa_1 + node[1::] + ":" + node[0]
            node_labels[node] = label

        deltanw = nw2 - nw1
        size_map.append(abs(deltanw) * 100 + 500)
        if deltanw > 0:
            node_borders.append("green")
        elif deltanw < 0:
            node_borders.append("red")
        else:
            node_borders.append("gray")
    #
    #    node_labels = {node: node for node in net.nodes} #TO DO: labels as in normal networks, but with mutations

    #    edge_color_map = [edge_colors[edge] for edge in added_edges]

    return (
        net,
        node_labels,
        size_map,
        color_map,
        edge_colors,
        node_borders,
        edges_1only,
        edges_2only,
        common_edges,
    )


def create_ego_perturbation_network(
    net, central_node, all_labels, size_map, color_map, node_borders, edges_color_map
):
    ego = nx.ego_graph(net, central_node)
    labels = {node: all_labels[node] for node in ego.nodes}
    neighbors = ego.copy()
    neighbors.remove_node(central_node)

    node_borders_dict = {node: node_borders[i] for i, node in enumerate(net.nodes)}
    color_map_dict = {node: color_map[i] for i, node in enumerate(net.nodes)}
    size_map_dict = {node: size_map[i] for i, node in enumerate(net.nodes)}

    node_borders_ego = {node: node_borders_dict[node] for node in ego.nodes}
    color_map_ego = [color_map_dict[node] for node in ego.nodes]
    size_map_ego = [size_map_dict[node] for node in ego.nodes]

    edges_color_map_dict = {
        edge: edges_color_map[i] for i, edge in enumerate(net.edges)
    }
    edges_color_map_ego = []
    for u, v in ego.edges:
        try:
            edges_color_map_ego.append(edges_color_map_dict[(u, v)])
        except KeyError:
            edges_color_map_ego.append(edges_color_map_dict[(v, u)])

    if not node_borders:
        neighbors_structure = {}
        node_borders = {}
        node_borders[central_node] = "grey"
        for neighbor in neighbors:
            struct, dist = link_relation(central_node, neighbor)  # to modify
            neighbors_structure[neighbor] = (struct, dist)
            if struct == "1D":
                node_borders[neighbor] = "turquoise"
            elif struct == "2D":
                node_borders[neighbor] = "turquoise"
            elif struct == "3D":
                node_borders[neighbor] = "cadetblue"
            else:
                node_borders[neighbor] = "cadetblue"

    pos_original = nx.circular_layout(neighbors)
    pos = {}
    for i in range(len(neighbors.nodes)):
        pos[list(neighbors.nodes)[i]] = list(pos_original.values())[i]
    pos[central_node] = np.array([0, 0])
    width = nx.get_edge_attributes(ego, "weight")

    return (
        ego,
        labels,
        pos,
        size_map_ego,
        width,
        color_map_ego,
        node_borders_ego,
        edges_color_map_ego,
    )


def draw_neighborhood_perturbation(
    net,
    pos,
    all_labels,
    size_map,
    color_map,
    edge_color_map,
    node_borders,
    save_fig=True,
    file_name="figure",
):
    plt.figure(figsize=(6, 6))
    try:
        int(pos[0])
        central_node = "A" + pos
    except ValueError:
        central_node = pos
    if central_node in net.nodes:
        ego, labels, pos, sizes, width, color_map, node_borders, edges_color_map_ego = create_ego_perturbation_network(
            net,
            central_node,
            all_labels,
            size_map,
            color_map,
            node_borders,
            edge_color_map,
        )
        nx.draw_networkx_nodes(
            ego,
            pos,
            node_size=sizes,
            node_color=color_map,
            edgecolors=[node_borders[n] for n in ego.nodes],
            linewidths=4,
        )
        nx.draw_networkx_labels(
            ego, pos, labels=labels, font_size=12, font_weight="bold"
        )
        nx.draw_networkx_edges(
            ego,
            pos,
            width=[w / 5 for w in width.values()],
            edge_color=edges_color_map_ego,
        )

        plt.axis("off")
        plt.tight_layout()
        if save_fig:
            plt.savefig(file_name + ".pdf")
        plt.close()
        return ego, labels, pos, sizes, width, color_map, node_borders, edge_color_map
    else:
        print(central_node, "not in the network")
        return None, None, None, None, None, None, None


def plot_legends(legends=["nw_color", "link_color"], folder_path=""):

    legends = set(legends)

    G1 = nx.Graph()
    colors1 = ["blue", "cyan", "greenyellow", "yellow", "orange", "red"]
    for n in range(1, 7):
        G1.add_node(str(n))

    pos1 = {}
    for i, node in enumerate(G1.nodes):
        pos1[node] = np.array([0, i + 0.15])

    G2 = nx.Graph()
    sizes2 = [1, 2, 3, 4]
    sizes2 = [s * 400 for s in sizes2]
    for n in range(1, 5):
        G2.add_node(str(n))

    pos2 = {}
    for i, node in enumerate(G2.nodes):
        pos2[node] = np.array([0, 3 * i + 0.15])

    G2b = nx.Graph()
    sizes2b = [1, 2, 3, 4, 5, 6]
    sizes2b = [s * 400 for s in sizes2]
    for n in range(1, 7):
        G2b.add_node(str(n))

    pos2b = {}
    for i, node in enumerate(G2b.nodes):
        pos2b[node] = np.array([0, 3 * i + 0.15])

    if "nw_color" in legends:

        plt.figure()
        nx.draw_networkx_nodes(G1, pos1, node_color=colors1, edgecolors="k")
        plt.axis("off")
        plt.xlim(-0.2, 3)
        plt.ylim(-1, 6)
        plt.rc("text", usetex=True)
        plt.rc("font", family="calibri")
        plt.text(0.5, 0, "W/k $<$ 5", fontsize=16)
        plt.text(0.5, 1, "5 $\leq$ W/k $<$ 10", fontsize=16)
        plt.text(0.5, 2, "10 $\leq$ W/k $<$ 15", fontsize=16)
        plt.text(0.5, 3, "15 $\leq$ W/k $<$ 20", fontsize=16)
        plt.text(0.5, 4, "20 $\leq$ W/k $<$ 25", fontsize=16)
        plt.text(0.5, 5, "W/k $\geq$ 25", fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(folder_path, "network_pictures", "legend_color.pdf"))
        plt.close()

    if "link_color" in legends:

        texts = [
            "$w_{ij} \  <$ 10",
            "10 $\leq \ w_{ij} \  <$ 20",
            "20 $\leq \ w_{ij} \  <$ 30",
            "30 $\leq \ w_{ij} \  <$ 40",
            "40 $\leq \ w_{ij} \  <$ 50",
            "$w_{ij} \  \geq$ 50",
        ]
        fig = plt.figure()
        plt.rc("text", usetex=True)
        plt.rc("font", family="calibri")
        for index, c in enumerate(colors1):
            subf = int("61%s" % (index + 1))
            ax1 = fig.add_subplot(subf)
            ax1.add_patch(
                patches.Rectangle(
                    (0.1, 0.1), 0.1, 0.5, facecolor=c, edgecolor="k", linewidth=0.1
                )
            )
            ax1.axis("off")
            ax1.text(0.3, 0, texts[index], fontsize=16)
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(folder_path, "network_pictures", "legend_link_color.pdf"))
        plt.close()

    if "k_size" in legends:  # TO REMOVE?

        plt.figure()
        nx.draw_networkx_nodes(
            G2, pos2, node_size=sizes2, node_color="w", edgecolors="k"
        )
        plt.axis("off")
        plt.xlim(-0.3, 3)
        plt.ylim(-1, 6)
        plt.rc("text", usetex=True)
        plt.rc("font", family="calibri")
        plt.text(0.5, 0, "k $<$ 5", fontsize=16)
        plt.text(0.5, 1, "5 $\leq$ k $<$ 9", fontsize=16)
        plt.text(0.5, 2, "9 $\leq$ k $<$ 13", fontsize=16)
        plt.text(0.5, 3, "k $\geq$ 13", fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(folder_path, "network_pictures", "legend_size_k.pdf"))
        plt.close()

    if "nw_size" in legends:  # TO REMOVE?
        plt.figure()
        nx.draw_networkx_nodes(
            G2b, pos2b, node_size=sizes2b, node_color="w", edgecolors="k"
        )
        plt.axis("off")
        plt.xlim(-0.3, 3)
        plt.ylim(-1, 17)
        plt.rc("text", usetex=True)
        plt.rc("font", family="calibri")
        plt.text(0.5, 0, "W/k $<$ 5", fontsize=16)
        plt.text(0.5, 3, "5 $\leq$ W/k $<$ 10", fontsize=16)
        plt.text(0.5, 6, "10 $\leq$ W/k $<$ 15", fontsize=16)
        plt.text(0.5, 9, "15 $\leq$ W/k $<$ 20", fontsize=16)
        plt.text(0.5, 12, "20 $\leq$ W/k $<$ 25", fontsize=16)
        plt.text(0.5, 15, "W/k $\geq$ 25", fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(folder_path, "network_pictures", "legend_size_nw.pdf"))
        plt.close()


def subnets_chains(net, chains):
    subnets = {chain: nx.Graph() for chain in chains}
    for n in net.nodes:
        ch = n[0]
        subnets[ch].add_node(n)

    for u, v in net.edges:
        chu = u[0]
        chv = v[0]

        if chu == chv:
            w = net.get_edge_data(u, v)["weight"]
            subnets[chu].add_edge(u, v, weight=w)
    return subnets


def nodes_position_perturbation(net, dim, chains, max_x=20, scale=3, columns=False):
    if dim == "" or dim == "4D" or dim == "3-4D":
        nodes_pos = {}
        num_chains = len(chains)

        if num_chains > 2 or columns:
            x_interval = 30 / num_chains
            x_values = [x * x_interval for x in range(num_chains)]
            y_values = [0 for x in range(num_chains)]
            x_shifts = [1 for x in range(num_chains)]

            for i, node in enumerate(sorted(net.nodes, key=lambda k: int(k[1::]))):
                chain = node[0]
                c = chains.index(chain)
                x = x_values[c]
                y = y_values[c]
                x_shift = x_shifts[c]
                nodes_pos[node] = (x + x_shift, y)
                y_values[c] += 1
                x_shifts[c] *= -1

        else:
            subnets = subnets_chains(net, chains)

            x_interval = max_x / num_chains
            x_values = [x * x_interval for x in range(num_chains)]

            for i, chain in enumerate(chains):
                subnet = subnets[chain]
                pos = nx.spring_layout(
                    subnet, k=7 / np.sqrt(len(subnet.nodes)), scale=scale
                )
                xvalue = x_values[i]
                for n in pos:
                    pos[n][0] += xvalue
                nodes_pos.update(pos)

    else:
        nodes_pos = nx.spring_layout(net, k=2)

    return nodes_pos


def induced_perturbation_network(net_p, source):

    try:
        assert (
            source in net_p.nodes
        ), "No induced perturbation network from source %s" % (source)
        tree = nx.bfs_tree(net_p, source)
    except AssertionError as ex:
        print(ex)
        return

    # add triangles
    for v in tree.nodes:
        for n in net_p.neighbors(v):
            if n in tree.nodes:
                tree.add_edge(v, n)

    tree = tree.to_undirected()

    return tree
