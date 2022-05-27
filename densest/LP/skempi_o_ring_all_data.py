from gurobipy import *
import os

class Node:
    def __init__(self, id, index):
        self.id = id
        self.index = index
        self.neighbor = []
        self.degree = 0
        self.weight = 0


class Edge:
    def __init__(self, n1, n2, length):
        self.n1 = n1
        self.n2 = n2
        self.weight = 0
        self.length = length

too_small_y_list = []

def skempi_gurobi(PDB_id):
    node_list = []
    node_id_dict = {}
    edge_list = []

    with open("../../data/SKEMPI/graph_hotspot_o_ring/{}.txt".format(PDB_id)) as f1:
        for n, i in enumerate(f1):
            if n == 0:
                continue
            i = i.split()
            if i[0] not in node_id_dict:
                temp_node = Node(i[0], len(node_list))

                node_id_dict[i[0]] = len(node_list)
                node_list.append(temp_node)

            if i[2] not in node_id_dict:
                temp_node = Node(i[2], len(node_list))
                node_id_dict[i[2]] = len(node_list)
                node_list.append(temp_node)

            node_list[node_id_dict[i[0]]].neighbor.append(i[2])
            node_list[node_id_dict[i[2]]].neighbor.append(i[0])
            temp_edge = Edge(i[0], i[2], i[-1])
            edge_list.append(temp_edge)

    model_1 = Model("test1")
    model_1.setParam("LogToConsole", 0)
    variables_x_dict = {}
    variables_y_dict = {}
    for i in edge_list:
        variables_x_dict["{}-{}".format(i.n1, i.n2)] = model_1.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, name="x-{}-{}".format(i.n1, i.n2))

    for i in node_list:
        variables_y_dict["{}".format(i.id)] = model_1.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, name="y-{}".format(i.id))

    variables_x_list = list(variables_x_dict[i] for i in variables_x_dict)
    variables_y_list = list(variables_y_dict[i] for i in variables_y_dict)
    model_1.setObjective(quicksum(variables_x_list), GRB.MAXIMIZE)

    model_1.addConstr(quicksum(variables_y_list) <= 1)
    for n, i in enumerate(edge_list):
        model_1.addConstr(variables_x_dict["{}-{}".format(i.n1, i.n2)] <= variables_y_dict["{}".format(i.n1)], "c-{}-{}".format(n, 0))
        model_1.addConstr(variables_x_dict["{}-{}".format(i.n1, i.n2)] <= variables_y_dict["{}".format(i.n2)], "c-{}-{}".format(n, 1))
    model_1.optimize()

    results = []
    with open("result_o_ring_log/{}_log.txt".format(PDB_id), "w") as f3:
        f3.write("n = {}\n".format(len(node_list)))
        f3.write("1/n = {:.10}\n\n".format(1.0/len(node_list)))

        f3.write("########## y section ##########\n")
        temp_y_list = []
        for v in model_1.getVars():
            if v.varName[0] == "y":
                temp_y_list.append([v.varName, v.x])
                if v.x > 0:
                    results.append(v.varName.split("-")[1])
        temp_y_list.sort(reverse=True, key=lambda x: x[1])
        y_positive_counter = 0
        y_section_counter = 1
        y_section = round(temp_y_list[0][1], 10)
        y_counter = 0
        for i in temp_y_list:
            if i[1] > 0.0:
                y_positive_counter += 1
            if round(i[1], 10) == y_section:
                y_counter += 1
            else:
                f3.write("y section {}: {}\n".format(y_section_counter, y_section))
                f3.write("Number of y: {}\n\n".format(y_counter))
                y_counter = 1
                y_section_counter += 1
                y_section = round(i[1], 10)
        f3.write("y section {}: {}\n".format(y_section_counter, y_section))
        f3.write("Number of y: {}\n\n".format(y_counter))


        if y_positive_counter != len(results):
            f3.write("Warning: found too small y value!\n")
            too_small_y_list.append(PDB_id)
        f3.write("Number of y > 0.0: {}\n".format(y_positive_counter))
        f3.write("Number of y > 1/n: {}\n\n".format(len(results)))

        f3.write("########## x section ##########\n")
        temp_x_list = []
        for v in model_1.getVars():
            if v.varName[0] == "x":
                temp_x_list.append([v.varName, v.x])
        temp_x_list.sort(reverse=True, key=lambda x: x[1])
        x_positive_counter = 0
        x_section_counter = 1
        x_section = round(temp_x_list[0][1], 10)
        x_counter = 0
        for i in temp_x_list:
            if i[1] > 0.0:
                x_positive_counter += 1
            if round(i[1], 10) == x_section:
                x_counter += 1
            else:
                if x_section == 0.0:
                    break
                f3.write("x section {}: {}\n".format(x_section_counter, x_section))
                f3.write("Number of x: {}\n\n".format(x_counter))
                x_counter = 1
                x_section_counter += 1
                x_section = round(i[1], 10)
        f3.write("x section {}: {}\n".format(x_section_counter, x_section))
        f3.write("Number of x: {}\n\n".format(x_counter))
        f3.write("Number of x > 0.0: {}\n\n".format(x_positive_counter))

        f3.write("########## details ##########\n")
        for i in temp_y_list:
            f3.write("{}\t{}\n".format(i[0], i[1]))
        for i in temp_x_list:
            f3.write("{}\t{}\n".format(i[0], i[1]))

    print("Number of selected nodes:", len(results))
    print(" ".join(results))
    with open("result_o_ring/{}_dense.txt".format(PDB_id), "w") as f2:
        for i in results:
            f2.write("{}\n".format(i))


# skempi_gurobi("1AO7_o")
# skempi_gurobi("1A22_o")

graph_files = os.walk("../../data/SKEMPI/graph_hotspot_o_ring")
for path, dir_list, file_list in graph_files:
    for file_name in file_list:
        temp_pdb_id = file_name.split(".")[0]
        skempi_gurobi(temp_pdb_id)


if len(too_small_y_list) > 0:
    with open("too_small_y_list.txt", "w") as f4:
        for i in too_small_y_list:
            f4.write("{}\n".format(i))
