import copy
import time

from gurobipy import *

class Graph:
    def __init__(self):
        self.node_list = []
        self.edge_list = []
        self.node_dict = {}
        self.node_number = 0
        self.max_density = -1
        self.densest_set = []

    def add_node(self, node_id):
        if node_id in self.node_dict:
            # print("Warning in adding node! Existing node")
            return
        else:
            self.node_number += 1
            temp_node = Node(node_id, len(self.node_list))
            self.node_list.append(temp_node)
            self.node_dict[node_id] = temp_node.index

    def add_edge(self, node_pair, length):
        # this method contains add_node
        for i in node_pair:
            if i not in self.node_dict:
                self.add_node(i)
        temp_edge = Edge(node_pair[0], node_pair[1], length)
        self.edge_list.append(temp_edge)
        self.FNBI(node_pair[0]).add_neighbor(node_pair[1])
        self.FNBI(node_pair[1]).add_neighbor(node_pair[0])

    def FNBI(self, node_id):    # Find node by node.id
        return self.node_list[self.node_dict[node_id]]

    def remove_node(self, node_id):
        self.densest_set = []
        self.max_density = -1
        for i in self.FNBI(node_id).neighbor:
            self.FNBI(i).remove_neighbor(node_id)
        for n, i in enumerate(self.edge_list):
            if i == None:
                continue
            if i.n1 == node_id or i.n2 == node_id:
                self.edge_list[n] = None
        self.node_list[self.node_dict[node_id]] = None
        self.node_number -= 1

    def update(self):
        # remove None data and update node_dict
        new_node_list = []
        new_node_dict = {}
        new_edge_list = []
        for i in self.node_list:
            if i == None:
                continue
            new_node_dict[i.id] = len(new_node_list)
            new_node_list.append(i)
        for i in self.edge_list:
            if i == None:
                continue
            new_edge_list.append(i)
        self.node_list = new_node_list
        self.node_dict = new_node_dict
        self.edge_list = new_edge_list

        self.max_density = -1
        self.densest_set = []

    def subgraph(self, node_ids):
        if node_ids == None:
            return None
        sg = Graph()
        for i in node_ids:
            sg.add_node(i)
        for i in self.edge_list:
            if i.n1 in node_ids and i.n2 in node_ids:
                sg.add_edge((i.n1, i.n2), i.length)
        return sg


class Node:
    def __init__(self, id, index=None):
        self.id = id
        self.index = index
        self.neighbor = []
        self.degree = 0
        self.weight = 0

    def __repr__(self):
        return "node " + self.id

    def add_neighbor(self, id):
        self.neighbor.append(id)
        self.degree += 1

    def remove_neighbor(self, id):
        try:
            self.neighbor.remove(id)
            self.degree -= 1
        except:
            print("Error! Node {} has no neighbor {}!".format(self.id, id))

    def copy_node(self):
        return copy.deepcopy(self)


class Edge:
    def __init__(self, n1, n2, length):
        self.n1 = n1
        self.n2 = n2
        self.weight = 0
        self.length = length

    def __repr__(self):
        return "edge {}-{}".format(self.n1, self.n2)


ori_g = Graph()

with open("../data/SKEMPI/graph_hotspot_o_ring/5C6T_o.txt") as f1:
# with open("data/test3.txt") as f1:
    for n, i in enumerate(f1):
        if n == 0:
            continue
        i = i.split()
        ori_g.add_edge((i[0], i[2]), i[4])


def LP(g, log_flag=0, var_details=0):
    n_list = g.node_list
    e_list = g.edge_list
    if len(e_list) == 0 or len(n_list) == 0:
        g.densest_set = []
        g.max_density = 0
        return
    model_1 = Model("test1")
    model_1.setParam("LogToConsole", log_flag)

    variables_x_dict = {}
    variables_y_dict = {}
    for i in e_list:
        variables_x_dict["{}-{}".format(i.n1, i.n2)] = model_1.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1,
                                                                      name="x-{}-{}".format(i.n1, i.n2))

    for i in n_list:
        variables_y_dict["{}".format(i.id)] = model_1.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, name="y-{}".format(i.id))

    variables_x_list = list(variables_x_dict[i] for i in variables_x_dict)
    variables_y_list = list(variables_y_dict[i] for i in variables_y_dict)
    model_1.setObjective(quicksum(variables_x_list), GRB.MAXIMIZE)

    model_1.addConstr(quicksum(variables_y_list) <= 1)
    for n, i in enumerate(e_list):
        model_1.addConstr(variables_x_dict["{}-{}".format(i.n1, i.n2)] <= variables_y_dict["{}".format(i.n1)],
                          "c-{}-{}".format(n, 0))
        model_1.addConstr(variables_x_dict["{}-{}".format(i.n1, i.n2)] <= variables_y_dict["{}".format(i.n2)],
                          "c-{}-{}".format(n, 1))
    model_1.optimize()

    g.max_density = int(model_1.objVal * 1000000) / 1000000     # floor
    for v in model_1.getVars():
        if v.varName[0] == "y":
            if log_flag and var_details:
                print("{}, {}".format(v.varName, v.x))
            if v.x >= 0.5 / len(n_list):
                g.densest_set.append(v.varName.split("-")[1])
    if log_flag and var_details:
        for v in model_1.getVars():
            if v.varName[0] == "x":
                print("{}, {}".format(v.varName, v.x))
                pass

def LP4(g, max_density, R):
    n_list = g.node_list
    e_list = g.edge_list
    epsilon = 1 / len(n_list) / 2
    max_density = int(max_density * 1000000) / 1000000
    model_1 = Model("test1")
    model_1.setParam("LogToConsole", 0)

    variables_x_dict = {}
    variables_y_dict = {}
    for i in e_list:
        variables_x_dict["{}-{}".format(i.n1, i.n2)] = model_1.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1,
                                                                      name="x-{}-{}".format(i.n1, i.n2))
    for i in n_list:
        variables_y_dict["{}".format(i.id)] = model_1.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, name="y-{}".format(i.id))

    variables_x_list = list(variables_x_dict[i] for i in variables_x_dict)
    variables_y_list = list(variables_y_dict[i] for i in variables_y_dict)


    model_1.addConstr(quicksum(variables_y_list) <= 1)
    model_1.addConstr(quicksum(variables_x_list) >= max_density)
    for i in R:
        model_1.addConstr(variables_y_dict[i] >= epsilon)

    y_list_2 = list(variables_y_dict["{}".format(i.id)] for i in n_list if i.id not in R)
    model_1.addConstr(quicksum(y_list_2) >= epsilon)

    model_1.setObjective(quicksum(variables_x_list), GRB.MAXIMIZE)

    for n, i in enumerate(e_list):
        model_1.addConstr(variables_x_dict["{}-{}".format(i.n1, i.n2)] <= variables_y_dict["{}".format(i.n1)],
                          "c-{}-{}".format(n, 0))
        model_1.addConstr(variables_x_dict["{}-{}".format(i.n1, i.n2)] <= variables_y_dict["{}".format(i.n2)],
                          "c-{}-{}".format(n, 1))
    model_1.optimize()

    if model_1.status != 2:
        return None
    else:
        g.max_density = max_density
        g.densest_set = []
        # print()
        for v in model_1.getVars():
            if v.varName[0] == "y":
                # print("{}, {}".format(v.varName, v.x))
                if v.x > 0:
                    g.densest_set.append(v.varName.split("-")[1])
        return g.densest_set


def preprocess(g):
    if g.max_density == -1:
        LP(g)
    max_density = g.max_density
    g_ = copy.deepcopy(g)
    while g_.node_number > 0:
        stop_flag = True
        for i in g_.node_list:
            if i == None:
                continue
            if len(i.neighbor) < max_density:
                g_.remove_node(i.id)
                stop_flag = False
        if stop_flag:
            break
    g_.update()
    LP(g_)
    return g_

def hotspot_read(path):
    g = Graph()
    with open(path) as f1:
        for n, i in enumerate(f1):
            if n == 0:
                continue
            i = i.split()
            g.add_edge((i[0], i[2]), i[4])
    return g

def find_maximal(PDB_id):
    g = hotspot_read("../data/SKEMPI/graph_hotspot_o_ring/{}.txt".format(PDB_id))
    g_ = preprocess(g)
    result = []
    max_density = g_.max_density
    while True:
        temp_result = LP4(g_, max_density, result)
        if temp_result == None:
            break
        result = temp_result

    with open("result2/{}_maximal_ver2.txt".format(PDB_id), "w") as f1:
        for i in result:
            f1.write("{}\n".format(i))


start_time = time.time()

graph_files = os.walk("../data/SKEMPI/graph_hotspot_o_ring")
for path, dir_list, file_list in graph_files:
    count = 1
    for file_name in file_list:
        temp_pdb_id = file_name.split(".")[0]
        print("Processing\t{}/{}:\t{}".format(count, len(file_list), temp_pdb_id))
        find_maximal(temp_pdb_id)
        count += 1

end_time = time.time()

print("Time used:", end_time - start_time, "s")
print("End of program")
# Note: the original LP always return an MDS