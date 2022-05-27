import copy
import math

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

    def subgraph(self, node_ids, max_density=-1):
        if node_ids == None:
            return None
        sg = Graph()
        for i in node_ids:
            sg.add_node(i)
        for i in self.edge_list:
            if i.n1 in node_ids and i.n2 in node_ids:
                sg.add_edge((i.n1, i.n2), i.length)
        sg.max_density = max_density
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

    g.max_density = floor(model_1.objVal)
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

def LP2(g, max_density, u):
    n_list = g.node_list
    e_list = g.edge_list
    max_density = floor(max_density)
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
    model_1.setObjective(variables_y_dict[u], GRB.MAXIMIZE)

    model_1.addConstr(quicksum(variables_y_list) <= 1)
    model_1.addConstr(quicksum(variables_x_list) >= max_density)
    model_1.addConstr(variables_y_dict[u] >= 1.0 / len(n_list) / 2)
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
        for v in model_1.getVars():
            if v.varName[0] == "y":
                # print("{}, {}".format(v.varName, v.x))
                lower_bound = floor(1.0 / len(n_list)) / 2
                if v.x >= lower_bound:
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
    if g_.node_number == 0:
        return None
    g_.update()
    LP(g_)
    return g_

def find_minimal(g):
    g_ = copy.deepcopy(g)
    g_.update()
    LP(g_)
    if g_.densest_set == []:
        return None
    # g_ = preprocess(g)
    # LP(g_)

    h = g_
    max_density = g_.max_density
    while True:
        u = h.densest_set[0]
        # print(" ".join(h.densest_set))
        h1 = try_remove(u, h)
        h2 = try_enhance(u, h)
        if h1 == None:
            # print(h2)
            return h.subgraph(h2, max_density=max_density)
        h_ = None
        if len(h1) > len(h2):
            h_ = h2
        else:
            h_ = h1
        h = h.subgraph(h_, max_density=max_density)
        LP(h)
        # print(h.densest_set)
        # h.densest_set = copy.deepcopy(h_)

def hotspot_read(path):
    g = Graph()
    with open(path) as f1:
        for n, i in enumerate(f1):
            if n == 0:
                continue
            i = i.split()
            g.add_edge((i[0], i[2]), i[4])
    return g

def find_all_minimal(PDB_id):
    tolerance = 0.91
    g = hotspot_read("../data/SKEMPI/graph_hotspot_o_ring/{}.txt".format(PDB_id))
    g_ = copy.deepcopy(g)
    g_.update()
    LP(g_)
    # g_ = preprocess(g)
    # LP(g_)
    result = []
    max_density = g_.max_density
    while True:
        temp_result = find_minimal(g_)
        if temp_result == None:
            break
        if temp_result.max_density < max_density * tolerance:
            break
        result = result + temp_result.node_list
        id_list1 = list(i.id for i in g_.node_list)
        id_list2 = list(i.id for i in temp_result.node_list)
        remained_nodes = list(i for i in id_list1 if i not in id_list2)
        if remained_nodes == None:
            break
        g_ = g_.subgraph(remained_nodes)
        g_.update()
        LP(g_)
        # g_ = preprocess(g_)


    with open("result/{}_minimal_sub.txt".format(PDB_id), "w") as f1:
        for i in result:
            f1.write("{}\n".format(i.id))

    # return result


def try_remove(u, g):
    g_ = copy.deepcopy(g)
    g_.remove_node(u)
    g_.update()
    LP(g_)
    # print(g_.densest_set)
    if approximate_equal(g_.max_density, g.max_density):
        return g_.densest_set
    else:
        return None


def try_enhance(u, g):
    g_ = copy.deepcopy(g)
    g_.max_density = -1
    g_.densest_set = []
    return LP2(g_, g.max_density, u)


def approximate_equal(a, b, precision=0.00001):
    if abs(a-b) < precision:
        return True
    else:
        return False


def floor(a, precision=0.0000000001):
    return int(a/precision) * precision


graph_files = os.walk("../data/SKEMPI/graph_hotspot_o_ring")
for path, dir_list, file_list in graph_files:
    count = 1
    for file_name in file_list:
        temp_pdb_id = file_name.split(".")[0]
        print("Processing\t{}/{}:\t{}".format(count, len(file_list), temp_pdb_id))
        find_all_minimal(temp_pdb_id)
        count += 1

print("End of program")

# Note: the original LP always return an MDS