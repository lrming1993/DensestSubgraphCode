import os


class Node:
    node_list = []
    node_id_dict = {}
    def __init__(self, id, index):
        self.id = id
        self.index = index
        self.neighbor = []
        self.degree = 0
        self.weight = 0

    @staticmethod
    def common_neighbor(node_ids):
        output = set(Node.find_node_by_id(node_ids[0]).neighbor)
        for i in node_ids:
            output.intersection_update(set(Node.find_node_by_id(i).neighbor))
        return list(output)

    @staticmethod
    def closure(node_ids):
        return Node.common_neighbor(Node.common_neighbor(node_ids))

    @staticmethod
    def support(node_ids):
        return len(Node.common_neighbor(node_ids))

    @staticmethod
    def tail(node_ids):
        max_index = -1
        for i in node_ids:
            temp_node = Node.find_node_by_id(i)
            if max_index < temp_node.index:
                max_index = temp_node.index
        if max_index >= len(Node.node_list) - 1:
            return []
        return [i.id for i in Node.node_list[max_index + 1:]]

    @staticmethod
    def is_closed(node_ids):
        closure = Node.closure(node_ids)
        if set(node_ids) == set(closure):
            return True
        return False

    @staticmethod
    def find_node_by_id(id):
        return Node.node_list[Node.node_id_dict[id]]

    @staticmethod
    def closed_pattern_pair(pattern):
        return Node.common_neighbor(pattern)

    @staticmethod
    def clear():
        Node.node_id_dict = {}
        Node.node_list = []

class Edge:
    edge_list = []
    def __init__(self, n1, n2):
        self.n1 = n1
        self.n2 = n2
        self.weight = 0

    @staticmethod
    def clear():
        Edge.edge_list = []

def hotspot_read(path):
    with open(path) as f1:
        for n, i in enumerate(f1):
            if n == 0:
                continue
            i = i.split()
            if i[0] not in Node.node_id_dict:
                temp_node = Node(i[0], len(Node.node_list))

                Node.node_id_dict[i[0]] = len(Node.node_list)
                Node.node_list.append(temp_node)

            if i[2] not in Node.node_id_dict:
                temp_node = Node(i[2], len(Node.node_list))
                Node.node_id_dict[i[2]] = len(Node.node_list)
                Node.node_list.append(temp_node)

            Node.find_node_by_id(i[0]).neighbor.append(i[2])
            Node.find_node_by_id(i[2]).neighbor.append(i[0])
            temp_edge = Edge(i[0], i[2])
            Edge.edge_list.append(temp_edge)

def LCM_MBC(X, tail_X, p, q):
    output = []
    new_tail_X = tail_X.copy()
    for v in tail_X:
        set_X = set(X)
        set_X.add(v)
        if Node.support(list(set_X)) < q:
            new_tail_X.remove(v)
    if len(X) + len(new_tail_X) < p:
        return
    for v in new_tail_X:
        node_v = Node.find_node_by_id(v)
        set_X = set(X)
        set_X.add(v)
        Y = list(set_X)
        tail_Y = [i for i in new_tail_X if Node.find_node_by_id(i).index > node_v.index]
        sup_Y = Node.support(Y)
        set_Y = set(Y)
        Z = []
        if len(Y) + len(tail_Y) >= p:
            Z = Y.copy()
            for u in tail_Y:
                new_set_Y = set_Y.copy()
                new_set_Y.add(u)
                if Node.support(list(new_set_Y)) == sup_Y:
                    Z.append(u)
            if Node.is_closed(Z):
                sup_Z = Node.support(Z)
                if len(Z) >= p and sup_Z >= len(Z):
                    # print(Z, Node.closed_pattern_pair(Z))
                    output.append([Z, Node.closed_pattern_pair(Z)])
                if sup_Z > len(Z):
                    output += LCM_MBC(Z, set(tail_Y) - set(Z), p, q)
    return output

def duplicate_removal(li):
    output = []
    check_list = []
    for i in li:
        if len(i) != 2:
            print("Error from duplicate_removal:")
            print("Invalid input list:", i)
        if len(i[0]) != len(i[1]):
            output.append(i)
        else:
            temp_list = i[0] + i[1]
            temp_list.sort()
            if temp_list in check_list:
                continue
            else:
                output.append(i)
                check_list.append(temp_list)
    return output

def clear_all():
    Node.clear()
    Edge.clear()

def sigle_read(PDB_id):
    clear_all()
    print("Processing:", PDB_id)
    hotspot_read("../data/SKEMPI/graph_hotspot_o_ring/{}.txt".format(PDB_id))
    test_output = LCM_MBC([], Node.tail([]), 1, 2)
    # print(test_output)
    test_output2 = duplicate_removal(test_output)
    with open("log/{}_log.txt".format(PDB_id), "w") as f2:
        for i in test_output2:
            f2.write("{} | {}\n".format(" ".join(i[0]), " ".join(i[1])))
    # print(test_output2)
    with open("result/{}_LCM_MBC.txt".format(PDB_id), "w") as f3:
        for i in test_output2:
            if len(i[0]) < 3:
                continue
            f3.write("{}\n".format(" ".join(i[0] + i[1])))

# sigle_read("1AO7_o")

graph_files = os.walk("../data/SKEMPI/graph_hotspot_o_ring")
for path, dir_list, file_list in graph_files:
    for file_name in file_list:
        temp_pdb_id = file_name.split(".")[0]
        sigle_read(temp_pdb_id)

print("\nEnd of process")
