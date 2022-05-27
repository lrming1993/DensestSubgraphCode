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
    def add_node_list(node):
        Node.node_id_dict[node.id] = len(Node.node_list)
        Node.node_list.append(node)

    @staticmethod
    def common_neighbor(node_ids):
        output = set(Node.find_node_by_id(node_ids[0]).neighbor)
        for i in node_ids:
            output.intersection_update(set(Node.find_node_by_id(i).neighbor))
        return list(output)

    @staticmethod
    def find_node_by_id(id):
        return Node.node_list[Node.node_id_dict[id]]

    @staticmethod
    def clear():
        Node.node_id_dict = {}
        Node.node_list = []


class Edge:
    edge_list = []

    def __init__(self, n1, n2, w):
        self.n1 = n1
        self.n2 = n2
        self.w = w
        self.id = "{}-{}".format(n1, n2)

    def __repr__(self):
        return "Edge %s :%s" % (self.id, self.w)

    @staticmethod
    def clear():
        Edge.edge_list = []


class Pipe:
    pipe_list = []
    pipe_dict = {}

    def __init__(self, source, sink, capacity):
        self.source = source
        self.sink = sink
        self.capacity = capacity
        self.residual = capacity
        self.rpipe = None
        self.id = "{}-{}".format(source, sink)
        self.flow = 0

    def __repr__(self):
        return "%s -> %s: %s / %s (%s)" % (self.source, self.sink, self.flow, self.capacity, self.residual)

    def add_flow(self, f):
        self.flow += f
        self.residual -= f
        self.rpipe.flow -= f
        self.rpipe.residual += f

    @staticmethod
    def add_pipe(u, v, capacity):
        pipe1 = Pipe(u, v, capacity)
        pipe2 = Pipe(v, u, capacity)
        pipe1.rpipe = pipe2
        pipe2.rpipe = pipe1

        Pipe.pipe_dict[pipe1.id] = len(Pipe.pipe_list)
        Pipe.pipe_list.append(pipe1)
        Pipe.pipe_dict[pipe2.id] = len(Pipe.pipe_list)
        Pipe.pipe_list.append(pipe2)

    @staticmethod
    def find_path(u, v, path=None, nodes=None):
        if path is None:
            path = []
        if nodes is None:
            nodes = []
        if u == v:
            return path
        for i in Node.find_node_by_id(u).neighbor:
            temp_pipe = Pipe.find_pipe_by_id("{}-{}".format(u, i))
            if temp_pipe.residual > 0 and temp_pipe.id not in path and temp_pipe.rpipe.id not in path and i not in nodes:
                result = Pipe.find_path(i, v, path + [temp_pipe.id], nodes + [i])
                if result is not None:
                    return result
        return None

    @staticmethod
    def find_path_breadth(u, v):
        check_list = [u]
        reach_list = [u]
        back_dict = {u: None}
        result = []
        while check_list:
            i = check_list.pop()
            end_flag = False
            for j in Node.find_node_by_id(i).neighbor:
                if Pipe.find_pipe_by_id("{}-{}".format(i, j)).residual == 0:
                    continue
                if j not in back_dict:
                    back_dict[j] = i
                    reach_list.append(j)
                    check_list.append(j)
                    if j == v:
                        end_flag = True
                        break
            if end_flag:
                break
        if v in reach_list:
            t = back_dict[v]
            result = ["{}-{}".format(t, v)] + result
            while t != u:
                result = ["{}-{}".format(back_dict[t], t)] + result
                t = back_dict[t]
            return result
        else:
            return None

    @staticmethod
    def clean_flow():
        for i in Pipe.pipe_list:
            i.flow = 0
            i.residual = i.capacity

    @staticmethod
    def min_cut(u, v):
        # output: 1. max flow between u and v; 2. edge set of cut; 3. groups after cut.
        if u == v:
            raise ValueError("min_cut error: u == v")
        path = Pipe.find_path_breadth(u, v)
        while path is not None:
            residuals = [Pipe.find_pipe_by_id(i).residual for i in path]
            flow = min(residuals)
            for i in path:
                Pipe.find_pipe_by_id(i).add_flow(flow)
            path = Pipe.find_path_breadth(u, v)

        max_flow = sum(Pipe.find_pipe_by_id("{}-{}".format(u, i)).flow for i in Node.find_node_by_id(u).neighbor)
        # print("Max flow from {} to {}: {}".format(u, v, max_flow))
        result = Pipe.find_cut(Pipe.reachable_nodes(u, [u]))

        Pipe.clean_flow()
        return max_flow, result[0], result[1]

    @staticmethod
    def reachable_nodes(u, result=None):
        if result is None:
            result = []
        for i in Node.find_node_by_id(u).neighbor:
            temp_pipe = Pipe.find_pipe_by_id("{}-{}".format(u, i))
            if temp_pipe.residual > 0 and i not in result:
                # result.append(i)
                result = Pipe.reachable_nodes(i, result + [i])
        return result

    @staticmethod
    def find_cut(g1):
        g2 = [i.id for i in Node.node_list if i.id not in g1]
        result = []
        for i in g1:
            for j in g2:
                if "{}-{}".format(i, j) in Pipe.pipe_dict:
                    result.append("{}-{}".format(i, j))
        return result, [g1, g2]

    @staticmethod
    def find_pipe_by_id(id):
        return Pipe.pipe_list[Pipe.pipe_dict[id]]

    @staticmethod
    def clear():
        Pipe.pipe_list = []
        Pipe.pipe_dict = {}

class Tree_node:
    def __init__(self, id):
        self.id = id
        self.neighbor = []
        self.parent = None
        self.fl = 0
        self.dW = 0

    def __repr__(self):
        return "Tree node %s" % self.id

    def modify_neighbor(self, new_n, old_n):
        if old_n is not None:
            self.neighbor.remove(old_n.id)
            old_n.neighbor.remove(self.id)
        if new_n is not None:
            self.neighbor.append(new_n.id)
            new_n.neighbor.append(self.id)


class Equivalent_flow_tree:
    def __init__(self):
        self.node_id_list = [i.id for i in Node.node_list]
        self.node_list = []
        self.node_dict = {}
        self.edge_list = []
        self.edge_dict = {}
        self.tree_type = "None"
        self.delta = 3
        self.Wt = None

    def equivalent_flow_tree(self):
        self.tree_type = "Equivalent flow tree"
        self.node_list = []
        self.node_dict = {}
        self.edge_list = []
        self.edge_dict = {}
        result = []

        root = Tree_node(self.node_id_list[0])
        self.add_node_list(root)

        for n, i in enumerate(self.node_id_list):
            if n == 0:
                continue
            temp_node = Tree_node(i)
            temp_node.modify_neighbor(root, None)
            self.add_node_list(temp_node)

        for n, i in enumerate(self.node_list):
            if n == 0:
                continue
            if len(i.neighbor) != 1:
                raise ValueError("Non-leaf node")
            t = self.find_node_by_id(i.neighbor[0])

            mincut = Pipe.min_cut(i.id, t.id)
            result.append([i.id, t.id, mincut[0]])
            X = None
            if i.id in mincut[2][0]:
                X = mincut[2][0]
            else:
                X = mincut[2][1]
            for j in self.node_list[n + 1:]:
                if len(j.neighbor) != 1:
                    raise ValueError("Non-leaf node")
                if j.id in X and j.neighbor[0] == t.id:
                    j.modify_neighbor(i, self.find_node_by_id(j.neighbor[0]))
        return result

    def mincut_tree(self):
        self.tree_type = "Mincut tree"
        self.node_list = []
        self.node_dict = {}
        self.edge_list = []
        self.edge_dict = {}
        self.Wt = 0
        result = []

        root = Tree_node(self.node_id_list[0])
        self.add_node_list(root)

        for n, i in enumerate(self.node_id_list):
            if n == 0:
                continue
            temp_node = Tree_node(i)
            temp_node.parent = root.id
            self.add_node_list(temp_node)

        for n, i in enumerate(self.node_list):
            if n == 0:
                continue
            t = self.find_node_by_id(i.parent)

            mincut = Pipe.min_cut(i.id, t.id)
            i.fl = mincut[0]
            X = None
            if i.id in mincut[2][0]:
                X = mincut[2][0]
            else:
                X = mincut[2][1]

            for j in self.node_list:
                if j.id != i.id and j.id in X and j.parent == t.id:
                    j.parent = i.id

            if t.parent in X:
                i.parent = t.parent
                t.parent = i.id
                i.fl = t.fl
                t.fl = mincut[0]

        self.tree_completion_by_parent()
        for i in self.node_list:
            if i.parent == None:
                i.dW = sum(self.find_node_by_id(j).fl for j in i.neighbor)
                continue
            temp_edge = Edge(i.id, i.parent, i.fl)
            self.add_edge_list(temp_edge)
            self.Wt += i.fl
            i.dW = sum(self.find_node_by_id(j).fl for j in i.neighbor) - self.find_node_by_id(i.parent).fl + i.fl
            result.append([i.id, i.parent, i.fl])
        self.Wt /= len(self.edge_list)
        # print("Mincut tree:")
        # for i in result:
        #     print("{} - {}: {}".format(i[0], i[1], i[2]))
        # print()
        return result

    def add_node_list(self, node):
        self.node_dict[node.id] = len(self.node_list)
        self.node_list.append(node)

    def add_edge_list(self, edge):
        self.edge_dict[edge.id] = len(self.edge_list)
        self.edge_list.append(edge)

    def find_node_by_id(self, id):
        return self.node_list[self.node_dict[id]]

    def find_edge_by_ids(self, id1, id2):
        id = "{}-{}".format(id1, id2)
        if id in self.edge_dict:
            return self.edge_list[self.edge_dict[id]]
        else:
            id = "{}-{}".format(id2, id1)
            return self.edge_list[self.edge_dict[id]]

    def tree_completion_by_parent(self):
        # use node parent data to complete node neighbors data
        for i in self.node_list:
            if i.parent == None:
                continue
            i.neighbor.append(i.parent)
            self.find_node_by_id(i.parent).neighbor.append(i.id)

    def critical_tree(self):
        if self.tree_type == "None":
            print("Tree is not built")
            return

        max_weight = -1
        max_weight_node = None
        for i in self.node_list:
            if i.dW > max_weight:
                max_weight = i.dW
                max_weight_node = i

        result = [max_weight_node.id]
        check_list = [max_weight_node.id]

        while check_list:
            i = check_list.pop()
            for j in self.find_node_by_id(i).neighbor:
                j = self.find_node_by_id(j)
                if j.dW > self.Wt and len(j.neighbor) >= self.delta and j.id not in result:
                    result.append(j.id)
                    check_list.append(j.id)
        print("Critical nodes:")
        print(" ".join(result))
        return result


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
            temp_edge = Edge(i[0], i[2], 1 / float(i[4]))
            Edge.edge_list.append(temp_edge)


def clear_all():
    Node.clear()
    Edge.clear()
    Pipe.clear()


def mincut_tree_processing(PDB_id):
    clear_all()
    hotspot_read("../data/SKEMPI/graph_hotspot_o_ring/{}.txt".format(PDB_id))

    for i in Edge.edge_list:
        Pipe.add_pipe(i.n1, i.n2, i.w)

    # print(Pipe.min_cut("a", "c"))
    t1 = Equivalent_flow_tree()
    mincut_tree = t1.mincut_tree()
    result = t1.critical_tree()

    with open("mincut tree/{}_mincut_tree.txt".format(PDB_id), "w") as f1:
        for i in mincut_tree:
            f1.write("{}\t{}\t{}\n".format(i[0], i[1], i[2]))

    with open("result/{}_mincut.txt".format(PDB_id), "w") as f1:
        for i in result:
            f1.write("{}\n".format(i))

# mincut_tree_processing("1AO7_o")

graph_files = os.walk("../data/SKEMPI/graph_hotspot_o_ring")
for path, dir_list, file_list in graph_files:
    for file_name in file_list:
        temp_pdb_id = file_name.split(".")[0]
        mincut_tree_processing(temp_pdb_id)
