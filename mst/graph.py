import numpy as np
import heapq
from typing import Union


#finds minimum distance between vertex v and mst mst_vertices according to adj_list
def findMinDistance(mst_vertices, v, adj_list):
    min_distance = float('inf')

    #finds value in column of v for all mst_vertices
    for i in range(len(mst_vertices)):
        if adj_list[mst_vertices[i]][v] < min_distance and adj_list[mst_vertices[i]][v] > 0:
            min_distance = adj_list[mst_vertices[i]][v]

    return min_distance

class Graph:

    def __init__(self, adjacency_mat: Union[np.ndarray, str]):
        """
    
        Unlike the BFS assignment, this Graph class takes an adjacency matrix as input. `adjacency_mat` 
        can either be a 2D numpy array of floats or a path to a CSV file containing a 2D numpy array of floats.

        In this project, we will assume `adjacency_mat` corresponds to the adjacency matrix of an undirected graph.
    
        """
        if type(adjacency_mat) == str:
            self.adj_mat = self._load_adjacency_matrix_from_csv(adjacency_mat)
        elif type(adjacency_mat) == np.ndarray:
            self.adj_mat = adjacency_mat
        else: 
            raise TypeError('Input must be a valid path or an adjacency matrix')
        self.mst = None

    def _load_adjacency_matrix_from_csv(self, path: str) -> np.ndarray:
        with open(path) as f:
            return np.loadtxt(f, delimiter=',')

    def construct_mst(self):
        """
    
        TODO: Given `self.adj_mat`, the adjacency matrix of a connected undirected graph, implement Prim's 
        algorithm to construct an adjacency matrix encoding the minimum spanning tree of `self.adj_mat`. 
            
        `self.adj_mat` is a 2D numpy array of floats. Note that because we assume our input graph is
        undirected, `self.adj_mat` is symmetric. Row i and column j represents the edge weight between
        vertex i and vertex j. An edge weight of zero indicates that no edge exists. 
        
        This function does not return anything. Instead, store the adjacency matrix representation
        of the minimum spanning tree of `self.adj_mat` in `self.mst`. We highly encourage the
        use of priority queues in your implementation. Refer to the heapq module, particularly the 
        `heapify`, `heappop`, and `heappush` functions.

        """

        #get number of vertices in graph
        num_v = len(self.adj_mat)
        
        #make list of vertices and edges for mst 
        mst_vertices = []
        mst_edges = []

        #make list of vertices to be added to the tree, 
        #the identity of the vertex is its index and its value in this list is a 2 item list of its
        #closest value to the existing mst, and a bool of whether its in the tree or not,
        #and its index identity
        vertices = []
        for i in range(0,num_v):
            vertices.append([float('inf'),False,i])

        #add first vertex to MST
        mst_vertices.append(0)
        vertices[0] = [0,True,0]

        #make empty priority queue
        pq = []

        #add vertices not in tree to priority queue
        for i in range(num_v):
            if vertices[i][1] == False:
                vertices[i][0] = findMinDistance(mst_vertices,i,self.adj_mat)
                heapq.heappush(pq,vertices[i])
        
        #main loop of prim algorithm
        while len(pq) > 0:

            #pop from priority queue
            v_to_add = heapq.heappop(pq)

            #add vertex and edge to mst
            mst_vertices.append(v_to_add[2])
            mst_edges.append(v_to_add[0])

            #set vertex as in mst
            vertices[v_to_add[2]][1] == True

            #for all vertices not in the tree, update their distance to the tree
            for i in range(num_v):
                if vertices[i][1] == False:
                    vertices[i][0] = findMinDistance(mst_vertices,i,self.adj_mat)

            #recalculate priority queue
            heapq.heapify(pq)

        #turn mst_vertices and mst_edges into a tree
        #make empty adjacency list
        self.mst = np.zeros((num_v,num_v))

        #add edges on both sides of the diagonal
        for i in range(len(mst_vertices) - 1):
            self.mst[mst_vertices[i]][mst_vertices[i+1]] = mst_edges[i]
            self.mst[mst_vertices[i+1]][mst_vertices[i]] = mst_edges[i]

        for i in range(len(mst_edges)):
            if mst_edges[i] == float('inf'):
                self.mst = 'Warning: Graph is not connected'


