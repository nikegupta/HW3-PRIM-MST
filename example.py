from mst import Graph
import numpy as np
from sklearn.metrics import pairwise_distances

def main():
    file_path = './data/small.csv'
    g = Graph(file_path)
    g.construct_mst()
    print(g.mst)


if __name__ == "__main__":
    main()