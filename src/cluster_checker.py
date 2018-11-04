""" checks integrity of clusters """

from fn3client import fn3Client
import networkx as nx                

class cluster_checker():
    """ tests network generation and integrity  """
    def __init__(self):
        
        self.fn3c =fn3Client()         # expect success
    def run(self):
        """ checks all clusters from a fn3 server.
        Tests whether all nodes can be reached from each other
        """

        clustering = self.fn3c.clustering()

        for algorithm in clustering['algorithms']:
            cluster_ids = self.fn3c.cluster_ids(algorithm)
            for cluster_id in cluster_ids:
                network = self.fn3c.network(algorithm, cluster_id)
                #print(">>",algorithm, network['nNodes'], network['nEdges'])
                N = []
                E = []
                for item in network['elements']:
                    if item['group']=='nodes':
                        #print(algorithm, cluster_id, "node",item['data']['id'])
                        N.append(item['data']['id'])
                        
                    if item['group']=='edges':
                        #print(algorithm, cluster_id,'edge',item['data']['source'],item['data']['target'])
                        E.append((item['data']['source'],item['data']['target']))
                
                G= nx.Graph()
                G.add_nodes_from(N)
                G.add_edges_from(E)
                items = list()
                for item in nx.connected_components(G):
                    items.append(item)
                if len(items)==1:
                    pass
                    #print(algorithm, cluster_id, network['nNodes'], network['nEdges'], "OK")
                else:
                    print(algorithm, cluster_id, network['nNodes'], network['nEdges'], "NOT OK: {0} items".format(len(items)))
                    for item in items:
                        print(len(item))

if __name__ == '__main__':
    checker = cluster_checker()
    checker.run()