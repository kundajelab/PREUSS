#Converts a bpRNA annotation to a 2D graph representation
#Generates a nodes and edges output as well as the shortest path between each pair of nodes
import argparse
import pickle
from Graph import *
from dijkstra import *
import pdb 

def parse_args():
    parser=argparse.ArgumentParser(description="Converts a bpRNA annotation to a 2D graph representation.Generates a nodes and edges output as well as the shortest path between each pair of nodes")
    parser.add_argument("--bpRNA_pickle")
    parser.add_argument("--out_prefix")
    parser.add_argument("--approach",default="inferred",choices=["inferred","bootstrap"])
    return parser.parse_args()


def get_graph_representation(bpRNA_entry,approach):
    annotation=bpRNA_entry[approach].split('\n')
    features=annotation[5]
    cur_sequence=annotation[3]
    cur_structure=annotation[4]
    cur_structure_names=annotation[5] 
        
    #populate all edges in the graph
    complement_nodes=get_matched_brackets(cur_structure) 
    g=Graph()
    for node in complement_nodes:
        #annotate the edges 
        g.weights[(node,complement_nodes[node])]=1
        g.edges[node]=[complement_nodes[node]]    
    for i in range(len(cur_sequence)):
        if i < (len(cur_sequence)-1): 
            downstream_node=i+1
            if i not in g.edges:
                g.edges[i]=[downstream_node]
            else:
                g.edges[i].append(downstream_node)
            g.weights[tuple([i,downstream_node])]=1
            g.weights[tuple([downstream_node,i])]=1
        if i > 0:
            upstream_node=i-1
            if i not in g.edges:
                g.edges[i]=[upstream_node]
            else: 
                g.edges[i].append(upstream_node)
            g.weights[tuple([i,upstream_node])]=1
            g.weights[tuple([upstream_node,i])]=1
        #annotate the nodes
        g.nodes[i]=(cur_sequence[i],cur_structure_names[i])
    return g

def get_shortest_paths_for_all_nodes(g,cur_id,distance_dict,edge_dict):
    edge_dict[cur_id]=[]
    distance_dict[cur_id]=dict() 
    #get shortest path from each node to all other nodes;
    for i in list(g.nodes.keys()):        
        connections_i=g.edges[i]
        for c in connections_i:
            edge_dict[cur_id].append((i,c))
        for j in list(g.nodes.keys()):
            if i==j:
                continue
            shortest_path_i_j,length_shortest_path_i_j=dijkstra(g,i,j)
            distance_dict[cur_id][(i,j)]=length_shortest_path_i_j
    return distance_dict,edge_dict

def get_matched_brackets(cur_structure):
    edge_dict=dict()
    open=[]
    for i in range(len(cur_structure)):
        cur_character=cur_structure[i]
        if cur_character==".":
            continue 
        if cur_character=="(":
            open.append(i)
        elif cur_character==")":
            last_open=open.pop()
            edge_dict[last_open]=i
            edge_dict[i]=last_open
    assert len(open)==0
    return edge_dict 

def main():
    args=parse_args()
    bprna_data=pickle.load(open(args.bpRNA_pickle,'rb'))
    #for each isoform, store the shortest distance between each pair of bases 
    distance_dict=dict()
    #for each isoform, generate a graph indicating all connections between bases 
    edge_dict=dict()
    #for each isoform, generate a dictionary of all nodes : node_id -> (sequence,structure)
    node_dict=dict() 
    
    for cur_id in bprna_data:
        print("analyzing id:"+str(cur_id))
        cur_graph=get_graph_representation(bprna_data[cur_id],args.approach)
        node_dict[cur_id]=cur_graph.nodes        
        distance_dict,edge_dict=get_shortest_paths_for_all_nodes(cur_graph,cur_id,distance_dict,edge_dict)
        
    print("writing outputs") 
    #write the outputs
    outf_distance=open(args.out_prefix+".distance.tsv",'w')
    outf_edges=open(args.out_prefix+".edges.tsv",'w')
    outf_nodes=open(args.out_prefix+".nodes.tsv",'w')

    for cur_id in distance_dict:
        for entry in distance_dict[cur_id]:
            outf_distance.write(str(cur_id)+
                                '\t'+
                                str(entry[0]+1)+
                                '\t'+
                                str(entry[1]+1)+
                                '\t'+
                                str(distance_dict[cur_id][entry])+
                                '\n')
    for cur_id in edge_dict:
        for entry in edge_dict[cur_id]:
            outf_edges.write(str(cur_id)+'\t'+str(entry[0]+1)+'\t'+str(entry[1]+1)+'\n')
                
    for cur_id in node_dict:
        for entry in node_dict[cur_id]:
            cur_seq_base=node_dict[cur_id][entry][0]
            cur_structure=node_dict[cur_id][entry][1] 
            outf_nodes.write(str(cur_id)+'\t'+str(entry+1)+'\t'+cur_seq_base+'\t'+cur_structure+'\n')
            

if __name__=="__main__":
    main() 

