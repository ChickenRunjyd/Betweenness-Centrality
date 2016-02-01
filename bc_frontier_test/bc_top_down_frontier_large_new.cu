/**
 File name: bc_top_down_frontier_large.cu
 Author: Yuede Ji
 Last update: 22:52 01-12-2016
 Description: CPU bc on small graph
    (1) read begin position, csr, weight value from binary file
    (2) betweenness centrality from 0 to others
    (3) top down sharing frontiers test
    (4) one step, one frontier file

**/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"

using namespace std;

const char output_file[] = "/home/yuede/small_graph/result_bc/bc_sa_cpu.result";
const char debug_file[] = "/home/yuede/small_graph/result_bc/bc_sa_cpu.debug";
const char frontier_file_1[] = "top_down_frontier.txt";
const char frontier_file_2[] = "top_down_frontier.txt";

const int INF = 0x7fffffff;
const unsigned int V = 61578405;
/*
path_t bc[V];
index_t sa[2][V];
index_t sp_count[2][V];
path_t dist[2][V];
path_t bc_root[2][V];

path_t frontier[2];
void print_result()
{
    FILE * fp = fopen(output_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%d %g\n", i, bc[i]);
    }
    fclose(fp);
}
*/
/*
void print_debug(graph *g)
{
    FILE * fp = fopen(debug_file, "w");
    for(index_t j=0; j<V; ++j)
    {
        fprintf(fp, "\n%u:\t", j);
        for(index_t i=g->beg_pos[j]; i<g->beg_pos[j+1]; ++i)
        {
            fprintf(fp, "%u %lf\t", g->csr[i], g->weight[i]); 
        }
    }
    fclose(fp);
}
*/
index_t sssp(index_t root1, index_t root2, graph *g, index_t *sa, index_t *sp_count, path_t *dist, path_t *bc_root, path_t *frontier)
{
    index_t step = g->vert_count;
    index_t level = 0;
    for(index_t i=0; i<g->vert_count; ++i)
    {
        dist[i] = INF;
        dist[step + i] = INF;
    }
    dist[root1] = 0;
    dist[step + root2] = 0;
    memset(sa, -1, 2*sizeof(index_t)*g->vert_count);
    memset(sp_count, 0, 2*sizeof(index_t)*g->vert_count);
    sp_count[root1] = 1;
    sp_count[step + root2] = 1;
    sa[root1] = 0;
    sa[step + root2] = 0;
    bool flag = true;
    path_t frontier_tmp = 0;
    
    while(flag)
    {
        //sharing frontiers check
        index_t frontier_joint = 0;
        index_t frontier_shared = 0;
        for(index_t j=0; j<g->vert_count; ++j)
        {
            if(sa[j] == level || sa[step + j] == level )
            {
                frontier_joint += 1;
                if(sa[j] == sa[step + j])
                    frontier_shared += 1;
            }
        }
        //printf("shared = %d, joint = %d\n", frontier_shared, frontier_joint);
        if(frontier_joint != 0)
            frontier_tmp += frontier_shared * 1.0/frontier_joint;
        printf("%g\n", frontier_shared * 1.0/frontier_joint);
        
        flag = false;
        
        for(index_t j=0; j<g->vert_count; ++j)
        {
            if(sa[j] != level)
                continue;
            for(index_t i=g->beg_pos[j]; i<g->beg_pos[j+1]; ++i)
            {
                index_t v = g->csr[i];
                
                //std::cout<<"v = "<<v<<"\n";
                //printf("weight[%u] = %g\n", v, g->weight[v]);
                if(dist[v] > dist[j] + 1)
                {
                    //printf("v = %u, j = %u\n", v, j);
                    dist[v] = dist[j] + 1;
                    flag = true;
                    sa[v] = level + 1;
                    sp_count[v] = 0;
                }
                //printf("v = %u\n", v);
                if(dist[v] == dist[j] + 1)
                    sp_count[v] += sp_count[j];
            }
        }

        for(index_t j=0; j<g->vert_count; ++j)
        {
            if(sa[step + j] != level)
                continue;
            for(index_t i=g->beg_pos[j]; i<g->beg_pos[j+1]; ++i)
            {
                index_t v = g->csr[i];
                //printf("weight[%u] = %g\n", v, g->weight[v]);
                if(dist[step + v] > dist[step + j] + 1)
                {
                    dist[step + v] = dist[step + j] + 1;
                    flag = true;
                    sa[step + v] = level + 1;
                    sp_count[step + v] = 0;
                }
                if(dist[step + v] == dist[step + j] + 1)
                    sp_count[step + v] += sp_count[step + j];
            }
        }
//        printf("level = %d, frontier_tmp = %g\n", level, frontier_tmp);
        ++level;
    }
    printf("level = %d, frontier_tmp = %g\n", level, frontier_tmp);
    frontier[0] += frontier_tmp/level;
    return level;
}
void bc_one(index_t root1, index_t root2, graph *g, index_t level, index_t *sa, index_t *sp_count, path_t *dist, path_t *bc_root, path_t *frontier, path_t *bc)
{
    index_t step = g->vert_count;
    memset(bc_root, 0, 2*sizeof(path_t)*step);
    //reverse the order
    path_t frontier_tmp = 0;
    for(index_t cur=level-1; cur>=0; --cur)
    {
        index_t frontier_joint = 0;
        index_t frontier_shared = 0;
        for(index_t j=0; j<g->vert_count; ++j)
        {
            if(sa[j] == cur || sa[step + j] == cur)
            {
                frontier_joint += 1;
                if(sa[j] == sa[step + j])
                    frontier_shared += 1;
            }
        }
        if(frontier_joint != 0)
            frontier_tmp += frontier_shared * 1.0/frontier_joint;
        
        printf("%g\n", frontier_shared * 1.0/frontier_joint);
        for(index_t i=0; i<g->vert_count; ++i)
        {
            if(sa[i] == cur)
            {
                index_t w = i;
                //printf("g->beg_pos[%u] = %u\n", w, g->beg_pos[w]);
                //Undirected graph
                for(index_t j=g->beg_pos[w]; j<g->beg_pos[w+1]; ++j)
                {
                    index_t v = g->csr[j];
                    //printf("v = %u\n", v);
                    if(dist[w] == dist[v] + 1)
                    {
                        if(sp_count[w] != 0)
                            bc_root[v] += sp_count[v]*1.0/sp_count[w]*(1+bc_root[w]); 
                    }
                }
                if(w != root1)///cur > 0 can avoid this branch?
                    bc[w] += bc_root[w];
            }
        }
        for(index_t i=0; i<g->vert_count; ++i)
        {
            if(sa[step + i] == cur)
            {
                index_t w = i;
                //printf("g->beg_pos[%u] = %u\n", w, g->beg_pos[w]);
                //Undirected graph
                for(index_t j=g->beg_pos[w]; j<g->beg_pos[w+1]; ++j)
                {
                    index_t v = g->csr[j];
                    //printf("v = %u\n", v);
                    if(dist[step + w] == dist[step + v] + 1)
                    {
                        if(sp_count[step + w] != 0)
                            bc_root[step + v] += sp_count[step + v]*1.0/sp_count[step + w]*(1+bc_root[step + w]); 
                    }
                }
                if(w != root2)///cur > 0 can avoid this branch?
                    bc[w] += bc_root[step + w];
            }
        }
    }
    if(level > 0) 
        frontier[1] += frontier_tmp/level;
}

void bc_all(graph *g)
{
    path_t *bc = new path_t[g->vert_count+1];
    index_t *sa = new index_t[(g->vert_count+1)*2];
    index_t *sp_count = new index_t[(g->vert_count+1)*2];
    path_t *dist = new path_t[(g->vert_count+1)*2];
    path_t *bc_root = new path_t[(g->vert_count+1)*2];
    path_t *frontier = new path_t[2];
    
    memset(bc, 0, sizeof(path_t)*g->vert_count);
    frontier[0] = 0;
    frontier[1] = 0;
    index_t combi = 0;
    index_t root1 = 0;
    srand(time(NULL));
    index_t root2 = rand() % g->vert_count;
    printf("root1 = %u, root2 = %u\n", root1, root2);
    index_t depth = sssp(root1, root2, g, sa, sp_count, dist, bc_root, frontier);
    bc_one(root1, root2, g, depth, sa, sp_count, dist, bc_root, frontier, bc);
    ++combi;
    printf("step 1 frontier sharing ratio is %g\n", frontier[0]/combi);
    printf("step 2 frontier sharing ratio is %g\n", frontier[1]/combi);
    
}
int main(int args, char ** argv)
{
    printf("Input: ./bc /path/to/beg /path/to/csr vertex-count\n");

  //  if(args != 5)
  //      exit(-1);
    const char *beg_filename = argv[1];
    const char *csr_filename = argv[2];
  //  const char *weight_filename = argv[3];
  //  const int thd_count = atoi(argv[4]);
    graph *g = new graph(beg_filename, csr_filename);
//    for(int i=0; i<10; ++i)
//        std::cout<<g->csr[i]<<"\n";
//    for(int i=0; i<10; ++i)
//        std::cout<<g->beg_pos[i]<<"\n";
    bc_all(g);
// g->vert_count, g->edge_count);
//    
//    print_result();
//
//    print_debug(g);
//
    return 0;
}
