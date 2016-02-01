/**
 File name: bc_bottom_up_frontier.cu
 Author: Yuede Ji
 Last update: 20:40 01-14-2016
 Description: Bottom up CPU bc
    (1) read begin position, csr, weight value from binary file
    (2) betweenness centrality
    (3) average sharing frontiers between two different BCs
    (4) the frontiers are counted for two split steps

**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wtime.h"
#include "graph.h"

using namespace std;

const char output_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_cpu.result";
const char sp_count_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_cpu.sp_count";
const char dist_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_cpu.dist";
const char sa_file[] = "/home/yuede/small_graph/result_bc/bc_bottom_up_cpu.sa";

const char fontier_step1[] = "frontier_shared_step1.txt";
const char fontier_step2[] = "frontier_shared_step2.txt";

const int INF = 0x7fffffff;
const int V = 218;

path_t dist[2][V];
index_t sa[2][V];
index_t sp_count[2][V];
path_t bc[V];
path_t frontier[2];
path_t bc_tmp[2][V];

void print_result()
{
    FILE * fp = fopen(output_file, "w");
    for(int i=0; i<V; ++i)
    {
        fprintf(fp, "%d %g\n", i, bc[i]);
        //fprintf(fp, "%g\n", dist[i]);
    }
    fclose(fp);
}
void print_debug(graph *g)
{
    FILE * fp_count = fopen(sp_count_file, "w");
    for(index_t i=0; i<V; ++i)
        fprintf(fp_count, "%u\n", sp_count[i]);
    fclose(fp_count);

    FILE * fp_sa = fopen(sa_file, "w");
    for(index_t i=0; i<V; ++i)
        fprintf(fp_sa, "%u\n", sa[i]);
    fclose(fp_sa);

    FILE * fp_dist = fopen(dist_file, "w");
    for(index_t i=0; i<V; ++i)
        fprintf(fp_dist, "%g\n", dist[i]);
    fclose(fp_dist);
}
index_t sssp(index_t root1, index_t root2, graph *g)
{
    //printf("root1 = %u, root2 = %u\n", root1, root2);
    //cout<<"v = "<<v<<", e = "<<e<<endl;
    memset(sp_count, 0, 2*sizeof(index_t) * V);
    sp_count[0][root1] = 1;
    sp_count[1][root2] = 1;
    for(index_t i=0; i<g->vert_count; ++i)
    {
        dist[0][i] = INF;
        sa[0][i] = INF;
        dist[1][i] = INF;
        sa[1][i] = INF;
    }
    dist[0][root1] = 0;
    sa[0][root1] = 0;
    dist[1][root1] = 0;
    sa[1][root2] = 0;

    index_t level = 0;
    bool flag = true;
    path_t frontier_tmp = 0;
    while(flag)
    {
        flag = false; 
        index_t frontier_joint = 0;
        index_t frontier_shared = 0;
        for(index_t i=0; i<g->vert_count; ++i)
        {
            if(sa[0][i] == INF || sa[1][i] == INF)
            {
                frontier_joint += 1;
                if(sa[0][i] == sa[1][i])
                    frontier_shared += 1;
            }
        }
        frontier_tmp += frontier_shared*1.0/frontier_joint;
        //printf("%u, %u, %g\n", frontier_shared, frontier_joint, frontier_shared*1.0/frontier_joint);
        for(index_t i=0; i<g->vert_count; ++i)
        {
            bool flag_one = false;
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(dist[0][g->csr[j]] < INF)
                {
                    if(dist[0][i] > dist[0][g->csr[j]] + 1)
                    {
                        dist[0][i] = dist[0][g->csr[j]] + 1;
                        sp_count[0][i] = 0;
                        sa[0][i] = sa[0][g->csr[j]] + 1;
                      //  if(sa[0][i] > level)
                      //      level = sa[0][i];
                        if(!flag_one)
                            flag_one = true;
                    }

                }
            }
            if(flag_one)
            {
                flag = true;
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                    if(dist[0][i] == dist[0][g->csr[j]] + 1)
                    {
                        sp_count[0][i] += sp_count[0][g->csr[j]];
                    }
            }

        }
        for(index_t i=0; i<g->vert_count; ++i)
        {
            bool flag_one = false;
            for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
            {
                if(dist[1][g->csr[j]] < INF)
                {
                    if(dist[1][i] > dist[1][g->csr[j]] + 1)
                    {
                        dist[1][i] = dist[1][g->csr[j]] + 1;
                        sp_count[1][i] = 0;
                        sa[1][i] = sa[1][g->csr[j]] + 1;
               //         if(sa[1][i] > level)
               //             level = sa[1][i];
                        if(!flag_one)
                            flag_one = true;
                    }

                }
            }
            if(flag_one)
            {
                flag = true;
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                    if(dist[1][i] == dist[1][g->csr[j]] + 1) 
                    {
                        sp_count[1][i] += sp_count[1][g->csr[j]];
                    }
            }

        }
        ++level;
        //printf("level = %u\n", level);
        //printf("sa[%u] = %lf\n", v-1, sa[v-1]);
    }
    frontier[0] += frontier_tmp/level;
    //printf("%g, %g\n", frontier_tmp, frontier[0]);
    return level;
}
void bc_one(index_t root1, index_t root2, graph *g, index_t level)
{ 
    path_t frontier_tmp = 0;
    memset(bc_tmp, 0, 2*sizeof(path_t)*V);
    for(index_t cur=level; cur>=0; --cur)
    {
        index_t frontier_joint = 0;
        index_t frontier_shared = 0;
        for(index_t i=0; i<g->vert_count; ++i)
        {
            if(sa[0][i] > level || sa[1][i] > level)
            {
                frontier_joint += 1;
                frontier_shared += 1;
            }
            if(sa[0][i] < cur || sa[1][i] < cur)
            {
                frontier_joint += 1;
                if(sa[0][i] == sa[1][i])
                    frontier_shared += 1;
            }
        }
        if(frontier_joint != 0)
            frontier_tmp += frontier_shared*1.0/frontier_joint;
        
        for(index_t i=0; i<g->vert_count; ++i)
        {
            if(sa[0][i] == cur)
            {
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                {
                    index_t w = g->csr[j];
                    if(dist[0][w] == dist[0][i] + 1)
                    {
                        //if(sp_count[w] != 0)
                        bc_tmp[0][i] += sp_count[0][i]*1.0*(1+bc_tmp[0][w])/sp_count[0][w];
                    }

                }

            }
        }
        for(index_t i=0; i<g->vert_count; ++i)
        {
            if(sa[1][i] == cur)
            {
                for(index_t j=g->beg_pos[i]; j<g->beg_pos[i+1]; ++j)
                {
                    index_t w = g->csr[j];
                    if(dist[1][w] == dist[1][i] + 1)
                    {
                        //if(sp_count[w] != 0)
                        bc_tmp[1][i] += sp_count[1][i]*1.0*(1+bc_tmp[1][w])/sp_count[1][w];
                    }

                }

            }
        }
    }
    bc_tmp[0][root1] = 0;
    bc_tmp[1][root2] = 0;
    for(index_t i=0; i<V; ++i)
    {
        bc[i] += bc_tmp[0][i] + bc_tmp[1][i];
    }
    frontier[1] += frontier_tmp/(level+1);
    
    //printf("%g\n", bc[1]);
}
void bc_all(graph *g)
{
    memset(bc, 0, sizeof(path_t)*V);
    frontier[0] = 0;
    frontier[1] = 0;
    index_t combi = 0;
    //index_t level = sssp(0, g);
    for(index_t i=0; i<g->vert_count; ++i)
    {
        for(index_t j=i+1; j<g->vert_count; ++j)
        {
            index_t level = sssp(i, j, g);
            bc_one(i, j, g, level);
            //printf("%u, %u, %u\n", i, j, level);
            ++combi;
            //break;
        }
        //break;
    }
    printf("step 1 frontier sharing ratio is %g\n", frontier[0]/combi);
    printf("step 2 frontier sharing ratio is %g\n", frontier[1]/combi);
}
int main(int args, char ** argv)
{
    printf("Input: ./bfs_small_graph /path/to/beg /path/to/csr thread-count\n");

//    if(args != 5)
//        exit(-1);
    const char *beg_filename = argv[1];
    const char *csr_filename = argv[2];
    const char *weight_filename = argv[3];
    const int thd_cound = atoi(argv[4]);
    
    graph *g = new graph(beg_filename, csr_filename);//, weight_filename);
    
    bc_all(g);
    //sssp(0, g);// g->vert_count, g->edge_count);
    
    print_result();

    print_debug(g);

    return 0;
}
