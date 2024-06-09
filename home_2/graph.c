/*
graph.c

Set of vertices and edges implementation.

Implementations for helper functions for graph construction and manipulation.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2024
*/
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "graph.h"
#include "utils.h"
#include "pq.h"

#define INITIALEDGES 32

struct edge;

/* Definition of a graph. */
struct graph {
  int numVertices;
  int numEdges;
  int allocedEdges;
  struct edge **edgeList;
};

/* Definition of an edge. */
struct edge {
  int start;
  int end;
  int cost;
};

struct graph *newGraph(int numVertices){
  struct graph *g = (struct graph *) malloc(sizeof(struct graph));
  assert(g);
  /* Initialise edges. */
  g->numVertices = numVertices;
  g->numEdges = 0;
  g->allocedEdges = 0;
  g->edgeList = NULL;
  return g;
}

/* Adds an edge to the given graph. */
void addEdge(struct graph *g, int start, int end, int cost){
  assert(g);
  struct edge *newEdge = NULL;
  /* Check we have enough space for the new edge. */
  if((g->numEdges + 1) > g->allocedEdges){
    if(g->allocedEdges == 0){
      g->allocedEdges = INITIALEDGES;
    } else {
      (g->allocedEdges) *= 2;
    }
    g->edgeList = (struct edge **) realloc(g->edgeList,
      sizeof(struct edge *) * g->allocedEdges);
    assert(g->edgeList);
  }

  /* Create the edge */
  newEdge = (struct edge *) malloc(sizeof(struct edge));
  assert(newEdge);
  newEdge->start = start;
  newEdge->end = end;
  newEdge->cost = cost;

  /* Add the edge to the list of edges. */
  g->edgeList[g->numEdges] = newEdge;
  (g->numEdges)++;
}

/* Returns a new graph which is a deep copy of the given graph (which must be 
  freed with freeGraph when no longer used). */
struct graph *duplicateGraph(struct graph *g){
  struct graph *copyGraph = (struct graph *) malloc(sizeof(struct graph));
  assert(copyGraph);
  copyGraph->numVertices = g->numVertices;
  copyGraph->numEdges = g->numEdges;
  copyGraph->allocedEdges = g->allocedEdges;
  copyGraph->edgeList = (struct edge **) malloc(sizeof(struct edge *) * g->allocedEdges);
  assert(copyGraph->edgeList || copyGraph->numEdges == 0);
  int i;
  /* Copy edge list. */
  for(i = 0; i < g->numEdges; i++){
    struct edge *newEdge = (struct edge *) malloc(sizeof(struct edge));
    assert(newEdge);
    newEdge->start = (g->edgeList)[i]->start;
    newEdge->end = (g->edgeList)[i]->end;
    newEdge->cost = (g->edgeList)[i]->cost;
    (copyGraph->edgeList)[i] = newEdge;
  }
  return copyGraph;
}

/* Frees all memory used by graph. */
void freeGraph(struct graph *g){
  int i;
  for(i = 0; i < g->numEdges; i++){
    free((g->edgeList)[i]);
  }
  if(g->edgeList){
    free(g->edgeList);
  }
  free(g);
}

struct solution *graphSolve(struct graph *g, enum problemPart part,
  int numLocations, int startingLocation, int finalLocation){
  struct solution *solution = (struct solution *)
    malloc(sizeof(struct solution));
  assert(solution);
  if(part == PART_A){
    /* IMPLEMENT 2A SOLUTION HERE */
    solution->damageTaken = -1;

    // Dijkstra's algorithm with each edge's cost set to 1
    int i;
    int *verticesCost = (int *) malloc(sizeof(int) * g->numVertices);

    for (i = 0;i < g->numVertices;++i)
      *(verticesCost + i) = -1;

    *(verticesCost + startingLocation) = 0;

    struct pq *edgeCostPQ = newPQ();

    for (i = 0;i < g->numEdges;++i) {
      int start = (*(g->edgeList + i))->start;
      int end   = (*(g->edgeList + i))->end;
      if (start == startingLocation || end == startingLocation)
        enqueue(edgeCostPQ, g->edgeList[i], 0);
    }

    // use a prior which increment once after used can make the pq become a FIFO list.
    int prior = 0;
    while (!empty(edgeCostPQ)) {
      struct edge *pEdge = (struct edge *) deletemin(edgeCostPQ);
      int edgeStart = pEdge->start;
      int edgeEnd = pEdge->end;
      int edgeCost = 1;
      if (*(verticesCost + edgeEnd) == -1) {
        *(verticesCost + edgeEnd) = *(verticesCost + edgeStart) + edgeCost;
        for (i = 0;i < g->numEdges;++i) {
          int start = (*(g->edgeList + i))->start;
          int end   = (*(g->edgeList + i))->end;
          if (start == edgeEnd || end == edgeEnd)
            enqueue(edgeCostPQ, g->edgeList[i], prior++);
        }
      }
      else if (*(verticesCost + edgeStart) == -1) {
        *(verticesCost + edgeStart) = *(verticesCost + edgeEnd) + edgeCost;
        for (i = 0;i < g->numEdges;++i) {
          int start = (*(g->edgeList + i))->start;
          int end   = (*(g->edgeList + i))->end;
          if (start == edgeStart || end == edgeStart)
            enqueue(edgeCostPQ, g->edgeList[i], prior++);
        }
      }
      else {
        if (*(verticesCost + edgeStart) < *(verticesCost + edgeEnd))
          *(verticesCost + edgeEnd) = *(verticesCost + edgeStart) + 1;
        else if (*(verticesCost + edgeStart) > *(verticesCost + edgeEnd))
          *(verticesCost + edgeStart) = *(verticesCost + edgeEnd) + 1;
      }
    }

    solution->damageTaken = *(verticesCost + finalLocation);
    free(verticesCost);
  } else if(part == PART_B) {
    /* IMPLEMENT 2B SOLUTION HERE */
    solution->totalCost = -1;

    // Dijkstra's algorithm with each edge's cost set to 1
    int i;
    int *verticesCost = (int *) malloc(sizeof(int) * g->numVertices);

    for (i = 0;i < g->numVertices;++i)
      *(verticesCost + i) = -1;

    *(verticesCost + startingLocation) = 0;

    struct pq *edgeCostPQ = newPQ();

    for (i = 0;i < g->numEdges;++i) {
      int start = (*(g->edgeList + i))->start;
      int end   = (*(g->edgeList + i))->end;
      if (start == startingLocation || end == startingLocation)
        enqueue(edgeCostPQ, g->edgeList[i], 0);
    }

    // use a prior which increment once after used can make the pq become a FIFO list.
    int prior = 1;
    while (!empty(edgeCostPQ)) {
      struct edge *pEdge = (struct edge *) deletemin(edgeCostPQ);
      int edgeStart = pEdge->start;
      int edgeEnd = pEdge->end;
      int edgeCost = pEdge->cost;
      if (*(verticesCost + edgeEnd) == -1) {
        *(verticesCost + edgeEnd) = *(verticesCost + edgeStart) + edgeCost;
        for (i = 0;i < g->numEdges;++i) {
          int start = (*(g->edgeList + i))->start;
          int end   = (*(g->edgeList + i))->end;
          if (start == edgeEnd || end == edgeEnd)
            enqueue(edgeCostPQ, g->edgeList[i], prior++);
        }
      }
      else if (*(verticesCost + edgeStart) == -1) {
        *(verticesCost + edgeStart) = *(verticesCost + edgeEnd) + edgeCost;
        for (i = 0;i < g->numEdges;++i) {
          int start = (*(g->edgeList + i))->start;
          int end   = (*(g->edgeList + i))->end;
          if (start == edgeStart || end == edgeStart)
            enqueue(edgeCostPQ, g->edgeList[i], prior++);
        }
      }
      else {
        if (*(verticesCost + edgeStart) + edgeCost < *(verticesCost + edgeEnd)) {
          *(verticesCost + edgeEnd) = *(verticesCost + edgeStart) + edgeCost;
          for (i = 0;i < g->numEdges;++i) {
            int start = (*(g->edgeList + i))->start;
            int end   = (*(g->edgeList + i))->end;
            if (start == edgeEnd || end == edgeEnd)
              enqueue(edgeCostPQ, g->edgeList[i], prior++);
          }
        } else if (*(verticesCost + edgeStart) > *(verticesCost + edgeEnd) + edgeCost) {
          *(verticesCost + edgeStart) = *(verticesCost + edgeEnd) + edgeCost;
          for (i = 0;i < g->numEdges;++i) {
            int start = (*(g->edgeList + i))->start;
            int end   = (*(g->edgeList + i))->end;
            if (start == edgeStart || end == edgeStart)
              enqueue(edgeCostPQ, g->edgeList[i], prior++);
          }
        }
      }
    }

    solution->totalCost = *(verticesCost + finalLocation);
    free(verticesCost);

  } else if(part == PART_C) {
    /* IMPLEMENT 2C SOLUTION HERE */
    solution->artisanCost = 0;

    struct pq *edgeCostPQ = newPQ();

    int sumVerticesAdd = 1;
    int i;

    // use the field item to store the edge
    for (i = 0;i < g->numEdges;++i) {
      int start = (*(g->edgeList + i))->start;
      int end   = (*(g->edgeList + i))->end;
      int cost  = (*(g->edgeList + i))->cost;
      if (start == startingLocation || end == startingLocation) {
        enqueue(edgeCostPQ, g->edgeList[i], cost);
      }
    }

    // record whether a vertex is added
    int *verticesAdd = (int *) malloc(g->numVertices * sizeof(int));
    for (i = 0;i < g->numVertices;++i) {
      *(verticesAdd + i) = 0;
    }
    *(verticesAdd + startingLocation) = 1;

    while (sumVerticesAdd < g->numVertices) {
      struct edge * edge = (struct edge *)deletemin(edgeCostPQ);

      int edgeStart = edge->start;
      int edgeEnd = edge->end;

      if (verticesAdd[edgeStart] == 1 && verticesAdd[edgeEnd] == 0) {
        for (i = 0;i < g->numEdges;++i) {
          int start = (*(g->edgeList + i))->start;
          int end   = (*(g->edgeList + i))->end;
          int cost  = (*(g->edgeList + i))->cost;
          if (edge != g->edgeList[i] && (start == edgeEnd || end == edgeEnd)) {
            enqueue(edgeCostPQ, g->edgeList[i], cost);
          }
        }
        verticesAdd[edgeEnd] = 1;
      } else if (verticesAdd[edgeStart] == 0 && verticesAdd[edgeEnd] == 1) {
        for (i = 0;i < g->numEdges;++i) {
          int start = (*(g->edgeList + i))->start;
          int end   = (*(g->edgeList + i))->end;
          int cost  = (*(g->edgeList + i))->cost;
          if (edge != g->edgeList[i] && (start == edgeStart || end == edgeStart)) {
            enqueue(edgeCostPQ, g->edgeList[i], cost);
          }
        }
        verticesAdd[edgeStart] = 1;
      } else {
        // if both vertices are added, continue;
        continue;
      }

      solution->artisanCost += edge->cost;
      sumVerticesAdd++;
    }

  } else {
    /* IMPLEMENT 2D SOLUTION HERE */
    solution->totalPercentage = -1;
  }
  return solution;
}

