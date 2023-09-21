#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double computePathwayScore(DataFrame expr,
                           NumericVector bfs,
                           List edges,
                           List weights) {
  
  // set the names of nodes in the graph
  CharacterVector nodes = bfs.names();
  
  // extract the number of nodes
  int n = nodes.size();
  
  // extract differential expression results
  CharacterVector allIds = expr["ID"];
  NumericVector lfc = expr["logFC"];
  NumericVector nW = expr["weights"];
  CharacterVector types = expr["type"];
  lfc.attr("names") = allIds;
  nW.attr("names") = allIds;
  types.attr("names") = allIds;
  
  // define weights, edges, expression values and node types for pathway nodes
  NumericVector lfcNodes;
  NumericVector weightNodes;
  CharacterVector nodeTypes;
  List u(n);
  List w(n);
  for (int i = 0; i < n; i++) {
    std::string name = Rcpp::as<std::string>(nodes[i]);
    u[i] = edges[name];
    w[i] = weights[name];
    lfcNodes.push_back(lfc[name], name);
    weightNodes.push_back(nW[name], name);
    nodeTypes.push_back(types[name], name);
  }
  
  // define expression change direction for each node
  IntegerVector sgnNodes = sign(lfcNodes);
  sgnNodes.attr("names") = nodes;
  
  // define pathway score
  double pS = 0;
  
  // define a vector that stores node scores
  NumericVector nodeScores (n);
  nodeScores.attr("names") = nodes;
  
  // iterate over each node
  for (int i = 0; i < n; i++) {
    
    // extract node name
    String node = nodes[i];
    
    // define node score
    double nS = weightNodes[node];
    
    // extract interactions for this node
    CharacterVector interactors = u[i];
    
    // skip nodes without upstream genes
    if (interactors.length() == 0) {
      nodeScores[node] = nS;
      pS += nS;
      continue;
    }
    
    // extract edge weights for this node
    NumericVector edgeW = w[i];
    edgeW.attr("names") = interactors;
    
    // iterate over each interaction
    for (int j = 0; j < interactors.size(); j++) {
      
      // get the name of the upstream gene
      String dest = interactors[j];
      
      // compute the contribution factor
      int beta = edgeW[j];
      int gamma = 0;
      int sgnDown = sgnNodes[node];
      int sgnUp = sgnNodes[dest];
      if (sgnDown == sgnUp) {
        gamma = beta;
      } else {
        gamma = -beta;
      }
      
      // get the weight of the upstream gene
      double k = 0;
      if (nodeScores[dest] == 0) {
        k = weightNodes[dest];
      } else {
        k = nodeScores[dest];
      }
      
      // add the contribution of the upstream gene to the node score
      double upCon = gamma * k;
      nS += upCon;
      
    }
    
    // add the calculated node score to the score vector
    nodeScores[node] = nS;
    
    // add the contribution of the node score to the pathway score
    pS += nS;
    
  }
  
  // determine the proportion of miRNAs
  IntegerVector numTypes = table(nodeTypes);
  double miRs = 0;
  if (numTypes.length() > 1) {
    miRs = numTypes["miRNA"];
  }
  double M = miRs / n;
  
  // normalize pathway score for miRNA proportion and the number of nodes
  pS = pS * ((1 - M) / n);
  
  // return pathway score
  return pS;
  
}

