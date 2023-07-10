#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double computePathwayScore(NumericVector lfc,
                           List edges,
                           List weights,
                           int nN) {
  
  // extract the names of nodes in the graph
  CharacterVector nodes = edges.names();
  
  // remove miRNAs and genes not present in the pathway
  CharacterVector allIds = lfc.names();
  NumericVector lfcNodes;
  for (int i = 0; i < lfc.size(); i++) {
    std::string elementName = as<std::string>(allIds[i]);
    if (std::find(nodes.begin(), nodes.end(), elementName) != nodes.end()) {
      lfcNodes.push_back(lfc[i], elementName);
    }
  }
  
  // extract gene IDs and the number of nodes
  CharacterVector ids = lfcNodes.names();
  int n = lfcNodes.size();
  
  // reorder edges and weights according to IDs vector
  int m = edges.size();
  List u(m);
  List w(m);
  
  for (int i = 0; i < n; i++) {
    std::string name = Rcpp::as<std::string>(ids[i]);
    u[i] = edges[name];
    w[i] = weights[name];
  }
  
  // define pathway score
  double pS = 0;
  
  // iterate over each node
  for (int i = 0; i < n; i++) {
    
    // extract node name
    String node = ids[i];
    
    // extract interactions and edge weights for this node
    CharacterVector interactors = u[i];
    NumericVector edgeW = w[i];
    edgeW.attr("names") = interactors;
    
    // iterate over each node for the second time
    for (int j = 0; j < n; j++) {
      
      // get the name of the upstream gene
      String dest = ids[j];
      
      // determine if the two nodes interact in this pathway
      bool interaction = false;
      for (int nInt = 0; nInt < interactors.size(); nInt ++) {
        if (interactors[nInt] == dest) {
          interaction = true;
        }
      }
      
      // compute the score for each node
      if (interaction == true) {
        pS += lfcNodes[i] * lfcNodes[j] * edgeW[dest];
      }
    }
  }
  
  // normalize pathway score for the number of nodes
  pS /= nN;
  
  // return pathway score
  return pS;
  
}

