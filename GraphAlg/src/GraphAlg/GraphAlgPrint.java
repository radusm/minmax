/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GraphAlg;

import GraphAlg.GraphAlg.ColoredEdge;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.jgrapht.Graph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.interfaces.MatchingAlgorithm;

/**
 *
 * @author radu
 */
public class GraphAlgPrint {

    public static void printTime(long refTime) {
        double time = (System.currentTimeMillis() - refTime) / 1000.0;
        System.out.print("\t" + time);
    }
    public static void printTime(long refTime, long someTime) {
        double time = (someTime - refTime) / 1000.0;
        System.out.print("\t" + time);
    }
    public static void printStats(long startTime, Graph<Integer, GraphAlg.ColoredEdge> G, MatchingAlgorithm.Matching<Integer, GraphAlg.ColoredEdge> M, Graph<Integer, GraphAlg.ColoredEdge> MM) {
        System.out.print("\t" + ((MM == null) ? (M.getEdges().size()) : (MM.edgeSet().size())));
        ConnectivityInspector ci = new ConnectivityInspector(G);
        List<Set<Integer>> list = ci.connectedSets();
        int numempty = 0;
        boolean isEmpty;
        for (Set<Integer> s : list) {
            isEmpty = true;
            for (Integer n : s) {
                if (G.degreeOf(n) > 0) {
                    isEmpty = false;
                    break;
                }
            }
            if (isEmpty) {
                numempty++;
            }
        }
        System.out.print("\t" + (list.size()) /*+ "\t" + numempty*/);
        int max = list.get(0).size();
        int maxdeg = -1;
        for (Integer n : G.vertexSet()) {
            if (maxdeg < G.degreeOf(n)) {
                maxdeg = G.degreeOf(n);
            }
        }

        for (Set<Integer> s : list) {
            if (s.size() > max) {
                max = s.size();
            }
        }

        System.out.print("\t" + (((MM == null) ? (M.getEdges().size()) : (MM.edgeSet().size())) + list.size() - numempty));
        System.out.print("\t" + max);
        //System.out.print("\t" + maxdeg);
        printTime(startTime);
        System.out.print("\n");
    }

    public static void printMinMaxStats(long startTime, Graph<Integer, GraphAlg.ColoredEdge> G, Graph<Integer, GraphAlg.ColoredEdge> M) {
        ConnectivityInspector ci = new ConnectivityInspector(G);
        List<Set<Integer>> list = ci.connectedSets();
        int numConComp = list.size();
        int numempty = 0;
        boolean isEmpty;
        for (Set<Integer> s : list) {
            isEmpty = true;

            for (Integer n : s) {
                if (G.degreeOf(n) > 0) {
                    isEmpty = false;
                    //System.out.println();
                    HashSet<GraphAlg.ColoredEdge> hs = new HashSet<>();
                    for (Integer nn : s) {
                        for (GraphAlg.ColoredEdge edge : G.edgesOf(nn)) {
                            hs.add(edge);
                        }
                    }
                    //System.out.println(hs.size());
                    break;
                }
            }
            if (isEmpty) {
                numempty++;
            }
        }
        System.out.print("\t" + (list.size()) + "\t" + numempty);
        int max = list.get(0).size();
        int maxdeg = -1;
        for (Integer n : G.vertexSet()) {
            if (maxdeg < G.degreeOf(n)) {
                maxdeg = G.degreeOf(n);
            }
        }

        for (Set<Integer> s : list) {
            HashSet<GraphAlg.ColoredEdge> hs = new HashSet<>();
            for (Integer n : s) {
                for (GraphAlg.ColoredEdge e : G.edgesOf(n)) {
                    hs.add(e);
                }
            }

            if (hs.size() > max) {
                max = hs.size();
            }
        }
        int numComp = 0;
        ci = new ConnectivityInspector(M);
        list = ci.connectedSets();

        System.out.print("\t" + (list.size() + numConComp - numempty));
        System.out.print("\t" + max);
        System.out.print("\t" + maxdeg);
        printTime(startTime);
        System.out.print("\n");
    }

    public static void printMinMaxStatsMatchColoring(long startTime, Graph<Integer, GraphAlg.ColoredEdge> G, Graph<Integer, GraphAlg.ColoredEdge> M,DisjointSets forest) {
        ArrayList<Integer> setSizes=new ArrayList<>();
        ConnectivityInspector ci = new ConnectivityInspector(G);
        List<Set<Integer>> list = ci.connectedSets();
        int numConComp = list.size();
        int numempty = 0;
        boolean isEmpty;
        for (Set<Integer> s : list) {
            isEmpty = true;

            for (Integer n : s) {
                if (G.degreeOf(n) > 0) {
                    isEmpty = false;
                    //System.out.println();
                    HashSet<GraphAlg.ColoredEdge> hs = new HashSet<>();
                    for (Integer nn : s) {
                        for (GraphAlg.ColoredEdge edge : G.edgesOf(nn)) {
                            hs.add(edge);
                        }
                    }
                    setSizes.add(hs.size());
                    break;
                }
            }
            if (isEmpty) {
                numempty++;
            }
        }
        System.out.print("\t" + (list.size()- numempty));
        int max = list.get(0).size();

/*
        for (Set<Integer> s : list) {
            HashSet<GraphAlg.ColoredEdge> hs = new HashSet<>();
            for (Integer n : s) {
                for (GraphAlg.ColoredEdge e : G.edgesOf(n)) {
                    hs.add(e);
                }
            }

            if (hs.size() > max) {
                max = hs.size();
            }
        }*/
        int numComp = 0;
        ci = new ConnectivityInspector(M);
        list = ci.connectedSets();

        System.out.print("\t" + (list.size() + numConComp - numempty));
        //System.out.print("\t" + max);

       

        
        HashSet<Integer> hhh = new HashSet<Integer>();
        for (GraphAlg.ColoredEdge e : M.edgeSet()) {
            e.color = forest.find(e.color);
            hhh.add(e.color);
        }
        
        for (int i : hhh) {
            setSizes.add(forest.getNumber(i));
        }
        
        Collections.sort(setSizes);
        Collections.reverse(setSizes);
        System.out.print("\t"+setSizes.get(0));
        
        
        printTime(startTime);  
        System.out.println();
        
        
    } 
    public static void printMinMaxStatsClustering(long startTime, Graph<Integer, GraphAlg.ColoredEdge> G, DisjointSetsWeight forest) {


        //System.out.print("\t" + max);


        
        ArrayList<Integer> clusters=new ArrayList<>();
        
        for (Integer i:G.vertexSet())
            clusters.add(forest.getNumber(i));
        int undercluster=0;
        for (ColoredEdge e:G.edgeSet())
            undercluster+=e.color;
        clusters.add(undercluster);
        
        Collections.sort(clusters);
        Collections.reverse(clusters);
        System.out.print("\t"+clusters.size());
            System.out.print("\t"+clusters.get(0));
              printTime(startTime);
        System.out.println();
        
        
    }
}
