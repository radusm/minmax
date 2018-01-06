package GraphAlg;

/*
 * (C) Copyright 2017-2017, by Joris Kinable and Contributors.
 *
 * JGraphT : a free Java graph-theory library
 *
 * This program and the accompanying materials are dual-licensed under
 * either
 *
 * (a) the terms of the GNU Lesser General Public License version 2.1
 * as published by the Free Software Foundation, or (at your option) any
 * later version.
 *
 * or (per the licensee's choosing)
 *
 * (b) the terms of the Eclipse Public License v1.0 as published by
 * the Eclipse Foundation.
 */
import static GraphAlg.GraphAlgPrint.*;
import com.mxgraph.layout.*;

import com.mxgraph.swing.*;
import com.mxgraph.util.*;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.*;
import java.io.*;
import java.text.DecimalFormat;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.alg.interfaces.MatchingAlgorithm;
import org.jgrapht.alg.interfaces.MatchingAlgorithm.Matching;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.*;
import org.jgrapht.alg.matching.EdmondsMaximumCardinalityMatching;

import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.util.Pair;
import javax.swing.JFrame;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.jgrapht.ext.JGraphXAdapter;

import org.w3c.dom.Document;

/**
 * Tests for EdmondsMaximumCardinalityMatching
 *
 * @author Joris Kinable
 */
public final class GraphAlg {

    public static int readEdgesFromFile(Graph<Integer, ColoredEdge> g, String filename) {

        File file = new File(filename);
        Scanner sc = null;
        try {
            sc = new Scanner(file);
        } catch (FileNotFoundException ex) {
            Logger.getLogger(GraphAlg.class.getName()).log(Level.SEVERE, null, ex);
        }
        if (sc == null) {
            return -1;
        }
        while (sc.hasNext("#.*")) {
            sc.nextLine();
        }
        int vertexLabel = 1;
        HashMap<Integer, Integer> map = new HashMap<>();
        while (sc.hasNextInt()) {
            int a = sc.nextInt();
            if (!map.containsKey(a)) {
                map.put(a, vertexLabel);
                a = vertexLabel;
                vertexLabel++;
            } else {
                a = map.get(a);
            }
            int b = sc.nextInt();
            if (!map.containsKey(b)) {
                map.put(b, vertexLabel);
                b = vertexLabel;
                vertexLabel++;
            } else {
                b = map.get(b);
            }
            if (a != b) {
                g.addVertex(a);
                g.addVertex(b);
                g.addEdge(a, b/*,new ColoredEdge()*/);
            }
        }
        sc.close();

        ConnectivityInspector ci = new ConnectivityInspector(g);
        List<Set> l = ci.connectedSets();
        int max = -1;
        Set maxSet = null;
        for (Set s : l) {
            if (max < s.size()) {
                max = s.size();
                maxSet = s;
            }
        }
        if (maxSet == null) {
            return -1;
        }
        HashSet<Integer> copy = new HashSet<>();
        copy.addAll(g.vertexSet());
        for (int i : copy) {
            if (!maxSet.contains(i)) {
                g.removeVertex(i);
            }
        }

        int maxdeg = -1;
        double frac = 1.0 / g.vertexSet().size();
        double avgdeg = 0;
        for (Integer n : g.vertexSet()) {
            if (maxdeg < g.degreeOf(n)) {
                maxdeg = g.degreeOf(n);
            }
            avgdeg += g.degreeOf(n) * frac;
        }

        DecimalFormat df = new DecimalFormat("#.00");
        //System.out.print(filename/*.substring(0,filename.length()-4)*/ + "\t" + g.vertexSet().size());
        //System.out.print("\t" + g.edgeSet().size());
        //System.out.print("\t" + df.format(avgdeg) + "\t" + maxdeg);
        //System.out.println();
        //System.out.println(g.vertexSet());

        return vertexLabel;
    }

    public static class ColoredEdge extends DefaultEdge {

        //int id1;
        //int id2;
        public int color;

        public ColoredEdge() {
            color = 0;
        }

        @Override
        public boolean equals(Object v) {
            if (v == this) {
                return true;
            }
            if (v == null) {
                return false;
            }
            if (v instanceof ColoredEdge) {
                return super.equals(v);
                //return (((ColoredEdge)v).id1==id1)&&(((ColoredEdge)v).id2==id2);
            }
            return false;
        }

        //might not need this
        @Override
        public int hashCode() {
            /*int hash = 5;
            hash = 89 * hash + this.id1;
            hash = 89 * hash + this.id2;
            return hash;*/
            return super.hashCode();
        }

        @Override
        public String toString() {
            return /*super.toString() + " " +*/ String.valueOf(color);
        }
    };

    public static void testGraphStandard(String filename) {
        int currentColor = 1;
        Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
        readEdgesFromFile(graph, filename);
        long time = System.currentTimeMillis();

        MatchingAlgorithm<Integer, ColoredEdge> matcher = new EdmondsMaximumCardinalityMatching<>(graph);
        Matching<Integer, ColoredEdge> m = matcher.getMatching();

        for (ColoredEdge e : m.getEdges()) {
            e.color = currentColor;
            currentColor++;
            //graph.removeEdge(e);
        }
        printStats(time, graph, m, null);
        //displayGraph(graph, -1);
    }

    public static void testGraphGreedy(String filename) {
        int currentColor = 1;
        Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
        readEdgesFromFile(graph, filename);
        long time = System.currentTimeMillis();
        Graph<Integer, ColoredEdge> m = new SimpleGraph<>(ColoredEdge.class);
        Set<ColoredEdge> edgeset = graph.edgeSet();
        for (ColoredEdge e : edgeset) {
            if (!m.containsVertex(graph.getEdgeSource(e)) && !m.containsVertex(graph.getEdgeTarget(e))) {
                m.addVertex(graph.getEdgeSource(e));
                m.addVertex(graph.getEdgeTarget(e));
                m.addEdge(graph.getEdgeSource(e), graph.getEdgeTarget(e), e);
                //System.out.println(currentColor+": "+n+" "+nn + " "+ nn.color);
                currentColor++;
            }
        }
        for (ColoredEdge e : m.edgeSet()) {
            graph.removeEdge(e);
        }
        printStats(time, graph, null, m);
    }

    public static class EdgeDegComparator implements Comparator<ColoredEdge> {

        Graph<Integer, ColoredEdge> comparingGraph;

        public EdgeDegComparator(Graph<Integer, ColoredEdge> g) {
            comparingGraph = g;
        }

        @Override
        public int compare(ColoredEdge t, ColoredEdge t1) {
            return (comparingGraph.degreeOf(comparingGraph.getEdgeSource(t)) + comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t)))
                    - (comparingGraph.degreeOf(comparingGraph.getEdgeSource(t1)) + comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t1)));
            /*
            int mindeg1=comparingGraph.degreeOf(comparingGraph.getEdgeSource(t));
            mindeg1=comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t))<mindeg1?comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t)):mindeg1;
            
            int mindeg2=comparingGraph.degreeOf(comparingGraph.getEdgeSource(t1));
            mindeg2=comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t1))<mindeg2?comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t1)):mindeg2;
            
            return +mindeg1
                    - ( mindeg2 );
             */
        }

    };

    public static class EdgeDegForestComparator implements Comparator<ColoredEdge> {

        Graph<Integer, ColoredEdge> comparingGraph;
        DisjointSets forest;
        HashMap<Integer, Integer> colors;

        public EdgeDegForestComparator(Graph<Integer, ColoredEdge> g, DisjointSets d, HashMap<Integer, Integer> c) {
            comparingGraph = g;
            forest = d;
            colors = c;
        }

        @Override
        public int compare(ColoredEdge t, ColoredEdge t1) {
            int v1 = comparingGraph.getEdgeSource(t), v2 = comparingGraph.getEdgeTarget(t),
                    v3 = comparingGraph.getEdgeSource(t1), v4 = comparingGraph.getEdgeTarget(t1);
            Integer col1, col2, col3, col4;

            col1 = colors.get(v1);
            if (col1 == null) {
                col1 = 0;
            }
            col1 = forest.find(col1);
            col2 = colors.get(v2);
            if (col2 == null) {
                col2 = 0;
            }
            col2 = forest.find(col2);
            col3 = colors.get(v3);
            if (col3 == null) {
                col3 = 0;
            }
            col3 = forest.find(col3);
            col4 = colors.get(v4);
            if (col4 == null) {
                col4 = 0;
            }
            col4 = forest.find(col4);

            //PRIORITY1: UNMATCHED VTX (COLOR 0)
            int magicResult
                    = -(((col1 == 0) || (col2 == 0)) ? 1 : 0)
                    + (((col3 == 0) || (col4 == 0)) ? 1 : 0);
            if (magicResult != 0) {
                return magicResult;
            }
            //PRIORITY2: EDGES WITH ENDS IN SAME COLOR,smaller set first
            if (col1 == col2 && col3 == col4) {
                return +forest.getNumber(col1)
                        - forest.getNumber(col3);
            }

            magicResult
                    = -((col1 == col2) ? 100 : 0)
                    + ((col3 == col4) ? 100 : 0);
            if (magicResult != 0) {
                return magicResult;
            }
            //PRIORITY3: EDGES BINDING THE SMALLEST COLOR SETS
            magicResult
                    = +(forest.getNumber(col1) + forest.getNumber(col2))
                    - (forest.getNumber(col3) + forest.getNumber(col4));
            if (magicResult != 0) {
                return magicResult;
            }
            //PRIORITY4: EDGES WITH THE SMALLEST SUM OF DEGREES OF ENDS
            /*
            int mindeg1=comparingGraph.degreeOf(comparingGraph.getEdgeSource(t));
            mindeg1=comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t))<mindeg1?comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t)):mindeg1;
            
            int mindeg2=comparingGraph.degreeOf(comparingGraph.getEdgeSource(t1));
            mindeg2=comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t1))<mindeg2?comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t1)):mindeg2;
            
            return +mindeg1
                    - ( mindeg2 );
            
             */

            return +(comparingGraph.degreeOf(comparingGraph.getEdgeSource(t)) + comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t)))
                    - (comparingGraph.degreeOf(comparingGraph.getEdgeSource(t1)) + comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t1)));

        }

    };

    public static class EdgeDegWeightedComparator implements Comparator<ColoredEdge> {

        Graph<Integer, ColoredEdge> comparingGraph;
        DisjointSetsWeight weights;

        public EdgeDegWeightedComparator(Graph<Integer, ColoredEdge> g, DisjointSetsWeight d) {
            comparingGraph = g;
            weights = d;
        }

        @Override
        public int compare(ColoredEdge t, ColoredEdge t1) {
            //PRIORITY1:EDGES WITH ENDS IN THE SMALLEST CLUSTERS
            int magicResult
                    = +(weights.getNumber(comparingGraph.getEdgeSource(t)) + weights.getNumber(comparingGraph.getEdgeTarget(t)))
                    - (weights.getNumber(comparingGraph.getEdgeSource(t1)) + weights.getNumber(comparingGraph.getEdgeTarget(t1)));
            if (magicResult != 0) {
                return magicResult;
            }
            //PRIORITY2:EDGES THAT REPRESENT THE BIGGEST (?) NUMBER OF EDGES
            magicResult = -(t.color - t1.color);
            if (magicResult != 0) {
                return magicResult;
            }
            //PRIORITY3: EDGES WITH THE BIGGEST SUM OF DEGREES OF ENDS
            return -(comparingGraph.degreeOf(comparingGraph.getEdgeSource(t)) + comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t)))
                    + (comparingGraph.degreeOf(comparingGraph.getEdgeSource(t1)) + comparingGraph.degreeOf(comparingGraph.getEdgeTarget(t1)));
        }

    };

    public static void testGraphGreedySort(String filename) {
        int currentColor = 1;
        Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
        readEdgesFromFile(graph, filename);
        long time = System.currentTimeMillis();

        Graph<Integer, ColoredEdge> m = new SimpleGraph<>(ColoredEdge.class);

        Set<ColoredEdge> edgeset = graph.edgeSet();
        ArrayList<ColoredEdge> sortset = new ArrayList<ColoredEdge>();
        sortset.addAll(edgeset);
        Collections.sort(sortset, new EdgeDegComparator(graph));

        for (ColoredEdge e : sortset) {
            if (!m.containsVertex(graph.getEdgeSource(e)) && !m.containsVertex(graph.getEdgeTarget(e))) {
                m.addVertex(graph.getEdgeSource(e));
                m.addVertex(graph.getEdgeTarget(e));
                m.addEdge(graph.getEdgeSource(e), graph.getEdgeTarget(e), e);
                //System.out.println(currentColor+": "+n+" "+nn + " "+ nn.color);
                currentColor++;
            }
        }
        for (ColoredEdge e : m.edgeSet()) {
            graph.removeEdge(e);
        }
        printStats(time, graph, null, m);
    }
    static int it = 0;
    static int wh = 0;
    static JFrame frame;

    static float f1 = 0.2f;
    static float f2 = 0.8f;
    static Random r = new Random(123);

    public static String newColor() {

        float f3 = r.nextFloat() * 0.8f + 0.1f;
        if (Math.abs(f3 - 0.5) < 0.1) {
            if (f3 < 0.5) {
                f3 *= f3;
            } else {
                f3 = 1.0f - f3 * f3;
            }
        }

        if (f2 + f1 < 0.5) {
            f3 = 1.0f - f3 * f3;
        } else if (f3 + f2 < 0.5) {
            f3 = 0.05f;
        }

        String s1 = ("0" + Integer.toHexString((int) (f1 * 255)));
        s1 = s1.substring(s1.length() - 2);
        String s2 = ("0" + Integer.toHexString((int) (f2 * 255)));
        s2 = s2.substring(s2.length() - 2);
        String s3 = ("0" + Integer.toHexString((int) (f3 * 255)));
        s3 = s3.substring(s3.length() - 2);
        f1 = f2;
        f2 = f3;
        return "#" + s1 + s2 + s3;
    }
    static HashMap<Integer, String> cols = new HashMap<>();

    public static void displayGraph(Graph<Integer, ColoredEdge> graph, int N, boolean svgOutput) {

        JGraphXAdapter<Integer, ColoredEdge> jgxAdapter = new JGraphXAdapter<>(graph);
        mxGraphComponent comp = new mxGraphComponent(jgxAdapter);
        frame.getContentPane().removeAll();
        frame.getContentPane().add(comp);

        if (N < 0) {
            mxFastOrganicLayout layout = new mxFastOrganicLayout(jgxAdapter);

            layout.setForceConstant(1500);
            layout.execute(jgxAdapter.getDefaultParent());
        } else {

            for (int i = 0; i < N; i++) {

                for (int j = 1; j <= N; j++) {
                    if (jgxAdapter.getVertexToCellMap().containsKey(i * N + j)) {
                        jgxAdapter.getVertexToCellMap().get(i * N + j).getGeometry().setX(100 * j);
                        jgxAdapter.getVertexToCellMap().get(i * N + j).getGeometry().setY(100 * i);
                        jgxAdapter.getVertexToCellMap().get(i * N + j).getGeometry().setHeight(35);
                        jgxAdapter.getVertexToCellMap().get(i * N + j).getGeometry().setWidth(35);
                    }
                }
            }
        }
        int currentColor = 0;
        HashMap<Integer, ArrayList<Object>> hm = new HashMap<>();

        for (ColoredEdge e : graph.edgeSet()) {
            if (hm.containsKey(e.color)) {
                hm.get(e.color).add(jgxAdapter.getEdgeToCellMap().get(e));
            } else {
                ArrayList<Object> a = new ArrayList<>();
                a.add(jgxAdapter.getEdgeToCellMap().get(e));
                hm.put(e.color, a);
                currentColor++;
            }
        }

        for (int i : hm.keySet()) {
            int ccc = jgxAdapter.getCellToEdgeMap().get(hm.get(i).get(0)).color;
            if (!cols.containsKey(ccc)) {
                cols.put(jgxAdapter.getCellToEdgeMap().get(hm.get(i).get(0)).color, ccc == 0 ? "#E4E4E4" : newColor());
            }
            jgxAdapter.setCellStyles(mxConstants.STYLE_STROKECOLOR, cols.get(jgxAdapter.getCellToEdgeMap().get(hm.get(i).get(0)).color), hm.get(i).toArray());
            jgxAdapter.setCellStyles(mxConstants.STYLE_STROKEWIDTH, "35", hm.get(i).toArray());
            //jgxAdapter.setCellStyles(mxConstants.STYLE_STARTARROW, mxConstants.NONE, hm.get(i).toArray());
            jgxAdapter.setCellStyles(mxConstants.STYLE_STARTARROW, "startSize=12;endSize=12;", hm.get(i).toArray());

            jgxAdapter.setCellStyles(mxConstants.STYLE_ENDARROW, mxConstants.NONE, hm.get(i).toArray());
            jgxAdapter.setCellStyles(mxConstants.STYLE_ENDARROW, "startSize=12;endSize=12;", hm.get(i).toArray());

        }
        comp.zoomTo(2, true);
        Object lol = new Object();
        //jgxAdapter.setCellStyles(mxConstants.,style, clrs)
        comp.getGraphControl().addMouseWheelListener(new MouseAdapter() {
            int lastx = 0, lasty = 0;
            Object cell;

            @Override
            public void mouseWheelMoved(MouseWheelEvent e) {
                int x = e.getX(), y = e.getY();
                if (x != lastx || y != lasty) {
                    cell = comp.getCellAt(x, y);
                }
                if (e.getWheelRotation() > 0) {
                    comp.zoomOut();
                } else {
                    comp.zoomIn();
                }

                if (cell != null) {
                    comp.scrollCellToVisible(cell, true);
                }
                lastx = x;
                lasty = y;
            }
        });
        comp.getGraphControl().addMouseListener(new MouseAdapter() {
            @Override
            public void mouseReleased(MouseEvent e) {
                synchronized (lol) {
                    lol.notify();
                }
            }
        });
        if (svgOutput) {
            Document doc = mxCellRenderer.createSvgDocument(jgxAdapter, null, 1, Color.white, jgxAdapter.getBoundingBoxFromGeometry(jgxAdapter.getChildVertices(jgxAdapter.getDefaultParent())));
            TransformerFactory tFactory
                    = TransformerFactory.newInstance();
            Transformer transformer = null;
            try {
                transformer = tFactory.newTransformer();
            } catch (TransformerConfigurationException ex) {
                Logger.getLogger(GraphAlg.class.getName()).log(Level.SEVERE, null, ex);
            }

            DOMSource source = new DOMSource(doc);
            Random rr = new Random(System.currentTimeMillis());
            BufferedWriter bf = null;
            try {
                bf = new BufferedWriter(new FileWriter("out" + it + "." + rr.nextInt(1000) + ".svg"));
            } catch (IOException ex) {
                Logger.getLogger(GraphAlg.class.getName()).log(Level.SEVERE, null, ex);
            }
            StreamResult result = new StreamResult(bf);
            try {
                transformer.transform(source, result);
            } catch (TransformerException ex) {
                Logger.getLogger(GraphAlg.class.getName()).log(Level.SEVERE, null, ex);
            }
            try {
                bf.close();
            } catch (IOException ex) {
                Logger.getLogger(GraphAlg.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        frame.repaint();
        frame.setTitle("Graph Visualisation" + (it++));
        frame.pack();
        frame.setVisible(true);
        try {
            synchronized (lol) {
                lol.wait();
            }
        } catch (InterruptedException ex) {
            Logger.getLogger(GraphAlg.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
    static int sortCounter = 0;

    public static int buildMaximalMatching(Graph<Integer, ColoredEdge> graph, Graph<Integer, ColoredEdge> result, HashMap<Integer, Integer> colors) {//initial block
        int currentColor = 1;
        ArrayList<ColoredEdge> sortset = new ArrayList<ColoredEdge>();
        sortset.addAll(graph.edgeSet());

        Collections.sort(sortset, new EdgeDegComparator(graph));
        Graph<Integer, ColoredEdge> m = new SimpleGraph<>(ColoredEdge.class);
        for (ColoredEdge e : sortset) {
            Integer src = graph.getEdgeSource(e);
            Integer tgt = graph.getEdgeTarget(e);

            if ((!m.containsVertex(src)) && (!m.containsVertex(tgt))) {
                colors.put(src, currentColor);
                colors.put(tgt, currentColor);
                m.addVertex(src);
                m.addVertex(tgt);
                e.color = currentColor;
                m.addEdge(src, tgt).color = currentColor;
                graph.removeEdge(e);
                currentColor++;
            }
        }

        for (ColoredEdge e : m.edgeSet()) {
            ColoredEdge n = new ColoredEdge();
            n.color = e.color;
            result.addEdge(m.getEdgeSource(e), m.getEdgeTarget(e), n);
        }

        //System.out.println();
        //System.out.print("\t" + m.edgeSet().size());
        //System.out.println(graph.edgeSet());
        return currentColor;
    }//end initial block

    public static void displayGr(Graph<Integer, ColoredEdge> graph, Graph<Integer, ColoredEdge> result, DisjointSets forest, int currentColor, int N) {

        int cc = currentColor + 1;
        Graph<Integer, ColoredEdge> copy = new SimpleGraph<>(ColoredEdge.class);
        for (int i : graph.vertexSet()) {
            copy.addVertex(i);
        }
        for (ColoredEdge e : graph.edgeSet()) {
            ColoredEdge n = new ColoredEdge();
            n.color = e.color;
            copy.addEdge(graph.getEdgeSource(e), graph.getEdgeTarget(e), n);
        }

        ConnectivityInspector ci = new ConnectivityInspector(copy);
        List<Set<Integer>> list = ci.connectedSets();

        for (Set<Integer> s : list) {
            for (int i : s) {
                for (ColoredEdge e : copy.edgesOf(i)) {
                    e.color = cc;
                }
            }
            cc++;
        }

        HashSet<Integer> hhh = new HashSet<Integer>();
        for (ColoredEdge e : result.edgeSet()) {
            e.color = forest.find(e.color);
            hhh.add(e.color);
        }
        System.out.println();
        for (int i : hhh) {
            System.out.println(forest.getNumber(i) + " at color " + i);
        }
        System.out.println();
        for (ColoredEdge e : result.edgeSet()) {
            ColoredEdge n = new ColoredEdge();
            n.color = e.color;
            copy.addEdge(result.getEdgeSource(e), result.getEdgeTarget(e), n);
        }
        displayGraph(copy, N, true);
    }

    public static void displayGrFinal(Graph<Integer, ColoredEdge> graph, Graph<Integer, ColoredEdge> result, int currentColor, int N) {
        Graph<Integer, ColoredEdge> copy = new SimpleGraph<>(ColoredEdge.class);
        for (int i : graph.vertexSet()) {
            copy.addVertex(i);
        }
        for (ColoredEdge e : graph.edgeSet()) {
            ColoredEdge n = new ColoredEdge();
            n.color = e.color;
            copy.addEdge(graph.getEdgeSource(e), graph.getEdgeTarget(e), n);
        }

        ConnectivityInspector ci = new ConnectivityInspector(copy);
        List<Set<Integer>> list = ci.connectedSets();

        for (Set<Integer> s : list) {
            for (int i : s) {
                for (ColoredEdge e : copy.edgesOf(i)) {
                    e.color = currentColor;
                }
            }
            currentColor++;
        }

        Graphs.addGraph(copy, result);
        displayGraph(copy, N, true);
    }

    public static void testMinMaxMatchColoring(String filename, boolean displaySteps, boolean display, int N) {
        /*
        builds a maximal matching and uses it as an initial solution;
        then: sorts the edges after a heuristic attractiveness of edges;
        the attractiveness was the attempt to select edges that are likely
        to disconnect the uncolored graph: but this is actually hard you need
        something like separators to do reliably
        after that: cycle through the sorted edges and find first one that
        does not increase the maximum color size; use this edge to unite 2 colors
         */

        sortCounter = 0;

        Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
        readEdgesFromFile(graph, filename);
        long time = System.currentTimeMillis();
        int M = -1;
        for (int i : graph.vertexSet()) {
            if (i > M) {
                M = i;
            }
        }
        HashMap<Integer, Integer> colors = new HashMap<>(2 * M);

        Graph<Integer, ColoredEdge> result = new SimpleGraph<>(ColoredEdge.class);
        Graphs.addAllVertices(result, graph.vertexSet());

        //build initial solution
        int currentColor = buildMaximalMatching(graph, result, colors);
        //start from a maximal matching

        DisjointSets forest = new DisjointSets(currentColor);

        ArrayList<ColoredEdge> sortset = new ArrayList<ColoredEdge>();
        sortset.addAll(graph.edgeSet());

        boolean somethingChanged = true;
        while (somethingChanged) {
            if (displaySteps) {
                displayGr(graph, result, forest, currentColor, N);
            }

            if (sortCounter % 5 == 0) {
                Collections.sort(sortset, new EdgeDegForestComparator(graph, forest, colors));
            }
            sortCounter++;
            somethingChanged = false;

            for (ColoredEdge e : sortset) {
                Integer src = graph.getEdgeSource(e);
                Integer tgt = graph.getEdgeTarget(e);
                //System.out.println("&&&&&" + src + " " + tgt);
                int srcColor = (colors.get(src) == null) ? 0 : colors.get(src);
                int tgtColor = (colors.get(tgt) == null) ? 0 : colors.get(tgt);
                int color1 = forest.find(srcColor);
                int color2 = forest.find(tgtColor);

                //System.out.println(";;" + color1 + " " + color2);
                //System.out.println(".." + srcColor + " " + tgtColor);
                int max;
                HashSet<ColoredEdge> hs = new HashSet<>();
                ConnectivityInspector ci = new ConnectivityInspector(graph);
                Set<Integer> s = ci.connectedSetOf(src);
                for (Integer n : s) {
                    hs.addAll(graph.edgesOf(n));
                }

                max = hs.size();
                max = forest.maxnumber > max ? forest.maxnumber : max;
                if (color1 == color2) {
                    if (forest.getNumber(color1) + 1 >= max) {
                        continue;
                    }
                } else {
                    if (color1 == 0) {
                        if (forest.getNumber(color2) + 1 >= max) {
                            continue;
                        }
                    }
                    if (color2 == 0) {
                        if (forest.getNumber(color1) + 1 >= max) {
                            continue;
                        }
                    }
                    if (forest.getNumber(color1) + forest.getNumber(color2) + 1 >= max) {
                        continue;
                    }
                }
                //System.out.println("\tADD combining successful!" +forest.getNumber(color1)+" "+forest.getNumber(color2)+ " "+max);
                int resultColor;
                if (color1 == 0) {
                    resultColor = forest.union(color2, color2);//this handles the number inside
                } else if (color2 == 0) {
                    resultColor = forest.union(color1, color1);//this handles the number inside
                } else {
                    resultColor = forest.union(color1, color2);//this handles the number inside
                }
                colors.put(src, resultColor);
                colors.put(tgt, resultColor);

                e.color = resultColor;
                result.addEdge(src, tgt, e);
                sortset.remove(e);
                graph.removeEdge(e);
                somethingChanged = true;
                //System.out.println(nn);
                break;

            }

        }

        printMinMaxStatsMatchColoring(time, graph, result, forest);

        if (display) {
            displayGrFinal(graph, result, currentColor, N);
        }

    }

    public static void testMinMaxMatchColoringV2(String filename, boolean displaySteps, boolean display, int N) {
        /*
        builds a maximal matching and uses it as an initial solution;
        then: sorts the edges after a heuristic attractiveness of edges;
        the attractiveness was the attempt to select edges that are likely
        to disconnect the uncolored graph: but this is actually hard you need
        something like separators to do reliably
        after that: cycle through the sorted edges and find the one that disconnects
        the underlying connected component into balanced pieces
         */

        sortCounter = 0;

        Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
        readEdgesFromFile(graph, filename);
        long time = System.currentTimeMillis();
        int M = -1;
        for (int i : graph.vertexSet()) {
            if (i > M) {
                M = i;
            }
        }
        HashMap<Integer, Integer> colors = new HashMap<>(2 * M);

        Graph<Integer, ColoredEdge> result = new SimpleGraph<>(ColoredEdge.class);
        Graphs.addAllVertices(result, graph.vertexSet());

        //build initial solution
        int currentColor = buildMaximalMatching(graph, result, colors);
        //start from a maximal matching

        DisjointSets forest = new DisjointSets(currentColor);

        ArrayList<ColoredEdge> sortset = new ArrayList<ColoredEdge>();
        sortset.addAll(graph.edgeSet());

        float best = 1;
        ColoredEdge bestEdge = null;
        do {
            if (displaySteps) {
                displayGr(graph, result, forest, currentColor, N);
            }

            //if (sortCounter % 5 == 0) {
            //Collections.sort(sortset, new EdgeDegForestComparator(graph, forest, colors));
            //}
            sortCounter++;

            best = 0;
            bestEdge = null;
            for (ColoredEdge e : sortset) {
                float b = 0;
                Integer src = graph.getEdgeSource(e);
                Integer tgt = graph.getEdgeTarget(e);
                //System.out.println("&&&&&" + src + " " + tgt);
                int srcColor = (colors.get(src) == null) ? 0 : colors.get(src);
                int tgtColor = (colors.get(tgt) == null) ? 0 : colors.get(tgt);
                int color1 = forest.find(srcColor);
                int color2 = forest.find(tgtColor);

                //System.out.println(";;" + color1 + " " + color2);
                //System.out.println(".." + srcColor + " " + tgtColor);
                int max;
                HashSet<ColoredEdge> hs = new HashSet<>();
                ConnectivityInspector ci = new ConnectivityInspector(graph);
                Set<Integer> s = ci.connectedSetOf(src);
                for (Integer n : s) {
                    hs.addAll(graph.edgesOf(n));
                }

                max = hs.size();
                max = forest.maxnumber > max ? forest.maxnumber : max;
                if (color1 == color2) {
                    if (forest.getNumber(color1) + 1 >= max) {
                        continue;
                    }
                    b = 10.0f - 0.001f * forest.getNumber(color1);
                    if (b > best) {
                        best = b;
                        bestEdge = e;
                    }
                } else {
                    if (color1 == 0) {
                        if (forest.getNumber(color2) + 1 >= max) {
                            continue;
                        }
                        bestEdge = e;
                        break;
                    }
                    if (color2 == 0) {
                        if (forest.getNumber(color1) + 1 >= max) {
                            continue;
                        }
                        bestEdge = e;
                        break;
                    }
                    if (forest.getNumber(color1) + forest.getNumber(color2) + 1 >= max) {
                        continue;
                    }

                }
                graph.removeEdge(e);

                if (!ci.pathExists(src, tgt)) {
                    /* hs.clear();
                    s = ci.connectedSetOf(src);
                    for (Integer n : s) {
                        hs.addAll(graph.edgesOf(n));
                    }
                    int max1 = hs.size();

                    hs.clear();
                    s = ci.connectedSetOf(tgt);
                    for (Integer n : s) {
                        hs.addAll(graph.edgesOf(n));
                    }
                    int max2 = hs.size();
                    if (max1 == 0) {
                        b = 1000000.0f-forest.getNumber(color1) - forest.getNumber(color2);
                    } else if (max2 == 0) {
                        b = 1000000.0f-forest.getNumber(color1) - forest.getNumber(color2);
                    } else if (max1 > max2) {
                        b = (float) max2 * ((float) max2 / max1)*10000000.0f-(forest.getNumber(color1) - forest.getNumber(color2));
                    } else {*/
                    b = /*(float) max1 * ((float) max1 / max2)*/ 10000000.0f - (forest.getNumber(color1) - forest.getNumber(color2));
                    //}
                } else {
                    b = 1000000.0f - forest.getNumber(color1) - forest.getNumber(color2);
                }

                graph.addEdge(src, tgt, e);
                if (b > best) {
                    best = b;
                    bestEdge = e;
                }

            }

            if (bestEdge != null) {
                Integer src = graph.getEdgeSource(bestEdge);
                Integer tgt = graph.getEdgeTarget(bestEdge);
                //System.out.println("&&&&&" + src + " " + tgt);
                int srcColor = (colors.get(src) == null) ? 0 : colors.get(src);
                int tgtColor = (colors.get(tgt) == null) ? 0 : colors.get(tgt);
                int color1 = forest.find(srcColor);
                int color2 = forest.find(tgtColor);

                int resultColor;
                if (color1 == 0) {
                    resultColor = forest.union(color2, color2);//this handles the number inside
                } else if (color2 == 0) {
                    resultColor = forest.union(color1, color1);//this handles the number inside
                } else {
                    resultColor = forest.union(color1, color2);//this handles the number inside
                }
                colors.put(src, resultColor);
                colors.put(tgt, resultColor);

                bestEdge.color = resultColor;
                result.addEdge(src, tgt, bestEdge);
                sortset.remove(bestEdge);
                graph.removeEdge(bestEdge);

            }

        } while (bestEdge != null);

        printMinMaxStatsMatchColoring(time, graph, result, forest);

        if (display) {
            displayGrFinal(graph, result, currentColor, N);
        }

    }

    public static void testMinMaxColoringClustering(String filename, boolean displaySteps, boolean display, int N) {

        Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
        readEdgesFromFile(graph, filename);
        long time = System.currentTimeMillis();

        int M = -1;
        for (int i : graph.vertexSet()) {
            if (i > M) {
                M = i;
            }
        }
        DisjointSetsWeight forest = new DisjointSetsWeight(M + 1);

        Graph<Integer, ColoredEdge> backup = new SimpleGraph<>(ColoredEdge.class);
        for (int i : graph.vertexSet()) {
            backup.addVertex(i);
        }
        for (ColoredEdge e : graph.edgeSet()) {
            ColoredEdge n = new ColoredEdge();
            n.color = e.color;
            backup.addEdge(graph.getEdgeSource(e), graph.getEdgeTarget(e), n);
        }

        int totalWeight = graph.edgeSet().size();
        for (ColoredEdge edge : graph.edgeSet()) {
            edge.color = 1;
        }
        boolean OK = true;

        PriorityQueue<ColoredEdge> queue = new PriorityQueue<>(new EdgeDegWeightedComparator(graph, forest));
        queue.addAll(graph.edgeSet());
        while (OK) {
            PriorityQueue<ColoredEdge> queue2 = new PriorityQueue<>(new EdgeDegWeightedComparator(graph, forest));
            for (ColoredEdge eeee : queue) {
                queue2.add(eeee);
            }
            queue = queue2;
            if (displaySteps) {
                Graph<Integer, ColoredEdge> copy = new SimpleGraph<>(ColoredEdge.class);
                for (int i : backup.vertexSet()) {
                    copy.addVertex(i);
                }
                for (ColoredEdge e : backup.edgeSet()) {
                    ColoredEdge n = new ColoredEdge();
                    n.color = e.color;
                    copy.addEdge(backup.getEdgeSource(e), backup.getEdgeTarget(e), n);
                }

                int currColor = 1;
                HashMap<Integer, Integer> hm = new HashMap<>();
                for (ColoredEdge e : copy.edgeSet()) {
                    int a = copy.getEdgeSource(e);
                    int b = copy.getEdgeTarget(e);
                    a = forest.find(a);
                    if (!hm.containsKey(a)) {
                        hm.put(a, currColor);
                        a = currColor;
                        currColor++;
                    } else {
                        a = hm.get(a);
                    }
                    b = forest.find(b);
                    if (!hm.containsKey(b)) {
                        hm.put(b, currColor);
                        b = currColor;
                        currColor++;
                    } else {
                        b = hm.get(b);
                    }
                    if (a != b) {
                        e.color = 0;
                    } else {
                        e.color = a;
                    }
                }

                displayGraph(copy, N, true);
            }

            //displayGraph(graph);
            //ArrayList<ColoredEdge> sortset = new ArrayList<ColoredEdge>();
            ColoredEdge e = queue.poll();
            //if (graph.vertexSet().size()%100==0) System.out.println(graph.vertexSet().size());
            if (e == null) {
                break;
            }

//            System.out.print(graph.getEdgeSource(e) + ":" + graph.getEdgeTarget(e) + ")\t" + (forest.getNumber(graph.getEdgeSource(e)) + forest.getNumber(graph.getEdgeTarget(e)) + e.color) + "##");
//            for (Object o : queue.toArray()) {
//                ColoredEdge eee = (ColoredEdge) o;
//                System.out.print("[" +/*graph.getEdgeSource(eee)+ ":"+ graph.getEdgeTarget(eee)+")\t"*/ +(forest.getNumber(graph.getEdgeSource(eee)) + forest.getNumber(graph.getEdgeTarget(eee)) + eee.color) + "]");
//            }
            //System.out.println();
            Integer src = graph.getEdgeSource(e);
            Integer tgt = graph.getEdgeTarget(e);
            int val = forest.getNumber(src) + forest.getNumber(tgt) + e.color;
            if (val >= totalWeight) {
                OK = false;
                /*for (Integer i : graph.vertexSet()) {
                    System.out.println("cluster edge count: " + forest.getNumber(i));
                }
                System.out.println("undercluster edge count: " + totalWeight);
                 */
                break;
            }
            if (displaySteps) {
                for (Integer i : graph.vertexSet()) {
                    System.out.println("cluster edge count: " + forest.getNumber(i));
                }
                System.out.println("undercluster edge count: " + totalWeight);
                System.out.println();
            }
            //here
            // collapse the edge [src,tgt]
            Set<ColoredEdge> srcAdj = graph.outgoingEdgesOf(src);
            HashSet<Integer> srcNodeAdj = new HashSet<>();
            for (ColoredEdge edge : srcAdj) {
                srcNodeAdj.add(Graphs.getOppositeVertex(graph, edge, src));
            }
            srcNodeAdj.remove(tgt);
            HashSet<ColoredEdge> tgtAdj = new HashSet<>();
            tgtAdj.addAll(graph.outgoingEdgesOf(tgt));
            for (ColoredEdge edge : tgtAdj) {
                int op = Graphs.getOppositeVertex(graph, edge, tgt);
                if (op == src) {
                    continue;
                }
                if (srcNodeAdj.contains(op)) {
                    ColoredEdge toModify = graph.getEdge(src, op);
                    queue.remove(toModify);
                    toModify.color += edge.color;
                    queue.offer(graph.getEdge(src, op));
                    graph.removeEdge(edge);
                    queue.remove(edge);
                } else {
                    ColoredEdge newEdge = new ColoredEdge();
                    newEdge.color = edge.color;
                    graph.addEdge(src, op, newEdge);
                    graph.removeEdge(edge);
                    queue.remove(edge);
                    queue.offer(graph.getEdge(src, op));
                }
            }
            //System.out.println(queue);
            graph.removeEdge(e);
            graph.removeVertex(tgt);
            forest.union(src, tgt, val);
            totalWeight -= e.color;

        }
        printMinMaxStatsClustering(time, graph, forest);
        if (display) {
            Graph<Integer, ColoredEdge> copy = new SimpleGraph<>(ColoredEdge.class);
            for (int i : backup.vertexSet()) {
                copy.addVertex(i);
            }
            for (ColoredEdge e : backup.edgeSet()) {
                ColoredEdge n = new ColoredEdge();
                n.color = e.color;
                copy.addEdge(backup.getEdgeSource(e), backup.getEdgeTarget(e), n);
            }

            int currColor = 1;
            HashMap<Integer, Integer> hm = new HashMap<>();
            for (ColoredEdge e : copy.edgeSet()) {
                int a = copy.getEdgeSource(e);
                int b = copy.getEdgeTarget(e);
                a = forest.find(a);
                if (!hm.containsKey(a)) {
                    hm.put(a, currColor);
                    a = currColor;
                    currColor++;
                } else {
                    a = hm.get(a);
                }
                b = forest.find(b);
                if (!hm.containsKey(b)) {
                    hm.put(b, currColor);
                    b = currColor;
                    currColor++;
                } else {
                    b = hm.get(b);
                }
                if (a != b) {
                    e.color = 0;
                } else {
                    e.color = a;
                }
            }

            displayGraph(copy, N, true);
        }
        //displayGraph(graph);
    }
    static int currColor = 1;

    public static class ret {

        Graph<Integer, ColoredEdge> g;
        int r;
        int bigcluster;
    };

    public static ret testMinMaxColoringClusteringV2(String filename, Graph<Integer, ColoredEdge> G, boolean mybool, boolean display, int N, int zeroCount) {
        ret rrr = new ret();
        Graph<Integer, ColoredEdge> graph;
        if (G == null) {
            graph = new SimpleGraph<>(ColoredEdge.class);
            readEdgesFromFile(graph, filename);
        } else {
            graph = G;
        }

        long time = System.currentTimeMillis();

        int M = -1;
        for (int i : graph.vertexSet()) {
            if (i > M) {
                M = i;
            }
        }
        DisjointSetsWeight forest = new DisjointSetsWeight(M + 1);

        Graph<Integer, ColoredEdge> backup = new SimpleGraph<>(ColoredEdge.class);

        for (int i : graph.vertexSet()) {
            backup.addVertex(i);
        }
        for (ColoredEdge e : graph.edgeSet()) {
            ColoredEdge n = new ColoredEdge();
            n.color = e.color;
            backup.addEdge(graph.getEdgeSource(e), graph.getEdgeTarget(e), n);
        }

        int totalWeight = graph.edgeSet().size();
        for (ColoredEdge edge : graph.edgeSet()) {
            edge.color = 1;
        }
        boolean OK = true;

        PriorityQueue<ColoredEdge> queue = new PriorityQueue<>(new EdgeDegWeightedComparator(graph, forest));
        queue.addAll(graph.edgeSet());
        while (OK) {
            PriorityQueue<ColoredEdge> queue2 = new PriorityQueue<>(new EdgeDegWeightedComparator(graph, forest));
            for (ColoredEdge eeee : queue) {
                queue2.add(eeee);
            }
            queue = queue2;

            //displayGraph(graph);
            //ArrayList<ColoredEdge> sortset = new ArrayList<ColoredEdge>();
            ColoredEdge e = queue.poll();
            //if (graph.vertexSet().size()%100==0) System.out.println(graph.vertexSet().size());
            if (e == null) {
                break;
            }

//            System.out.print(graph.getEdgeSource(e) + ":" + graph.getEdgeTarget(e) + ")\t" + (forest.getNumber(graph.getEdgeSource(e)) + forest.getNumber(graph.getEdgeTarget(e)) + e.color) + "##");
//            for (Object o : queue.toArray()) {
//                ColoredEdge eee = (ColoredEdge) o;
//                System.out.print("[" +/*graph.getEdgeSource(eee)+ ":"+ graph.getEdgeTarget(eee)+")\t"*/ +(forest.getNumber(graph.getEdgeSource(eee)) + forest.getNumber(graph.getEdgeTarget(eee)) + eee.color) + "]");
//            }
            //System.out.println();
            Integer src = graph.getEdgeSource(e);
            Integer tgt = graph.getEdgeTarget(e);
            int val = forest.getNumber(src) + forest.getNumber(tgt) + e.color;
            if (graph.edgeSet().size() == 1) {
                rrr.bigcluster = forest.getNumber(src) > forest.getNumber(tgt) ? forest.getNumber(src) : forest.getNumber(tgt);
                OK = false;
                /*for (Integer i : graph.vertexSet()) {
                    System.out.println("cluster edge count: " + forest.getNumber(i));
                }
                System.out.println("undercluster edge count: " + totalWeight);
                 */
                break;
            }

            //here
            // collapse the edge [src,tgt]
            Set<ColoredEdge> srcAdj = graph.outgoingEdgesOf(src);
            HashSet<Integer> srcNodeAdj = new HashSet<>();
            for (ColoredEdge edge : srcAdj) {
                srcNodeAdj.add(Graphs.getOppositeVertex(graph, edge, src));
            }
            srcNodeAdj.remove(tgt);
            HashSet<ColoredEdge> tgtAdj = new HashSet<>();
            tgtAdj.addAll(graph.outgoingEdgesOf(tgt));
            for (ColoredEdge edge : tgtAdj) {
                int op = Graphs.getOppositeVertex(graph, edge, tgt);
                if (op == src) {
                    continue;
                }
                if (srcNodeAdj.contains(op)) {
                    ColoredEdge toModify = graph.getEdge(src, op);
                    queue.remove(toModify);
                    toModify.color += edge.color;
                    queue.offer(graph.getEdge(src, op));
                    graph.removeEdge(edge);
                    queue.remove(edge);
                } else {
                    ColoredEdge newEdge = new ColoredEdge();
                    newEdge.color = edge.color;
                    graph.addEdge(src, op, newEdge);
                    graph.removeEdge(edge);
                    queue.remove(edge);
                    queue.offer(graph.getEdge(src, op));
                }
            }
            //System.out.println(queue);
            graph.removeEdge(e);
            graph.removeVertex(tgt);
            forest.union(src, tgt, val);
            totalWeight -= e.color;

        }

        printMinMaxStatsClustering(time, graph, forest);
        if (display) {
            Graph<Integer, ColoredEdge> copy = new SimpleGraph<>(ColoredEdge.class);
            for (int i : backup.vertexSet()) {
                copy.addVertex(i);
            }
            for (ColoredEdge e : backup.edgeSet()) {
                ColoredEdge n = new ColoredEdge();
                n.color = e.color;
                copy.addEdge(backup.getEdgeSource(e), backup.getEdgeTarget(e), n);
            }

            HashMap<Integer, Integer> hm = new HashMap<>();
            for (ColoredEdge e : copy.edgeSet()) {
                int a = copy.getEdgeSource(e);
                int b = copy.getEdgeTarget(e);
                a = forest.find(a);
                if (!hm.containsKey(a)) {
                    hm.put(a, currColor);
                    a = currColor;
                    currColor++;
                } else {
                    a = hm.get(a);
                }
                b = forest.find(b);
                if (!hm.containsKey(b)) {
                    hm.put(b, currColor);
                    b = currColor;
                    currColor++;
                } else {
                    b = hm.get(b);
                }
                if (a != b) {
                    e.color = 0;
                } else {
                    e.color = a;
                }
            }

            //displayGraph(copy, N);
            HashMap<Integer, List<ColoredEdge>> hhhhhhh = new HashMap();
            int ctrr = 0;
            for (ColoredEdge e : copy.edgeSet()) {
                if (e.color == 0) {
                    ctrr++;
                    continue;
                }
                if (!hhhhhhh.containsKey(e.color)) {
                    hhhhhhh.put(e.color, new ArrayList<ColoredEdge>());
                    hhhhhhh.get(e.color).add(e);
                } else {
                    hhhhhhh.get(e.color).add(e);
                }
            }

            for (int i : hhhhhhh.keySet()) {
                Graph<Integer, ColoredEdge> temp = new SimpleGraph<>(ColoredEdge.class);
                for (ColoredEdge e : hhhhhhh.get(i)) {
                    temp.addVertex(temp.getEdgeSource(e));
                    temp.addVertex(temp.getEdgeTarget(e));
                    temp.addEdge(temp.getEdgeSource(e), temp.getEdgeTarget(e));
                }
                if (temp.edgeSet().size() > 40) {
                    ret rr = testMinMaxColoringClusteringV2("", temp, false, true, 20, ctrr);

                    //displayGraph(res, N);
                    //split only if the biggest cluster is bigger than the sum of 0 edges, not very interesting, looks like first version
                    //if (zeroCount+rr.r<rr.bigcluster){
                    for (ColoredEdge e : rr.g.edgeSet()) {
                        copy.getEdge(rr.g.getEdgeSource(e), rr.g.getEdgeTarget(e)).color = e.color;
                    }
                    zeroCount += rr.r;
                    //}
                }
            }
            if (mybool) {
                displayGraph(copy, N, true);

            }

            rrr.g = copy;
            rrr.r = zeroCount;
            return rrr;
        }
        return null;
        //displayGraph(graph);
    }

    public static int DFC(int src, Graph<Integer, ColoredEdge> graph, int currentColor, int targetColor) {
        int t = 0;
        Stack<Integer> s = new Stack();
        s.push(src);
        while (!s.empty()) {
            src = s.pop();
            for (ColoredEdge e : graph.edgesOf(src)) {
                if (e.color == targetColor) {
                    int other = Graphs.getOppositeVertex(graph, e, src);
                    e.color = currentColor;
                    t++;
                    s.push(other);
                }
            }
        }
        return t;
    }

    public static int DFS(int src, int tgt, Graph<Integer, ColoredEdge> graph, int targetColor, DisjointSets forest, int N) {
        Stack<Integer> s = new Stack();

        char viz[] = new char[N + 2];
        for (int i = 0; i < N + 2; i++) {
            viz[i] = 0;
        }

        targetColor = forest.find(targetColor);
        int r = 0;
        int test = 1;
        for (ColoredEdge e : graph.edgesOf(src)) {
            if (forest.find(e.color) == targetColor) {
                test = 0;
                break;
            }
        }
        if (test == 1) {
            r = -1;
            //return -1;
        }
        test = 1;
        for (ColoredEdge e : graph.edgesOf(tgt)) {
            if (forest.find(e.color) == targetColor) {
                test = 0;
                break;
            }
        }

        if (test == 1) {
            if (r == 0) {
                return -2;
            }
            return -3;
        } else if (r == -1) {
            return -1;
        }

        s.push(src);
        viz[src] = 1;
        while (!s.empty()) {
            src = s.pop();
            for (ColoredEdge e : graph.edgesOf(src)) {
                if (forest.find(e.color) == targetColor) {
                    int other = Graphs.getOppositeVertex(graph, e, src);
                    if (other == tgt) {
                        return 1;
                    }
                    if (viz[other] == 0) {
                        viz[other] = 1;
                        s.push(other);
                    }
                }
            }
        }
        return 0;
    }

    public static int DFC2(int src, Graph<Integer, ColoredEdge> graph, int currentColor, int targetColor, DisjointSets forest, HashMap<Integer, Pair<Integer, Integer>> colors) {
        int t = 0;
        Stack<Integer> s = new Stack();
        s.push(src);
        targetColor = forest.find(targetColor);
        while (!s.empty()) {
            src = s.pop();
            Pair<Integer, Integer> p = colors.get(src);
            Integer i1 = forest.find(p.getKey());
            Integer i2 = forest.find(p.getValue());
            if (i1 == targetColor) {
                colors.put(src, new Pair(currentColor, i2));
            } else if (i2 == targetColor) {
                colors.put(src, new Pair(i1, currentColor));
            }
            for (ColoredEdge e : graph.edgesOf(src)) {
                if (forest.find(e.color) == targetColor) {
                    int other = Graphs.getOppositeVertex(graph, e, src);
                    e.color = currentColor;

                    t++;
                    s.push(other);
                }
            }
        }
        return t;
    }

    public static int buildMaximalMatching2(Graph<Integer, ColoredEdge> graph, HashMap<Integer, Pair<Integer, Integer>> colors, DisjointSets forest) {
        int currentColor = 1;
        ArrayList<ColoredEdge> sortset = new ArrayList<ColoredEdge>();
        sortset.addAll(graph.edgeSet());

        Collections.sort(sortset, new EdgeDegComparator(graph));
        Graph<Integer, ColoredEdge> m = new SimpleGraph<>(ColoredEdge.class);

        HashMap<Integer, Integer> hm = new HashMap();

        for (ColoredEdge e : sortset) {
            Integer src = graph.getEdgeSource(e);
            Integer tgt = graph.getEdgeTarget(e);

            if (!hm.containsKey(src) && !hm.containsKey(tgt)) {
                hm.put(src, currentColor);
                hm.put(tgt, currentColor);
                e.color = currentColor;
                currentColor++;
            }
        }

        for (Integer v : graph.vertexSet()) {
            if (hm.containsKey(v)) {
                int r = DFC(v, graph, currentColor, 0);
                if (r != 0) {
                    colors.put(v, new Pair(hm.get(v), currentColor));
                    forest.number[currentColor] = r;
                    if (r > forest.maxnumber) {
                        forest.maxnumber = r;
                    }
                    currentColor++;
                } else {
                    int color1 = hm.get(v);
                    int color2 = 0;
                    for (ColoredEdge e : graph.edgesOf(v)) {
                        if (e.color != color1) {
                            color2 = e.color;
                            break;
                        }
                    }
                    colors.put(v, new Pair(color1, color2));
                }
            } else {
                int color2 = 0;
                for (ColoredEdge e : graph.edgesOf(v)) {
                    color2 = e.color;
                    break;
                }
                colors.put(v, new Pair(color2, 0));
            }

        }
        //System.out.println();
        //System.out.println(colors);

        /*
                        HashSet<ColoredEdge> hs = new HashSet<>();

        ConnectivityInspector ci = new ConnectivityInspector(graph);
                Set<Integer> s = ci.connectedSetOf(src);
                for (Integer n : s) {
                    hs.addAll(graph.edgesOf(n));
                }
         */
        //System.out.println();
        //System.out.print("\t" + (hm.keySet().size() / 2));
        //System.out.println(graph.edgeSet());
        return currentColor;
    }//end initial block

    public static int testMinMaxMovesetAnnealing(
            boolean useMatching, int cutoff,
            float temp, float tempfactor, float b1, float b2, float b3, float b4, float w1, float w2, float w3, float w4,
            Graph<Integer, ColoredEdge> graph, int NN, boolean displaySteps, boolean display, int N) {
        /*
        builds a maximal matching and uses it as an initial solution;
        then: sorts the edges after a heuristic attractiveness of edges;
        the attractiveness was the attempt to select edges that are likely
        to disconnect the uncolored graph: but this is actually hard you need
        something like separators to do reliably
        after that: cycle through the sorted edges and find the one that disconnects
        the underlying connected component into balanced pieces
         */

        long time = System.currentTimeMillis();
        long timebest =0;
        int M = -1;
        for (int i : graph.vertexSet()) {
            if (i > M) {
                M = i;
            }
        }
        HashMap<Integer, Pair<Integer, Integer>> colors = new HashMap<>(2 * M);

        Graph<Integer, ColoredEdge> result = new SimpleGraph<>(ColoredEdge.class);
        Graphs.addAllVertices(result, graph.vertexSet());

        DisjointSets forest = new DisjointSets(graph.edgeSet().size() + 5001);
        //build initial solution
        int currentColor = 0;
        if (useMatching) {
            for (ColoredEdge e : graph.edgeSet()) {
                e.color = 0;
            }
            currentColor = buildMaximalMatching2(graph, colors, forest);//start from a maximal matching
        } else {
            for (ColoredEdge e : graph.edgeSet()) {
                e.color = 1;
            }
            for (Integer v : graph.vertexSet()) {
                colors.put(v, new Pair(1, 0));
            }

            forest.maxnumber = graph.edgeSet().size();
            forest.number[1] = graph.edgeSet().size();
            currentColor = 2;
        }

        //System.out.print(forest);
        Random rng = new Random(111);

        //float temp = 0.2f;
        float best = 0.0f;
        int c = 0;
        ColoredEdge bestEdge = null;
        //System.out.println();
        int test = 0;
        int cleanupCounter = 0;
        int timeCounter = 0;

        ColoredEdge[] edgeList = new ColoredEdge[graph.edgeSet().size()];
        {
            int t = 0;
            for (ColoredEdge e : graph.edgeSet()) {
                edgeList[t++] = e;
            }
        }
        int[] indices = new int[graph.edgeSet().size()];
        {
            int t = 0;
            for (t = 0; t < graph.edgeSet().size(); t++) {
                indices[t] = t;
            }
        }
        int bestnumber = graph.edgeSet().size();
        int iterationbest = 0;

        do {
            timeCounter++;
            {//Knuth Permutation of edges for better SA
                int n = graph.edgeSet().size();
                int i, j;
                for (i = 0; i <= n - 2; i++) {
                    j = r.nextInt(n - i);
                    int aux = indices[i];
                    indices[i] = indices[i + j];
                    indices[i + j] = aux;
                }
            }

            best = 0.0f;
            c = 0;
            bestEdge = null;

            for (int edgeIndex : indices/*ColoredEdge e : graph.edgeSet()*/) {
                ColoredEdge e = edgeList[edgeIndex];
                int i = graph.getEdgeSource(e);
                int j = graph.getEdgeTarget(e);
                int c1, c2, c3, c4;
                c1 = forest.find(colors.get(i).getKey());
                c2 = forest.find(colors.get(i).getValue());
                c3 = forest.find(colors.get(j).getKey());
                c4 = forest.find(colors.get(j).getValue());
                int ecolor = forest.find(e.color);
                int col1 = (ecolor == c1) ? c2 : c1;
                int col2 = (ecolor == c3) ? c4 : c3;

                float b = 0.0f;

                if (col1 == 0 && col2 == 0) {
                    //both fully same color
                    //can add a newcolor on this edge
                    b = b1 + w1 * ((float) (forest.number[ecolor]) / forest.maxnumber);

                    if (b > best || rng.nextFloat() * 1000 < Math.exp((b - best) / temp)) {
                        bestEdge = e;
                        best = b;
                        c = 1;
                    }
                } else if ((col1 == 0) || (col2 == 0) || (col1 == col2)) {
                    //first node is only 1 color (=ecolor)
                    //or second node is only 1 color
                    //or the first node and second both have ecolor and col1
                    //edge can take the color that is not ecolor
                    int otherColor = (col1 == 0) ? col2 : col1;
                    b = b2 + w2 * ((float) forest.number[ecolor] * ((float) (forest.maxnumber - forest.number[otherColor]) / forest.maxnumber));

                    //^ that value should be 0 iff otherColor is actually the largest
                    if (b > best || rng.nextFloat() < Math.exp((b - best) / temp)) {
                        bestEdge = e;
                        best = b;
                        c = 2;
                    }

                    b = b4 + w4 * ((float) (forest.maxnumber - forest.number[ecolor] - forest.number[otherColor]) / forest.maxnumber);
                    //^ this value should be 0 or less when im trying to bind things that exceed maxnumber
                    if (b > best | rng.nextFloat() * 10000 < Math.exp((b - best) / temp)) {
                        //the move effect is too powerful so it is called randomly less often
                        bestEdge = e;
                        best = b;
                        c = 4;
                    }

                } else {
                    //nodes got 3 different colors
                    //can try to bind the different ones
                    b = b3 + w3 * ((float) forest.number[ecolor] * ((float) (forest.maxnumber - forest.number[col1] - forest.number[col2] - 1) / forest.maxnumber));

                    //^ this value should be 0 or less when im trying to bind things that exceed maxnumber
                    if (b > best || rng.nextFloat() * 100 < Math.exp((b - best) / temp)) {
                        //the move effect is too powerful so it is called randomly less often
                        bestEdge = e;
                        best = b;
                        c = 3;
                    }

                }

            }

            if (bestEdge != null) {
                //System.out.print(best);
                int i = graph.getEdgeSource(bestEdge);
                int j = graph.getEdgeTarget(bestEdge);
                int c1, c2, c3, c4;
                c1 = forest.find(colors.get(i).getKey());
                c2 = forest.find(colors.get(i).getValue());
                c3 = forest.find(colors.get(j).getKey());
                c4 = forest.find(colors.get(j).getValue());
                //System.out.println(" ("+c1+" "+c2+") ("+c3+" "+c4+")");
                int bcolor = forest.find(bestEdge.color);
                int col1 = (bcolor == c1) ? c2 : c1;
                int col2 = (bcolor == c3) ? c4 : c3;
                if (c == 1) {

                    forest.number[bcolor]--;
                    bestEdge.color = currentColor;

                    if (graph.degreeOf(i) == 1) {
                        colors.put(i, new Pair(currentColor, 0));
                    } else {
                        colors.put(i, new Pair(bcolor, currentColor));
                    }
                    if (graph.degreeOf(j) == 1) {
                        colors.put(j, new Pair(currentColor, 0));
                    } else {
                        colors.put(j, new Pair(bcolor, currentColor));
                    }

                    currentColor++;
                }
                if (c == 2) {
                    //System.out.println("case2 [" + i + " " + j);
                    int otherColor = (col1 == 0) ? col2 : col1;

                    forest.number[bcolor]--;
                    forest.number[otherColor]++;

                    bestEdge.color = otherColor;

                    int cnt = 0;
                    for (ColoredEdge ee : graph.edgesOf(i)) {
                        if (forest.find(ee.color) == bcolor) {
                            cnt++;
                            break;
                        }
                    }
                    if (cnt == 0) {
                        colors.put(i, new Pair(otherColor, 0));
                    } else {
                        colors.put(i, new Pair(bcolor, otherColor));
                    }
                    cnt = 0;
                    for (ColoredEdge ee : graph.edgesOf(j)) {
                        if (forest.find(ee.color) == bcolor) {
                            cnt++;
                            break;
                        }
                    }
                    if (cnt == 0) {
                        colors.put(j, new Pair(otherColor, 0));
                    } else {
                        colors.put(j, new Pair(bcolor, otherColor));
                    }

                } else if (c == 3) {
                    //System.out.println("case3 [" + i + " " + j);
                    forest.number[bcolor]--;
                    bestEdge.color = forest.union(col1, col2);

                    int cnt = 0;
                    for (ColoredEdge ee : graph.edgesOf(i)) {
                        if (forest.find(ee.color) == bcolor) {
                            cnt++;
                            break;
                        }
                    }
                    if (cnt == 0) {
                        colors.put(i, new Pair(bestEdge.color, 0));
                    } else {
                        colors.put(i, new Pair(bcolor, bestEdge.color));
                    }
                    cnt = 0;
                    for (ColoredEdge ee : graph.edgesOf(j)) {
                        if (forest.find(ee.color) == bcolor) {
                            cnt++;
                            break;
                        }
                    }
                    if (cnt == 0) {
                        colors.put(j, new Pair(bestEdge.color, 0));
                    } else {
                        colors.put(j, new Pair(bcolor, bestEdge.color));
                    }

                } else if (c == 4) {
                    bestEdge.color = forest.union(bcolor, col1);
                    forest.number[bestEdge.color]--;
                    for (Integer vt : graph.vertexSet()) {
                        Pair<Integer, Integer> p = colors.get(vt);
                        if (forest.find(p.getKey()) == forest.find(p.getValue())) {
                            colors.put(vt, new Pair(forest.find(p.getKey()), 0));
                        }
                    }
                }

                if (displaySteps) {
                    for (ColoredEdge e : graph.edgeSet()) {
                        e.color = forest.find(e.color);
                    }
                    displayGraph(graph, N, false);
                }

                //System.out.println("DFS");
                int dfsresult = DFS(i, j, graph, bcolor, forest, NN);
                if (dfsresult == 0) {
                    //System.out.println("ZERO"+bcolor);
                    forest.number[currentColor] = DFC2(i, graph, currentColor, bcolor, forest, colors);
                    forest.number[bcolor] -= forest.number[currentColor];
                    currentColor++;
                } else if (dfsresult == -1) {
                    //quite messy, other cases of recoloring are handled in DFC2
                    //i need to remove bcolor from the source vtx
                    Pair<Integer, Integer> values = colors.get(i);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(i, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(i, new Pair(values.getKey(), 0));
                    }
                } else if (dfsresult == -2) {
                    //i need to remove bcolor from the target vtx
                    Pair<Integer, Integer> values = colors.get(j);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(j, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(j, new Pair(values.getKey(), 0));
                    }
                } else if (dfsresult == -3) {
                    //i need to remove bcolor from the target vtx
                    Pair<Integer, Integer> values = colors.get(i);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(i, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(i, new Pair(values.getKey(), 0));
                    }
                    values = colors.get(j);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(j, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(j, new Pair(values.getKey(), 0));
                    }
                }

                int iii = 0;
                int maxxx = 0;
                for (int ii = 1; ii < currentColor; ii++) {
                    int iiiiii = forest.find(ii);
                    if (forest.number[iiiiii] > maxxx) {
                        maxxx = forest.number[iiiiii];
                        iii = iiiiii;
                    }
                    forest.maxnumber = maxxx;
                }
                //System.out.println("\nmaxnumber:" + forest.maxnumber + " at color " + iii);

            }

            temp *= tempfactor;
            /*
            for (Integer vv : graph.vertexSet()) {
                int c11 = 0, c22 = 0, c33 = 0;
                for (ColoredEdge e : graph.edgesOf(vv)) {
                    if (c11 == 0) {
                        c11 = forest.find(e.color);
                    } else if (c22 == 0 && forest.find(e.color) != c11) {
                        c22 = forest.find(e.color);
                    } else if (forest.find(e.color) != c11 && forest.find(e.color) != c22) {
                        System.out.println("!!!!!!!!" + vv + " " + forest.find(e.color));
                        for (ColoredEdge ee : graph.edgeSet()) {
                            ee.color = forest.find(ee.color);
                        }

                        //displayGraph(graph, N);
                        System.exit(0);
                    }
                }
            }*/

            if (temp < 5f) {
                test++;
            }
            if (forest.maxnumber < bestnumber) {
                test = 0;
                bestnumber = forest.maxnumber;
                iterationbest = timeCounter;
                timebest=System.currentTimeMillis();
            }

            if (temp < 1f && test > 10000) {
                break;
            }

            /*
            for (Integer vtx:graph.vertexSet())
            {Pair<Integer,Integer> p =colors.get(vtx);
                int aaa=forest.find(p.getKey());
                int bbb=forest.find(p.getValue());
                boolean one=false,two=false;
                if (aaa==0) one=true;
                if (bbb==0) two=true;
                for (ColoredEdge e:graph.edgesOf(vtx))
                    if (forest.find(e.color)==aaa) one =true;
                    else if (forest.find(e.color)==bbb) two =true;
                
                if (one&&two)continue;
                System.out.println("major fuckup "+vtx+ " "+aaa+" "+bbb+ " "+cleanupCounter);
                
                for (ColoredEdge e:graph.edgesOf(vtx))
                    System.out.println(forest.find(e.color));
                
                
                                for (ColoredEdge e : graph.edgeSet()) {
                    e.color = forest.find(e.color);
                }
                
                displayGraph(graph, N);
                
                System.exit(15);
                
            }*/
            cleanupCounter++;

            if (cleanupCounter == 2000) {
                cleanupCounter = 0;
                HashSet<Integer> set = new HashSet<>();
                int idx = 1;
                int nrs[] = new int[M + 1];
                HashMap<Integer, Integer> mapping = new HashMap<>();

                for (ColoredEdge e : graph.edgeSet()) {
                    int ecol = forest.find(e.color);
                    if (!mapping.containsKey(ecol)) {
                        mapping.put(ecol, idx);
                        e.color = idx;
                        idx++;
                    } else {
                        e.color = mapping.get(ecol);
                    }
                }

                currentColor = idx;

                for (Integer vtx : graph.vertexSet()) {
                    Pair<Integer, Integer> p = colors.get(vtx);
                    int aaa = forest.find(p.getKey());
                    //System.out.println("aaa" + aaa + " " + mapping.containsKey(aaa) + " " + vtx);

                    if (aaa != 0) {
                        aaa = mapping.get(aaa);
                    }
                    int bbb = forest.find(p.getValue());
                    if (bbb != 0) {
                        bbb = mapping.get(bbb);
                    }
                    colors.put(vtx, new Pair(aaa, bbb));
                }
                for (Integer elem : mapping.keySet()) {
                    nrs[mapping.get(elem)] = forest.number[elem];
                }
                forest.makeSet();
                for (idx = 1; idx < currentColor; idx++) {
                    //System.out.print(nrs[idx] + " ");

                    forest.number[idx] = nrs[idx];
                }

                //displayGraph(graph, N);
            }
            
            /*if ((timeCounter % 10) == 1) 
                System.out.println(timeCounter + "\t"+ forest.maxnumber);
            */
            if (cutoff != 0 && timeCounter == cutoff) {
                //System.out.print("iteration: " + iterationbest + " best: " + bestnumber + " last: " + forest.maxnumber);
                System.out.print(iterationbest + "\t" + bestnumber );
                printTime(time,timebest);
                System.out.println();
                if (display) {

                    for (ColoredEdge e : graph.edgeSet()) {
                        e.color = forest.find(e.color);
                    }
                    displayGraph(graph, N, true);
                }

                return bestnumber;
            }

            if (displaySteps) {
                for (ColoredEdge e : graph.edgeSet()) {
                    e.color = forest.find(e.color);
                }
                displayGraph(graph, N, false);
            }
        } while (bestEdge != null);

        //System.out.print("iteration: " + iterationbest + " best: " + bestnumber + " last: " + forest.maxnumber);
        
        System.out.print(iterationbest + "\t" + bestnumber );
        printTime(time,timebest);
        System.out.println();
        //printMinMaxStatsMatchColoring(time, graph, result, forest);
        if (display) {
            for (ColoredEdge e : graph.edgeSet()) {
                e.color = forest.find(e.color);
            }
            displayGraph(graph, N, true);
        }

        return bestnumber;

    }

    public static int testMinMaxMovesetTabu(
            boolean useMatching, int cutoff, int tabuTenure, float freqPenalty,
            float b1, float b2, float b3, float b4, float w1, float w2, float w3, float w4,
            Graph<Integer, ColoredEdge> graph, int NN, boolean displaySteps, boolean display, int N) {
        /*
        builds a maximal matching and uses it as an initial solution;
        then: sorts the edges after a heuristic attractiveness of edges;
        the attractiveness was the attempt to select edges that are likely
        to disconnect the uncolored graph: but this is actually hard you need
        something like separators to do reliably
        after that: cycle through the sorted edges and find the one that disconnects
        the underlying connected component into balanced pieces
         */

        long time = System.currentTimeMillis();
        long timebest=0;
        int M = -1;
        for (int i : graph.vertexSet()) {
            if (i > M) {
                M = i;
            }
        }
        HashMap<Integer, Pair<Integer, Integer>> colors = new HashMap<>(2 * M);

        Graph<Integer, ColoredEdge> result = new SimpleGraph<>(ColoredEdge.class);
        Graphs.addAllVertices(result, graph.vertexSet());

        DisjointSets forest = new DisjointSets(graph.edgeSet().size() + 5001);
        //build initial solution
        int currentColor = 0;
        if (useMatching) {
            for (ColoredEdge e : graph.edgeSet()) {
                e.color = 0;
            }
            currentColor = buildMaximalMatching2(graph, colors, forest);//start from a maximal matching
        } else {
            for (ColoredEdge e : graph.edgeSet()) {
                e.color = 1;
            }
            for (Integer v : graph.vertexSet()) {
                colors.put(v, new Pair(1, 0));
            }

            forest.maxnumber = graph.edgeSet().size();
            forest.number[1] = graph.edgeSet().size();
            currentColor = 2;
        }

        //System.out.print(forest);
        HashMap<ColoredEdge, Integer> tabuList2 = new HashMap<>();

        HashMap<ColoredEdge, Integer> edgeFrequency = new HashMap<>(graph.edgeSet().size() * 2);
        for (ColoredEdge e : graph.edgeSet()) {
            edgeFrequency.put(e, 0);
        }

        float best = 0.0f;
        int c = 0;
        ColoredEdge bestEdge = null;
        //System.out.println();
        int test = 0;
        int bestnumber = graph.edgeSet().size();
        int iterationbest = 0;
        int cleanupCounter = 0;
        int timeCounter = 0;
        do {
            timeCounter++;
            best = 0.0f;
            c = 0;
            bestEdge = null;

            for (ColoredEdge e : graph.edgeSet()) {
                if (tabuList2.containsKey(e)) {
                    continue;
                }
                int i = graph.getEdgeSource(e);
                int j = graph.getEdgeTarget(e);
                int c1, c2, c3, c4;
                c1 = forest.find(colors.get(i).getKey());
                c2 = forest.find(colors.get(i).getValue());
                c3 = forest.find(colors.get(j).getKey());
                c4 = forest.find(colors.get(j).getValue());
                int ecolor = forest.find(e.color);
                int col1 = (ecolor == c1) ? c2 : c1;
                int col2 = (ecolor == c3) ? c4 : c3;

                float b = 0.0f;

                if (col1 == 0 && col2 == 0) {
                    //both fully same color
                    //can add a newcolor on this edge
                    b = b1 + w1 * ((float) (forest.number[ecolor]) / forest.maxnumber);
                    b -= freqPenalty * edgeFrequency.get(e);
                    //TODO TODO
                    if (b > best) {
                        bestEdge = e;
                        best = b;
                        c = 1;

                    }
                } else if ((col1 == 0) || (col2 == 0) || (col1 == col2)) {
                    //first node is only 1 color (=ecolor)
                    //or second node is only 1 color
                    //or the first node and second both have ecolor and col1
                    //edge can take the color that is not ecolor
                    int otherColor = (col1 == 0) ? col2 : col1;
                    b = b2 + w3 * ((float) forest.number[ecolor] * ((float) (forest.maxnumber - forest.number[otherColor]) / forest.maxnumber));
                    b -= freqPenalty * edgeFrequency.get(e);

                    //^ that value should be 0 iff otherColor is actually the largest
                    if (b > best) {
                        bestEdge = e;
                        best = b;
                        c = 2;
                    }

                    b = b4 + w4 * ((float) (forest.maxnumber - forest.number[ecolor] - forest.number[otherColor]) / forest.maxnumber);

                    //^ this value should be 0 or less when im trying to bind things that exceed maxnumber
                    if (b > best) {
                        bestEdge = e;
                        best = b;
                        c = 4;
                    }
                } else {
                    //nodes got 3 different colors
                    //can try to bind the different ones

                    b = b3 + w3 * ((float) forest.number[ecolor] * ((float) (forest.maxnumber - forest.number[col1] - forest.number[col2] - 1) / forest.maxnumber));
                    b -= freqPenalty * edgeFrequency.get(e);

                    //^ this value should be 0 or less when im trying to bind things that exceed maxnumber
                    if (b > best) {
                        bestEdge = e;
                        best = b;
                        c = 3;
                    }

                }

            }

            if (bestEdge != null) {
                edgeFrequency.put(bestEdge, edgeFrequency.get(bestEdge) + 1);
                //System.out.print(best);
                int i = graph.getEdgeSource(bestEdge);
                int j = graph.getEdgeTarget(bestEdge);
                int c1, c2, c3, c4;
                c1 = forest.find(colors.get(i).getKey());
                c2 = forest.find(colors.get(i).getValue());
                c3 = forest.find(colors.get(j).getKey());
                c4 = forest.find(colors.get(j).getValue());
                //System.out.println(" ("+c1+" "+c2+") ("+c3+" "+c4+")");
                int bcolor = forest.find(bestEdge.color);
                int col1 = (bcolor == c1) ? c2 : c1;
                int col2 = (bcolor == c3) ? c4 : c3;
                tabuList2.put(bestEdge, timeCounter);
                /*
                if (colorFrequency.containsKey(bcolor))
                    colorFrequency.put(bcolor, colorFrequency.get(bcolor)+1);
                else
                    colorFrequency.put(bcolor, 1);
                 */
                if (c == 1) {
                    {
                        //System.out.println("case1 [" + i + " " + j);
                        forest.number[bcolor]--;
                        bestEdge.color = currentColor;

                        if (graph.degreeOf(i) == 1) {
                            colors.put(i, new Pair(currentColor, 0));
                        } else {
                            colors.put(i, new Pair(bcolor, currentColor));
                        }
                        if (graph.degreeOf(j) == 1) {
                            colors.put(j, new Pair(currentColor, 0));
                        } else {
                            colors.put(j, new Pair(bcolor, currentColor));
                        }

                        currentColor++;
                    }
                } else if (c == 2) {

                    {
                        //System.out.println("case2 [" + i + " " + j);
                        int otherColor = (col1 == 0) ? col2 : col1;

                        forest.number[bcolor]--;
                        forest.number[otherColor]++;

                        bestEdge.color = otherColor;

                        int cnt = 0;
                        for (ColoredEdge ee : graph.edgesOf(i)) {
                            if (forest.find(ee.color) == bcolor) {
                                cnt++;
                            }
                        }
                        if (cnt == 0) {
                            colors.put(i, new Pair(otherColor, 0));
                        } else {
                            colors.put(i, new Pair(bcolor, otherColor));
                        }
                        cnt = 0;
                        for (ColoredEdge ee : graph.edgesOf(j)) {
                            if (forest.find(ee.color) == bcolor) {
                                cnt++;
                            }
                        }
                        if (cnt == 0) {
                            colors.put(j, new Pair(otherColor, 0));
                        } else {
                            colors.put(j, new Pair(bcolor, otherColor));
                        }
                        //colorFrequency.put(otherColor,colorFrequency.containsKey(otherColor)?colorFrequency.get(otherColor)+1:1);
                    }

                } else if (c == 3) {

                    {
                        //System.out.println(timeCounter+" case3 ");
                        //System.out.println("case3 [" + i + " " + j);
                        forest.number[bcolor]--;
                        bestEdge.color = forest.union(col1, col2);
                        /*
                        int v=0;
                        v+=colorFrequency.containsKey(col1)?colorFrequency.get(col1):0;
                        v+=colorFrequency.containsKey(col2)?colorFrequency.get(col2):0;
                        
                        colorFrequency.remove(col1);
                        colorFrequency.remove(col2);
                        if (v!=0)
                            colorFrequency.put(bestEdge.color,v);
                         */
                        
                        /*int cnt = 0;
                        for (ColoredEdge ee : graph.edgesOf(i)) {
                            if (forest.find(ee.color) == bcolor) {
                                cnt++;
                                break;
                            }
                        }
                        if (cnt == 0 || bcolor==bestEdge.color) {
                            colors.put(i, new Pair(bestEdge.color, 0));
                        } else {
                            colors.put(i, new Pair(bcolor, bestEdge.color));
                        }
                        cnt = 0;
                        for (ColoredEdge ee : graph.edgesOf(j)) {
                            if (forest.find(ee.color) == bcolor) {
                                cnt++;
                                break;
                            }
                        }
                        if (cnt == 0 || bcolor==bestEdge.color) {
                            colors.put(j, new Pair(bestEdge.color, 0));
                        } else {
                            colors.put(j, new Pair(bcolor, bestEdge.color));
                        }*/

                    }

                } else if (c == 4) {
                    //System.out.println("case4 [" + i + " " + j);
                    //System.out.println(timeCounter+" case4 ");
                    int otherColor = (col1 == 0) ? col2 : col1;
                    bestEdge.color = forest.union(bcolor, otherColor);
                    forest.number[bestEdge.color]--;
                    for (Integer vt : graph.vertexSet()) {
                        Pair<Integer, Integer> p = colors.get(vt);
                        if (forest.find(p.getKey()) == forest.find(p.getValue())) {
                            colors.put(vt, new Pair(forest.find(p.getKey()), 0));
                        }
                    }
                }

                if (displaySteps) {
                    for (ColoredEdge e : graph.edgeSet()) {
                        e.color = forest.find(e.color);
                    }
                    displayGraph(graph, N, false);
                }

                //System.out.println("DFS");
                int dfsresult = DFS(i, j, graph, bcolor, forest, NN);
                if (dfsresult == 0) {
                    //System.out.println("ZERO"+bcolor);
                    forest.number[currentColor] = DFC2(i, graph, currentColor, bcolor, forest, colors);
                    forest.number[bcolor] -= forest.number[currentColor];
                    currentColor++;
                } else if (dfsresult == -1) {
                    //quite messy, other cases of recoloring are handled in DFC2
                    //i need to remove bcolor from the source vtx
                    Pair<Integer, Integer> values = colors.get(i);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(i, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(i, new Pair(values.getKey(), 0));
                    }
                } else if (dfsresult == -2) {
                    //i need to remove bcolor from the target vtx
                    Pair<Integer, Integer> values = colors.get(j);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(j, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(j, new Pair(values.getKey(), 0));
                    }
                } else if (dfsresult == -3) {
                    //i need to remove bcolor from the target vtx
                    Pair<Integer, Integer> values = colors.get(i);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(i, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(i, new Pair(values.getKey(), 0));
                    }
                    values = colors.get(j);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(j, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(j, new Pair(values.getKey(), 0));
                    }
                }

                int iii = 0;
                int maxxx = 0;
                for (int ii = 1; ii < currentColor; ii++) {
                    int iiiiii = forest.find(ii);
                    if (forest.number[iiiiii] > maxxx) {
                        maxxx = forest.number[iiiiii];
                        iii = iiiiii;
                    }
                    forest.maxnumber = maxxx;
                }
                //System.out.println("\nmaxnumber:" + forest.maxnumber + " at color " + iii);
            }

            /*
            for (Integer vv : graph.vertexSet()) {
                int c11 = 0, c22 = 0, c33 = 0;
                for (ColoredEdge e : graph.edgesOf(vv)) {
                    if (c11 == 0) {
                        c11 = forest.find(e.color);
                    } else if (c22 == 0 && forest.find(e.color) != c11) {
                        c22 = forest.find(e.color);
                    } else if (forest.find(e.color) != c11 && forest.find(e.color) != c22) {
                        System.out.println("!!!!!!!!" + vv + " " + forest.find(e.color));
                        for (ColoredEdge ee : graph.edgeSet()) {
                            ee.color = forest.find(ee.color);
                        }

                        displayGraph(graph, N,false);
                        System.exit(0);
                    }
                }
            }*/
            /*if ((timeCounter % 10) == 1) 
                System.out.println(timeCounter + "\t"+ forest.maxnumber);
            */
            if (forest.maxnumber < bestnumber) {
                test = 0;
                bestnumber = forest.maxnumber;
                iterationbest = timeCounter;
                timebest=System.currentTimeMillis();
            }
            test++;
            if (test > 10000) {

                break;

            }

            /*
            for (Integer vtx:graph.vertexSet())
            {Pair<Integer,Integer> p =colors.get(vtx);
                int aaa=forest.find(p.getKey());
                int bbb=forest.find(p.getValue());
                boolean one=false,two=false;
                if (aaa==0) one=true;
                if (bbb==0) two=true;
                for (ColoredEdge e:graph.edgesOf(vtx))
                    if (forest.find(e.color)==aaa) one =true;
                    else if (forest.find(e.color)==bbb) two =true;
                
                if (one&&two)continue;
                System.out.println("major fuckup "+vtx+ " "+aaa+" "+bbb+ " "+cleanupCounter);
                
                for (ColoredEdge e:graph.edgesOf(vtx))
                    System.out.println(forest.find(e.color));
                
                
                                for (ColoredEdge e : graph.edgeSet()) {
                    e.color = forest.find(e.color);
                }
                
                displayGraph(graph, N,false);
                
                System.exit(15);
                
            }*/
            cleanupCounter++;
            if (cleanupCounter == 2000) {
                cleanupCounter = 0;
                HashSet<Integer> set = new HashSet<>();
                int idx = 1;
                int nrs[] = new int[M + 1];
                HashMap<Integer, Integer> mapping = new HashMap<>();

                for (ColoredEdge e : graph.edgeSet()) {
                    int ecol = forest.find(e.color);
                    if (!mapping.containsKey(ecol)) {
                        mapping.put(ecol, idx);
                        e.color = idx;
                        idx++;
                    } else {
                        e.color = mapping.get(ecol);
                    }
                }

                currentColor = idx;

                for (Integer vtx : graph.vertexSet()) {
                    Pair<Integer, Integer> p = colors.get(vtx);
                    int aaa = forest.find(p.getKey());
                    //System.out.println("aaa" + aaa + " " + mapping.containsKey(aaa) + " " + vtx);

                    if (aaa != 0) {
                        aaa = mapping.get(aaa);
                    }

                    int bbb = forest.find(p.getValue());
                    //System.out.println("bbb" + bbb + " " + mapping.containsKey(bbb)) ;
                    if (bbb != 0) {
                        bbb = mapping.get(bbb);
                    }
                    colors.put(vtx, new Pair(aaa, bbb));
                }
                for (Integer elem : mapping.keySet()) {
                    nrs[mapping.get(elem)] = forest.number[elem];
                }
                forest.makeSet();
                for (idx = 1; idx < currentColor; idx++) {
                    //System.out.print(nrs[idx] + " ");

                    forest.number[idx] = nrs[idx];
                }

                //displayGraph(graph, N);
            }

            if (cutoff != 0 && timeCounter == cutoff) {
                //System.out.print("iteration: " + iterationbest + " best: " + bestnumber + " last: " + forest.maxnumber);
                System.out.print(iterationbest + "\t" + bestnumber );
                printTime(time,timebest);
                System.out.println();
                if (display) {
                    for (ColoredEdge e : graph.edgeSet()) {
                        e.color = forest.find(e.color);
                    }
                    displayGraph(graph, N, true);
                }
                return bestnumber;
            }

            if (tabuTenure < 15 || (timeCounter & 7) == 0) {
                ArrayList<ColoredEdge> al = new ArrayList<>();
                for (ColoredEdge ce : tabuList2.keySet()) {
                    if (timeCounter - tabuList2.get(ce) > tabuTenure) {
                        al.add(ce);
                    }
                }
                for (ColoredEdge ce : al) {
                    tabuList2.remove(ce);
                }

            }

            if (displaySteps) {
                for (ColoredEdge e : graph.edgeSet()) {
                    e.color = forest.find(e.color);
                }
                displayGraph(graph, N, true);
            }
        } while (bestEdge != null);

        //System.out.print("iteration: " + iterationbest + " best: " + bestnumber + " last: " + forest.maxnumber);
        System.out.print(iterationbest + "\t" + bestnumber );
        printTime(time,timebest);
        System.out.println();

        if (display) {
            for (ColoredEdge e : graph.edgeSet()) {
                e.color = forest.find(e.color);
            }
            displayGraph(graph, N, true);
        }
        return bestnumber;

    }

    public static int countmax(ArrayList<ColoredEdge> edges) {
        HashMap<Integer, Integer> count = new HashMap();

        int max = 0;
        for (ColoredEdge e : edges) {
            if (!count.containsKey(e.color)) {
                count.put(e.color, 1);
            } else {
                int val = count.get(e.color) + 1;
                count.put(e.color, val);
                if (max < val) {
                    max = val;
                }
            }
        }
        return max;
    }

    //bad recursive solver does not explore entire space
    static int mincol = 99999;
    static int ccol = 2;

    public static void solverec(int index, Graph<Integer, ColoredEdge> graph,
            ArrayList<ColoredEdge> edges, HashMap<Integer, TreeSet<Integer>> colors, int N) {

        if (index == N) {
            int res = countmax(edges);
            if (mincol > res) {
                mincol = res;
                System.out.println(mincol);

                displayGraph(graph, 5, false);
            }
            //if (res==17) displayGraph(graph,5,false);
            return;
        }

        int i = graph.getEdgeSource(edges.get(index));
        int j = graph.getEdgeTarget(edges.get(index));
        TreeSet<Integer> seti = colors.get(i);
        TreeSet<Integer> setj = colors.get(j);

        for (Integer el1 : seti) {
            for (Integer el2 : setj) {
                if (el1.equals(el2)) {
                    edges.get(index).color = el1;
                    solverec(index + 1, graph, edges, colors, N);
                }
            }
        }

        if (setj.size() == 0
                || (setj.size() == 1 && (!seti.contains(setj.first())))) {
            for (Integer el1 : seti) {
                setj.add(el1);
                edges.get(index).color = el1;
                solverec(index + 1, graph, edges, colors, N);
                setj.remove(el1);
            }
        }

        if (seti.size() == 0
                || (seti.size() == 1 && (!setj.contains(seti.first())))) {
            for (Integer el2 : setj) {
                seti.add(el2);
                edges.get(index).color = el2;
                solverec(index + 1, graph, edges, colors, N);
                seti.remove(el2);
            }
        }

        if (seti.size() < 2 && setj.size() < 2) {
            edges.get(index).color = ccol;
            seti.add(ccol);
            setj.add(ccol);
            ccol++;
            solverec(index + 1, graph, edges, colors, N);
            ccol--;
            seti.remove(ccol);
            setj.remove(ccol);
        }

    }

    public static void solve5x5(String filename) {

        Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
        int NN = readEdgesFromFile(graph, filename);

        int M = -1;
        for (int i : graph.vertexSet()) {
            if (i > M) {
                M = i;
            }
        }
        HashMap<Integer, TreeSet<Integer>> colors = new HashMap<>(2 * M);
        ArrayList<ColoredEdge> edges = new ArrayList<>();
        edges.addAll(graph.edgeSet());
        int N = edges.size();

        for (Integer vt : graph.vertexSet()) {
            colors.put(vt, new TreeSet());
        }

        solverec(0, graph, edges, colors, N);

    }

    public static int testMinMaxMovesetHillclimb(
            boolean useMatching, int cutoff,
            float b1, float b2, float b3, float b4, float w1, float w2, float w3, float w4,
            Graph<Integer, ColoredEdge> graph, int NN, boolean displaySteps, boolean display, int N) {
        /*
        builds a maximal matching and uses it as an initial solution;
        then: sorts the edges after a heuristic attractiveness of edges;
        the attractiveness was the attempt to select edges that are likely
        to disconnect the uncolored graph: but this is actually hard you need
        something like separators to do reliably
        after that: cycle through the sorted edges and find the one that disconnects
        the underlying connected component into balanced pieces
         */

        long time = System.currentTimeMillis();
        long timebest=0;
        int M = -1;
        for (int i : graph.vertexSet()) {
            if (i > M) {
                M = i;
            }
        }
        HashMap<Integer, Pair<Integer, Integer>> colors = new HashMap<>(2 * M);

        Graph<Integer, ColoredEdge> result = new SimpleGraph<>(ColoredEdge.class);
        Graphs.addAllVertices(result, graph.vertexSet());

        DisjointSets forest = new DisjointSets(graph.edgeSet().size() + 5001);
        //build initial solution
        int currentColor = 0;
        if (useMatching) {
            for (ColoredEdge e : graph.edgeSet()) {
                e.color = 0;
            }
            currentColor = buildMaximalMatching2(graph, colors, forest);//start from a maximal matching
        } else {
            for (ColoredEdge e : graph.edgeSet()) {
                e.color = 1;
            }
            for (Integer v : graph.vertexSet()) {
                colors.put(v, new Pair(1, 0));
            }

            forest.maxnumber = graph.edgeSet().size();
            forest.number[1] = graph.edgeSet().size();
            currentColor = 2;
        }
        float best = 0.0f;
        int c = 0;
        ColoredEdge bestEdge = null;
        //System.out.println();
        int test = 0;
        int bestnumber = forest.maxnumber;
        int cleanupCounter = 0;
        int timeCounter = 0;
        int iterationbest=0;
        do {
            timeCounter++;
            best = 0.0f;
            c = 0;
            bestEdge = null;

            for (ColoredEdge e : graph.edgeSet()) {

                int i = graph.getEdgeSource(e);
                int j = graph.getEdgeTarget(e);
                int c1, c2, c3, c4;
                c1 = forest.find(colors.get(i).getKey());
                c2 = forest.find(colors.get(i).getValue());
                c3 = forest.find(colors.get(j).getKey());
                c4 = forest.find(colors.get(j).getValue());
                int ecolor = forest.find(e.color);
                int col1 = (ecolor == c1) ? c2 : c1;
                int col2 = (ecolor == c3) ? c4 : c3;

                float b = 0.0f;

                if (col1 == 0 && col2 == 0) {
                    //both fully same color
                    //can add a newcolor on this edge
                    b = b1 + w1 * ((float) (forest.number[ecolor]) / forest.maxnumber);

                    //TODO TODO
                    if (b > best) {
                        bestEdge = e;
                        best = b;
                        c = 1;

                    }
                } else if ((col1 == 0) || (col2 == 0) || (col1 == col2)) {
                    //first node is only 1 color (=ecolor)
                    //or second node is only 1 color
                    //or the first node and second both have ecolor and col1
                    //edge can take the color that is not ecolor
                    int otherColor = (col1 == 0) ? col2 : col1;
                    b = b2 + w2 * ((float) forest.number[ecolor] * ((float) (forest.maxnumber - forest.number[otherColor]) / forest.maxnumber));

                    //^ that value should be 0 iff otherColor is actually the largest
                    if (b > best) {
                        bestEdge = e;
                        best = b;
                        c = 2;
                    }

                    b = b4 + w4 * ((float) (forest.maxnumber - forest.number[ecolor] - forest.number[otherColor]) / forest.maxnumber);

                    //^ this value should be 0 or less when im trying to bind things that exceed maxnumber
                    if (b > best) {
                        bestEdge = e;
                        best = b;
                        c = 4;
                    }
                } else {
                    //nodes got 3 different colors
                    //can try to bind the different ones
                    if ((forest.maxnumber - forest.number[col1] - forest.number[col2] - 1) <= 0) {
                        continue;
                    }
                    b = b3 + w3 * ((float) forest.number[ecolor] * ((float) (forest.maxnumber - forest.number[col1] - forest.number[col2] - 1) / forest.maxnumber));

                    //^ this value should be 0 or less when im trying to bind things that exceed maxnumber
                    if (b > best) {
                        bestEdge = e;
                        best = b;
                        c = 3;
                    }

                }

            }

            if (bestEdge != null) {

                //System.out.print(best);
                int i = graph.getEdgeSource(bestEdge);
                int j = graph.getEdgeTarget(bestEdge);
                int c1, c2, c3, c4;
                c1 = forest.find(colors.get(i).getKey());
                c2 = forest.find(colors.get(i).getValue());
                c3 = forest.find(colors.get(j).getKey());
                c4 = forest.find(colors.get(j).getValue());
                //System.out.println(" ("+c1+" "+c2+") ("+c3+" "+c4+")");
                int bcolor = forest.find(bestEdge.color);
                int col1 = (bcolor == c1) ? c2 : c1;
                int col2 = (bcolor == c3) ? c4 : c3;

                /*
                if (colorFrequency.containsKey(bcolor))
                    colorFrequency.put(bcolor, colorFrequency.get(bcolor)+1);
                else
                    colorFrequency.put(bcolor, 1);
                 */
                if (c == 1) {
                    {
                        //System.out.println("case1 [" + i + " " + j);
                        forest.number[bcolor]--;
                        bestEdge.color = currentColor;

                        if (graph.degreeOf(i) == 1) {
                            colors.put(i, new Pair(currentColor, 0));
                        } else {
                            colors.put(i, new Pair(bcolor, currentColor));
                        }
                        if (graph.degreeOf(j) == 1) {
                            colors.put(j, new Pair(currentColor, 0));
                        } else {
                            colors.put(j, new Pair(bcolor, currentColor));
                        }

                        currentColor++;
                    }
                } else if (c == 2) {

                    {
                        //System.out.println("case2 [" + i + " " + j);
                        int otherColor = (col1 == 0) ? col2 : col1;

                        forest.number[bcolor]--;
                        forest.number[otherColor]++;

                        bestEdge.color = otherColor;

                        int cnt = 0;
                        for (ColoredEdge ee : graph.edgesOf(i)) {
                            if (forest.find(ee.color) == bcolor) {
                                cnt++;
                            }
                        }
                        if (cnt == 0) {
                            colors.put(i, new Pair(otherColor, 0));
                        } else {
                            colors.put(i, new Pair(bcolor, otherColor));
                        }
                        cnt = 0;
                        for (ColoredEdge ee : graph.edgesOf(j)) {
                            if (forest.find(ee.color) == bcolor) {
                                cnt++;
                            }
                        }
                        if (cnt == 0) {
                            colors.put(j, new Pair(otherColor, 0));
                        } else {
                            colors.put(j, new Pair(bcolor, otherColor));
                        }

                    }

                } else if (c == 3) {

                    {
                        //System.out.println("case3 [" + i + " " + j);
                        forest.number[bcolor]--;
                        bestEdge.color = forest.union(col1, col2);

                        /*int cnt = 0;
                        for (ColoredEdge ee : graph.edgesOf(i)) {
                            if (forest.find(ee.color) == bcolor) {
                                cnt++;
                                break;
                            }
                        }
                        if (cnt == 0 || bcolor==bestEdge.color) {
                            colors.put(i, new Pair(bestEdge.color, 0));
                        } else {
                            colors.put(i, new Pair(bcolor, bestEdge.color));
                        }
                        cnt = 0;
                        for (ColoredEdge ee : graph.edgesOf(j)) {
                            if (forest.find(ee.color) == bcolor) {
                                cnt++;
                                break;
                            }
                        }
                        if (cnt == 0 || bcolor==bestEdge.color) {
                            colors.put(j, new Pair(bestEdge.color, 0));
                        } else {
                            colors.put(j, new Pair(bcolor, bestEdge.color));
                        }*/
                    }

                } else if (c == 4) {
                    //System.out.println(" case4 [" + i + " " + j);
                    int otherColor = (col1 == 0) ? col2 : col1;
                    bestEdge.color = forest.union(bcolor, otherColor);
                    forest.number[bestEdge.color]--;
                    for (Integer vt : graph.vertexSet()) {
                        Pair<Integer, Integer> p = colors.get(vt);
                        if (forest.find(p.getKey()) == forest.find(p.getValue())) {
                            colors.put(vt, new Pair(forest.find(p.getKey()), 0));
                        }
                    }
                }

                if (displaySteps) {
                    for (ColoredEdge e : graph.edgeSet()) {
                        e.color = forest.find(e.color);
                    }
                    displayGraph(graph, N, false);
                }

                //System.out.println("DFS");
                int dfsresult = DFS(i, j, graph, bcolor, forest, NN);
                if (dfsresult == 0) {
                    //System.out.println("ZERO"+bcolor);
                    forest.number[currentColor] = DFC2(i, graph, currentColor, bcolor, forest, colors);
                    forest.number[bcolor] -= forest.number[currentColor];
                    currentColor++;
                } else if (dfsresult == -1) {
                    //quite messy, other cases of recoloring are handled in DFC2
                    //i need to remove bcolor from the source vtx
                    Pair<Integer, Integer> values = colors.get(i);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(i, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(i, new Pair(values.getKey(), 0));
                    }
                } else if (dfsresult == -2) {
                    //i need to remove bcolor from the target vtx
                    Pair<Integer, Integer> values = colors.get(j);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(j, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(j, new Pair(values.getKey(), 0));
                    }
                } else if (dfsresult == -3) {
                    //i need to remove bcolor from the target vtx
                    Pair<Integer, Integer> values = colors.get(i);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(i, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(i, new Pair(values.getKey(), 0));
                    }
                    values = colors.get(j);
                    if (forest.find(values.getKey()) == bcolor) {
                        colors.put(j, new Pair(values.getValue(), 0));
                    } else {
                        colors.put(j, new Pair(values.getKey(), 0));
                    }
                }

                int iii = 0;
                int maxxx = 0;
                for (int ii = 1; ii < currentColor; ii++) {
                    int iiiiii = forest.find(ii);
                    if (forest.number[iiiiii] > maxxx) {
                        maxxx = forest.number[iiiiii];
                        iii = iiiiii;
                    }
                    forest.maxnumber = maxxx;
                }
                //System.out.println("\nmaxnumber:" + forest.maxnumber + " at color " + iii);
            }

            /*
            for (Integer vv : graph.vertexSet()) {
                int c11 = 0, c22 = 0, c33 = 0;
                for (ColoredEdge e : graph.edgesOf(vv)) {
                    if (c11 == 0) {
                        c11 = forest.find(e.color);
                    } else if (c22 == 0 && forest.find(e.color) != c11) {
                        c22 = forest.find(e.color);
                    } else if (forest.find(e.color) != c11 && forest.find(e.color) != c22) {
                        System.out.println("!!!!!!!!" + vv + " " + forest.find(e.color));
                        for (ColoredEdge ee : graph.edgeSet()) {
                            ee.color = forest.find(ee.color);
                        }

                        displayGraph(graph, N,false);
                        System.exit(0);
                    }
                }
            }*/
            if (forest.maxnumber < bestnumber) {
                test = 0;
                bestnumber = forest.maxnumber;
                iterationbest=timeCounter;
                timebest=System.currentTimeMillis();
                
            }
/*
            if ((timeCounter % 10) == 1) 
                System.out.println(timeCounter + "\t"+ forest.maxnumber);*/
            test++;
            if (test > 2000) {
                break;

            }

            if (cutoff != 0 && timeCounter == cutoff) {
                //System.out.print("iteration: " + iterationbest + " last: " + bestnumber);
                System.out.print(iterationbest + "\t" + bestnumber);
                printTime(time,timebest);
                System.out.println();
                if (display) {
                    for (ColoredEdge e : graph.edgeSet()) {
                        e.color = forest.find(e.color);
                    }
                    displayGraph(graph, N, true);
                }
                return forest.maxnumber;
            }

            /*
            for (Integer vtx:graph.vertexSet())
            {Pair<Integer,Integer> p =colors.get(vtx);
                int aaa=forest.find(p.getKey());
                int bbb=forest.find(p.getValue());
                boolean one=false,two=false;
                if (aaa==0) one=true;
                if (bbb==0) two=true;
                for (ColoredEdge e:graph.edgesOf(vtx))
                    if (forest.find(e.color)==aaa) one =true;
                    else if (forest.find(e.color)==bbb) two =true;
                
                if (one&&two)continue;
                System.out.println("major fuckup "+vtx+ " "+aaa+" "+bbb+ " "+cleanupCounter);
                
                for (ColoredEdge e:graph.edgesOf(vtx))
                    System.out.println(forest.find(e.color));
                
                
                                for (ColoredEdge e : graph.edgeSet()) {
                    e.color = forest.find(e.color);
                }
                
                displayGraph(graph, N,false);
                
                System.exit(15);
                
            }*/
            cleanupCounter++;
            if (cleanupCounter == 2000) {
                cleanupCounter = 0;
                HashSet<Integer> set = new HashSet<>();
                int idx = 1;
                int nrs[] = new int[M + 1];
                HashMap<Integer, Integer> mapping = new HashMap<>();

                for (ColoredEdge e : graph.edgeSet()) {
                    int ecol = forest.find(e.color);
                    if (!mapping.containsKey(ecol)) {
                        mapping.put(ecol, idx);
                        e.color = idx;
                        idx++;
                    } else {
                        e.color = mapping.get(ecol);
                    }
                }

                currentColor = idx;

                for (Integer vtx : graph.vertexSet()) {
                    Pair<Integer, Integer> p = colors.get(vtx);
                    int aaa = forest.find(p.getKey());
                    //System.out.println("aaa" + aaa + " " + mapping.containsKey(aaa) + " " + vtx);

                    if (aaa != 0) {
                        aaa = mapping.get(aaa);
                    }

                    int bbb = forest.find(p.getValue());
                    //System.out.println("bbb" + bbb + " " + mapping.containsKey(bbb)) ;
                    if (bbb != 0) {
                        bbb = mapping.get(bbb);
                    }
                    colors.put(vtx, new Pair(aaa, bbb));
                }
                for (Integer elem : mapping.keySet()) {
                    nrs[mapping.get(elem)] = forest.number[elem];
                }
                forest.makeSet();
                for (idx = 1; idx < currentColor; idx++) {
                    //System.out.print(nrs[idx] + " ");

                    forest.number[idx] = nrs[idx];
                }

                //displayGraph(graph, N);
            }

            if (displaySteps) {
                for (ColoredEdge e : graph.edgeSet()) {
                    e.color = forest.find(e.color);
                }
                displayGraph(graph, N, false);
            }
        } while (bestEdge != null);
        //System.out.print("iteration: " + iterationbest + " last: " + bestnumber);
        System.out.print(iterationbest + "\t" + bestnumber);
        printTime(time,timebest);
        System.out.println();
        //printMinMaxStatsMatchColoring(time, graph, result, forest);
        if (display) {

            for (ColoredEdge e : graph.edgeSet()) {
                e.color = forest.find(e.color);
            }
            displayGraph(graph, N, true);
        }

        return forest.maxnumber;

    }

    public static class record {

        int i;
        ArrayList<String> names;
        ArrayList<Float> entries;

        record(int ii) {
            i = ii;
            names = new ArrayList<>();
            entries = new ArrayList<>();
        }

        public void addEntry(String s, float f) {
            names.add(s);
            entries.add(f);
        }

        public void print() {
            DecimalFormat df = new DecimalFormat("#0.0");
            DecimalFormat df2 = new DecimalFormat("#0.00000");
            int j;
            System.out.print(i + ": ");
            for (j = 0; j < names.size(); j++) {
                System.out.print(names.get(j) + ":" + ((entries.get(j)<1f && entries.get(j)>-1f)?df2.format(entries.get(j)):df.format(entries.get(j))) + " ");
            }
            System.out.println();
        }
    };

    public static class recordComparator implements Comparator<record> {

        @Override
        public int compare(record t, record t1) {
            return Integer.compare(t.i, t1.i);
        }
    };

    public static void main(String[] args) {
        frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setPreferredSize(new Dimension(1200, 1000));

        String[] testfiles = {"107.edges.txt", "as19991205.txt", "as20000102.txt", "1912.edges.txt", /*"Email-Enron.txt",*/ "CA-CondMat.txt", "CA-AstroPh.txt", "erdosrenyi.txt", "Cit-HepTh.txt"};
        String[] testfiles2 = {"dsjc250.5", "dsjc500.1", "dsjc500.5", "dsjr500.5", "flat300_28_0", "le450_25c", "le450_25d", "r250.5"};

        System.out.println("graph name\tN\tM\tavgdeg\tmaxdeg\t|matching|\t|connect.comp w edges|\tcolors\tlargest\ttime");
        System.out.println();
        /*
       for (String s:testfiles)
           testGraphStandard(s);       
       System.out.println("-------------------------------------------------------------------------------------------------------------------------------");
       for (String s:testfiles)
           testGraphGreedy(s);
       System.out.println("-------------------------------------------------------------------------------------------------------------------------------");
         *//*for (String s:testfiles)
           testGraphGreedySort(s);
         */
        //System.out.println("graph name\tsize\tN\tM\t|connect.comp|\tconcompwith no edges\tcolors\t|largest.c|\tmaxdeg\ttime");
        //System.out.println("graph name\tsize\tN\tM\tmindeg\tmaxdeg\ttime");
        //System.out.println();
        /*for (String s : testfiles2) {
            testMinMaxMatchColoring(s, false, false, -1);
        }
         */
        //testMinMaxColoringClustering("smallworld.txt", false,false,-1);
        //testMinMaxMatchColoring("smallworld.txt", false,false,-1);
        String ss[] = {"udg500.1.txt", "udg500.2.txt", "udg500.3.txt", "udg1000.1.txt", "udg1000.2.txt", "udg1000.3.txt"};
        String sss[] = {"qudg500.1.txt", "qudg500.2.txt", "qudg500.3.txt", "qudg1000.1.txt", "qudg1000.2.txt", "qudg1000.3.txt"};
        String ss2[] = {"udg10000.1.txt", "udg10000.2.txt", "udg10000.3.txt"};
        String sss2[] = {"qudg10000.1.txt", "qudg10000.2.txt", "qudg10000.3.txt"};
        String[] testfiles22 = {"dsjc1000.5", "dsjc1000.9", "dsjc500.5", "dsjc500.9", "dsjr500.1c", "dsjr500.5", "flat1000_50_0", "flat1000_60_0", "flat1000_76_0", "r1000.1c", "r1000.5"};
        
        Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
        int NN = readEdgesFromFile(graph, "gridtest.txt");
        
        //testMinMaxMovesetHillclimb(false, 0, 10000f, 1000f, 1000f, 1000f, 1f, 1f, 1f, 0f, graph, NN, false, true, 20);
        //testMinMaxMovesetAnnealing(false, 0, 100f,0.99999f, 10000f, 1000f, 1000f, 1000f, 1f, 1f, 1f, 0f, graph, NN, false, true, 20);
        testMinMaxMovesetTabu(false, 0, 543,10f, 10000f, 1000f, 1000f, 1000f, 10000f, 100f, 1f, 0f, graph, NN, false, true, 20);
        //testMinMaxMovesetHillclimb(false, 0, 10000f, 1000f, 1000f, 1000f, 1f, 1f, 1f, 0f, graph, NN, false, false, 20);
        //testMinMaxMovesetAnnealing(false, 0, 100f,0.99999f, 10000f, 1000f, 1000f, 1000f, 1f, 1f, 1f, 0f, graph, NN, false, false, 20);
        //testMinMaxMovesetTabu(false, 0, 7,10f, 10000f, 1000f, 1000f, 1000f, 1f, 1f, 1f, 0f, graph, NN, false, false, 20);

        
        
        /*for (String s:ss){
            Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
            int NN = readEdgesFromFile(graph, s);
            System.out.print(s+"\t");
            testMinMaxMovesetTabu(true, 0, graph.edgeSet().size()/10+1,10f, 10000f, 1000f, 1000f, 1000f, 1f, 1f, 1f, 0f, graph, NN, false, false, 20);
            
        }
        for (String s:sss){
            Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
            int NN = readEdgesFromFile(graph, s);
            System.out.print(s+"\t");
            testMinMaxMovesetTabu(true, 0,graph.edgeSet().size()/10+1,10f, 10000f, 1000f, 1000f, 1000f, 1f, 1f, 1f, 0f, graph, NN, false, false, 20);
            
        }*/
        /*for (String s:testfiles2){
            Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
            int NN = readEdgesFromFile(graph, s);
            System.out.print(s+"\t");
            testMinMaxMovesetTabu(true, 0,graph.edgeSet().size()/10+1,1000f, 1010000f, 101000f, 101000f, 101000f, 1f, 1f, 1f, 0f, graph, NN, false, false, 20);
            
        }*/
        
        //testMinMaxColoringClustering(s, false,false,-1);
        //testMinMaxMatchColoring("gridtest.txt", false, true, 20);
        //for (String s:ss)
        //testMinMaxMatchColoringV2(s,false,false,-1);
        //testMinMaxColoringClustering("5x5.txt", false,true,5);
        //testMinMaxMatchColoring("5x5.txt", false,true,5);
        //testMinMaxMatchColoring("qudg500.2.txt", false, false, 20);

        
        
        
        
        
        
        
        
        
        
        /*String fn[] = {"qudg500.3.txt", "qudg500.2.txt", "udg500.2.txt"};


        for (String filename : fn) {

            Graph<Integer, ColoredEdge> graph = new SimpleGraph<>(ColoredEdge.class);
            int NN = readEdgesFromFile(graph, filename);

            int min = 99999999;
            float temps[] = {10, 100};
            float tempfactors[] = {0.975f, 0.9995f};
            float biases[] = {1000, 10000};
            float weights[] = {0, 1};
            int tabutenures[] = {7, 12};
            float freqpenalties[] = {10};

            float progress = 0;
            int total = (biases.length * biases.length * biases.length * biases.length);
            int ct = 0;
            ArrayList<record> records = new ArrayList<>();
            recordComparator rc = new recordComparator();
            
            

            //testMinMaxMovesetAnnealing(true, 20000, 100, 0.9995f, 10000f, 1000f, 1000f, 1000f, 1f, 1f, 1f, 0f, graph, NN, false, false, 20);
            //testMinMaxMovesetTabu(true,15000, 500, 10, 1000f, 1000f, 1000f, 1000f, 0f, 1f, 1f, 10000f, graph, NN, false, false, 20);
            for (float temp : temps) {
                for (float tempfactor : tempfactors) {
                            for (float w1 : weights) {
                                for (float w4 : weights) {
                                    int v = testMinMaxMovesetAnnealing(false, 8000, temp, tempfactor, 10000f, 1000f, 1000f, 1000f, w1, 1f, 1f, w4, graph, NN, false, false, 20);
                                    if (v < graph.edgeSet().size() / 3) {
                                        record r = new record(v);
                                        r.addEntry("tmp", temp);
                                        r.addEntry("fac", tempfactor);
                                        r.addEntry("w1", w1);
                                        r.addEntry("w4", w4);
                                        records.add(r);
                                    }
                                }
                            
                        }
                    }
                }
            

            Collections.sort(records, rc);
            for (record rr : records) {
                rr.print();
            }
        }*/
        
        
        //testMinMaxMovesetAnnealing(0,100f,0.99999f,100,0,0,0,0,1,0,0,"qudg500.2.txt", false, false, 20);
        //solve5x5("5x5test.txt");

    }

}

