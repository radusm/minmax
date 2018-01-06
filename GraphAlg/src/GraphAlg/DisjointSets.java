/* DisjointSets.java */

package GraphAlg;


/**
 *  A disjoint sets ADT.  Performs union-by-rank and path compression.
 *  Implemented using arrays.  There is no error checking whatsoever.
 *  By adding your own error-checking, you might save yourself a lot of time
 *  finding bugs in your application code for Project 3 and Homework 9.
 *  Without error-checking, expect bad things to happen if you try to unite
 *  two elements that are not roots of their respective sets, or are not
 *  distinct.
 *
 *  Elements are represented by ints, numbered from zero.
 *
 *  @author Mark Allen Weiss
 **/


class DisjointSets
{
    int[] rank, parent,number;
    int n;
    int maxnumber=-1;

    // Constructor
    public DisjointSets(int n)
    {
        rank = new int[n];
        parent = new int[n];
        number = new int[n];

        
        this.n = n;
        makeSet();
    }

    @Override
    public String toString(){
    String str="";
    for (int i=0;i<n;i++)
    str+="{t["+i+"]="+parent[i]+";cnt="+number[i]+"\n";
    return str;
    }
    
    public int getNumber(int i){
        if (i<n) return number[find(i)];
        return -1;
    }
    
    // Creates n sets with single item in each
    void makeSet()
    {
        for (int i=0; i<n; i++)
        {
            // Initially, all elements are in
            // their own set.
            parent[i] = i;
            rank[i]=0;
            number[i]=1;
        }
    }

    // Returns representative of x's set
    int find(int x)
    {
        // Finds the representative of the set
        // that x is an element of
        if (parent[x]!=x)
        {
            // if x is not the parent of itself
            // Then x is not the representative of
            // his set,
            parent[x] = find(parent[x]);

            // so we recursively call Find on its parent
            // and move i's node directly under the
            // representative of this set
        }

        return parent[x];
    }

    // Unites the set that includes x and the set
    // that includes x
    int union(int x, int y)
    {
        // Find representatives of two sets
        int xRoot = find(x), yRoot = find(y);

        // Elements are in the same set, no need
        // to unite anything.
        if (xRoot == yRoot){
            number[xRoot]++;//mycomment: adding an edge inside the same connected comp
            if (number[xRoot]>maxnumber) maxnumber=number[xRoot];
            return xRoot;
        }
         // If x's rank is less than y's rank
        if (rank[xRoot] < rank[yRoot])

            // Then move x under y  so that depth
            // of tree remains less
            {parent[xRoot] = yRoot;
            number[yRoot]+=number[xRoot]+1;//mycomment: adding an edge between these 2 components
            if (number[yRoot]>maxnumber) maxnumber=number[yRoot];
            return yRoot;
        }

        // Else if y's rank is less than x's rank
        else if (rank[yRoot] < rank[xRoot])

            // Then move y under x so that depth of
            // tree remains less
            {parent[yRoot] = xRoot;
            number[xRoot]+=number[yRoot]+1;//mycomment: adding an edge between these 2 components
            if (number[xRoot]>maxnumber) maxnumber=number[xRoot];
            return xRoot;
        }

        // if ranks are the same

            // Then move y under x (doesn't matter
            // which one goes where)
            parent[yRoot] = xRoot;
            number[xRoot]+=number[yRoot]+1;//mycomment: adding an edge between these 2 components
            if (number[xRoot]>maxnumber) maxnumber=number[xRoot];
            // And increment the the result tree's
            // rank by 1
            rank[xRoot] = rank[xRoot] + 1;
            return xRoot;
        
    }
}
