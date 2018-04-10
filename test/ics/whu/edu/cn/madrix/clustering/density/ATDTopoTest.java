/**
 * 
 */
package ics.whu.edu.cn.madrix.clustering.density;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import ics.whu.edu.cn.madrix.clustering.wapper.ClusteringBenchTranslator;

/**
 * @author Administrator
 *
 */
public class ATDTopoTest {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        ClusteringBenchTranslator cbt =
                new ClusteringBenchTranslator("test/ics/whu/edu/cn/madrix/clustering/resources/Compound3.txt");
        ATD atd = new ATD(cbt.getData(), 6, 2, false);
        List<Integer> pq = atd.getPriority();
        int center = pq.get(15);
        @SuppressWarnings("unused")
        int[] centers = { 0, 1, 3, 7, 15, 26, 27 };
        //System.out.println(center);
        for (int i = 1; i < pq.size(); i++) {
            Set<Integer> nns = new HashSet<Integer>(atd.getNbs().get(center).subList(0, i + 1));
            double maxAvg = atd.entropy(nns);//atd.maxAvgDist(nns);
            System.out.println(maxAvg);
        }
    }

}
