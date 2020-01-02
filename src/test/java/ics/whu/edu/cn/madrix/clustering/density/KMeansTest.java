/**
 *
 */
package ics.whu.edu.cn.madrix.clustering.density;

import java.io.IOException;

import ics.whu.edu.cn.madrix.clustering.evaluation.ClusteringMetrics;
import ics.whu.edu.cn.madrix.clustering.kmeans.KMeans;
import ics.whu.edu.cn.madrix.clustering.wapper.ClusteringBenchTranslator;

/**
 * @author Administrator
 *
 */
public class KMeansTest {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        if (args.length != 4) {
            System.out.println("Usage: command path clusterNumber epsilon isDTW");
            System.exit(-1);
        }
        ClusteringBenchTranslator cbt = new ClusteringBenchTranslator(args[0]);
        KMeans ds = new KMeans(cbt.getData(), Integer.parseInt(args[1]), Double.parseDouble(args[2]),
                Boolean.parseBoolean(args[3]));
        ds.action();
        int[] labels = ds.export();
        ClusteringMetrics cm = new ClusteringMetrics(cbt.getLabels(), labels);
        System.out.println("Acc: " + cm.getAcc() + " Ari: " + cm.getAri() + " Ami: " + cm.getAmi());
    }

}
