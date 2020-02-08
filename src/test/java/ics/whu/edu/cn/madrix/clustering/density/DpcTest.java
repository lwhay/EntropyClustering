/**
 *
 */
package ics.whu.edu.cn.madrix.clustering.density;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import ics.whu.edu.cn.madrix.clustering.density.DpcDecision;
import ics.whu.edu.cn.madrix.clustering.evaluation.ClusteringMetrics;
import ics.whu.edu.cn.madrix.clustering.wapper.ClusteringBenchTranslator;
import ics.whu.edu.cn.madrix.common.exceptions.MadrixException;

/**
 * @author Administrator
 *
 */
public class DpcTest {
    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException, MadrixException {
        if (args.length != 4) {
            System.out.println("Usage: command path clusterNumber cutoff isDTW");
            System.exit(-1);
        }
        // TODO Auto-generated method stub
        ClusteringBenchTranslator cbt = new ClusteringBenchTranslator(args[0]);
        DpcDecision dd = new DpcDecision(cbt.getData(), Integer.parseInt(args[1]), Double.parseDouble(args[2]),
                Integer.parseInt(args[3]));
        dd.action();
        int[] labels = dd.export();

        Map<Integer, Integer> cnts = new HashMap<>();
        for (int i = 0; i < labels.length; i++) {
            if (cnts.containsKey(labels[i])) {
                cnts.put(labels[i], cnts.get(labels[i]) + 1);
            } else {
                cnts.put(labels[i], 1);
            }
        }
        for (int peak : dd.peaks()) {
            System.out.println(peak + "\t" + labels[peak] + "\t" + cnts.get(labels[peak]));
        }
        for (int i = 0; i < labels.length; i++) {
            System.out.println(cbt.getData()[i][0] + "\t" + cbt.getData()[i][1] + "\t" + labels[i]);
        }
        ClusteringMetrics cm = new ClusteringMetrics(cbt.getLabels(), labels);
        System.out.println("Cutoff: " + dd.getCutoff() + " Radius: " + dd.getRadius() + " MaxDist: "
                + dd.getMaximalDist() + " Acc: " + cm.getAcc() + " Ari: " + cm.getAri() + " Ami: " + cm.getAmi());
    }
}
