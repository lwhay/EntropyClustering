/**
 * 
 */
package ics.whu.edu.cn.madrix.clustering.density;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import ics.whu.edu.cn.madrix.clustering.evaluation.ClusteringMetrics;
import ics.whu.edu.cn.madrix.clustering.wapper.ClusteringBenchTranslator;

/**
 * @author Administrator
 *
 */
public class DBSCANTest {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
    	if (args.length != 4) {
    		System.out.println("Usage: command path cutoff mpts isDTW");
    		System.exit(-1);
    	}
        ClusteringBenchTranslator cbt = new ClusteringBenchTranslator(args[0]);
        DBSCAN ds = new DBSCAN(cbt.getData(), Double.parseDouble(args[1]), Integer.parseInt(args[2]),
                Boolean.parseBoolean(args[3]));
        ds.action();
        int[] labels = ds.export();
        Map<Integer, Integer> ccnt = new HashMap<>();
        for (int i = 0; i < labels.length; i++) {
            if (ccnt.containsKey(labels[i])) {
                ccnt.put(labels[i], ccnt.get(labels[i]) + 1);
            } else {
                ccnt.put(labels[i], 1);
            }
        }
        Set<Integer> output = new HashSet<>();
        for (Entry<Integer, Integer> pair : ccnt.entrySet()) {
            if (pair.getValue() >= Integer.parseInt(args[2]))
                output.add(pair.getKey());
        }
        Map<Integer, Integer> labelMap = new HashMap<>();
        for (int i = 0; i < labels.length; i++) {
            if (ccnt.get(labels[i]) >= Integer.parseInt(args[2])) {
                labelMap.put(i, labels[i]);
            } else {
                labelMap.put(i, -1);
            }
        }
        Map<Integer, Integer> groundTruth = new HashMap<>();
        int idx = 0;
        for (int label : cbt.getLabels()) {
            groundTruth.put(idx++, label);
        }
        for (int i = 0; i < labels.length; i++) {
            if (output.contains(labels[i]))
                System.out.println(cbt.getData()[i][0] + "\t" + cbt.getData()[i][1] + "\t" + labels[i]);
        }
        ClusteringMetrics cm = new ClusteringMetrics(groundTruth, labelMap);
        System.out.println("Cutoff: " + ds.getCutoff() + " Radius: " + ds.getRadius() + " MaxDist: "
                + ds.getMaximalDist() + "\nGenerated: " + output.size() + " Acc: " + cm.getAcc() + " Ari: "
                + cm.getAri() + " Ami: " + cm.getAmi());
    }
}
