/**
 *
 */
package ics.whu.edu.cn.madrix.clustering.density;

import ics.whu.edu.cn.madrix.clustering.evaluation.ClusteringMetrics;
import ics.whu.edu.cn.madrix.clustering.wapper.ClusteringBenchTranslator;
import ics.whu.edu.cn.madrix.stream.utils.OrderInformation;

import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;

/**
 * @author Administrator
 *
 */
public class ATDTest {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        if (args.length != 4) {
            System.out.println("Usage: command path clusterNumber k isDTW");
            System.exit(-1);
        }
        ClusteringBenchTranslator cbt = new ClusteringBenchTranslator(args[0]);
        int Ker = Integer.valueOf(args[1]).intValue();
        ATD atd = new ATD(cbt.getData(), Ker, Integer.parseInt(args[2]), Boolean.parseBoolean(args[3]));
        atd.action();
        List<Integer> priority = atd.getPriority();
        for (int i = 0; i < cbt.getData().length; i++) {
            System.out.println(cbt.getData()[priority.get(i)][0] + "\t" + cbt.getData()[priority.get(i)][1] + "\t"
                    + atd.getEntropies().get(priority.get(i)) + "\t" + priority.get(i));
        }
        Map<Integer, List<Integer>> derivedPath = atd.getPathMap();
        Map<Integer, Set<Integer>> pathExpansion = atd.getExpansion();
        for (Entry<Integer, List<Integer>> pair : derivedPath.entrySet()) {
            int key = pair.getKey();
            List<Integer> path = pair.getValue();
            System.out.println(key + 1);
            for (int kid : path) {
                System.out.print("\t" + (kid + 1) + "\t");
                for (int nid : pathExpansion.get(kid)) {
                    System.out.print(nid + 1 + "\t");
                }
                System.out.println();
            }
        }
        Collection<Integer> c = atd.getIndex().values();
        Object[] cindex = c.toArray();
        List<Integer> rhoK = new ArrayList<Integer>();
        for (Object nid : cindex) {
            rhoK.add((Integer) nid);
            System.out.print((Integer) nid + 1 + "\t");
        }
        System.out.print("\n");
        System.out.println("CTD");
        CTD ctd = new CTD(cbt.getData(), atd, rhoK, Ker);

        System.out.println("STD");
        STD std = new STD(ctd, atd, Ker, atd.getDH());
        Map<Integer, Double> resultK = OrderInformation.sortByValue(std.getreStack(), false);
        Collection<Integer> cs = std.getIndex().values();
        Object[] csindex = cs.toArray();
        List<Integer> rhosK = new ArrayList<Integer>();
        for (Object nid : csindex) {
            rhosK.add((Integer) nid);
            System.out.print((Integer) nid + 1 + "\t");
        }
        System.out.print("\n");

        System.out.println("LTD");
        LTD ltd = new LTD(ctd, std, atd, resultK, Ker);
        ltd.action();

        int[] labels = ltd.export();
        if (labels.length > 0) {
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
                if (pair.getValue() >= Integer.parseInt(args[1]))
                    output.add(pair.getKey());
            }
            Map<Integer, Integer> labelMap = new HashMap<>();
            for (int i = 0; i < labels.length; i++) {
                if (ccnt.get(labels[i]) >= Integer.parseInt(args[1])) {
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
            System.out.println("Acc: " + cm.getAcc() + " Ari: " + cm.getAri() + " Ami: " + cm.getAmi());
        }
    }
}
