/**
 * 
 */
package ics.whu.edu.cn.madrix.clustering.evaluation;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import ics.whu.edu.cn.madrix.stream.utils.OrderInformation;

/**
 * @author Administrator
 *
 */
public class ClusteringMetrics {
    private final Map<Integer, Integer> groundTruth;

    private final Map<Integer, Integer> clusters;

    private Map<Integer, Set<Integer>> testingMap = null;
    private Map<Integer, Set<Integer>> groundTruthMap = null;

    private Map<Integer, Integer> lTogMap = new HashMap<>();
    private Map<Integer, Integer> lTogTrue = new HashMap<>();

    private List<Integer> sequence = null;
    private List<Integer> labelSeq = null;

    private double acc = .0;

    // Specified for AMI.
    private double ami = .0;
    private int[][] commons;

    // Specified for ARI.
    private double ari = .0;
    private double[][] pp;
    private double[][] pn;
    private double[][] np;
    private double[][] nn;

    private int total = 0;

    public ClusteringMetrics(final Map<Integer, Integer> groundTruth, final Map<Integer, Integer> clusters) {
        this.groundTruth = groundTruth;
        this.clusters = clusters;
        initialize();
        computeAcc();
        computeAri();
        computeAmi();
    }

    public ClusteringMetrics(final int[] ground, final int[] labels) {
        groundTruth = new HashMap<>();
        clusters = new HashMap<>();
        for (int i = 0; i < labels.length; i++) {
            groundTruth.put(i, ground[i]);
            clusters.put(i, labels[i]);
        }
        initialize();
        computeAcc();
        computeAri();
        computeAmi();
    }

    private void computeAmi() {
        double HU = .0;
        double HV = .0;
        double HUV = .0;
        double HUCV = .0;
        double IUV = .0;

        for (Entry<Integer, Set<Integer>> cpair : testingMap.entrySet()) {
            double cperc = ((double) cpair.getValue().size()) / total;
            HV += cperc * Math.log(cperc);

        }
        for (Entry<Integer, Set<Integer>> gpair : groundTruthMap.entrySet()) {
            double gperc = ((double) gpair.getValue().size()) / total;
            HU += gperc * Math.log(gperc);
            for (Entry<Integer, Set<Integer>> cpair : testingMap.entrySet()) {
                int overlap = 0;
                for (int id : gpair.getValue()) {
                    if (cpair.getValue().contains(id)) {
                        overlap++;
                    }
                }
                double perc = (overlap == 0) ? 1 : ((double) overlap) / total;
                HUV += perc * Math.log(perc);
                HUCV += perc * ((overlap == 0) ? 0 : Math.log((double) overlap / cpair.getValue().size()));
                IUV += perc * ((overlap == 0) ? 0
                        : Math.log(((double) overlap * total) / (gpair.getValue().size() * cpair.getValue().size())));
            }
        }
        HU = -HU;
        HV = -HV;
        HUV = -HUV;
        HUCV = -HUCV;
        double EIUV1 = .0;
        double EIUV = .0;
        for (Entry<Integer, Set<Integer>> gpair : groundTruthMap.entrySet()) {
            for (Entry<Integer, Set<Integer>> cpair : testingMap.entrySet()) {
                int lsize = gpair.getValue().size();
                int rsize = cpair.getValue().size();
                double perc = ((double) lsize * rsize) / (total * total);
                double logp = (lsize <= 1 || rsize <= 1) ? 0
                        : Math.log(((double) total * (lsize - 1) * (rsize - 1)) / ((total - 1) * lsize * rsize)
                                + (double) total / (lsize * rsize));
                EIUV1 += perc * logp;
                /*int begin = Math.max(lsize + rsize - total, 1);
                int end = Math.min(lsize, rsize);
                for (int i = begin; i <= end; i++) {
                    double perc = (double) i / total;
                    double logp = Math.log(((double) total * i) / (lsize * rsize));
                    double nume = factorials(lsize) * factorials(rsize) * factorials(total - lsize)
                            * factorials(total - rsize);
                    double dnum = factorials(total) * factorials(i) * factorials(lsize - i) * factorials(rsize - i)
                            * factorials(total - lsize - rsize + i);
                    if (dnum == 0 || (perc * logp * nume / dnum) == Double.NaN) {
                        System.out.println("dnum equals to zero");
                    }
                    EIUV += perc * logp * nume / dnum;
                }*/
            }
        }
        EIUV1 = Math.log(
                ((double) total + groundTruthMap.size() * testingMap.size() - groundTruthMap.size() - testingMap.size())
                        / (total - 1));
        if (EIUV1 > IUV) {
            EIUV = Math.log(((double) total + groundTruthMap.size() * testingMap.size() - groundTruthMap.size()
                    - testingMap.size()) / (total - 1));
        } else {
            EIUV = EIUV1;
        }
        //IUV = -IUV;
        System.out.println("HU: " + HU + " HV: " + HV + "\nHUV: " + HUV + " HUCV: " + HUCV + "\nIUV: " + IUV + " EIUV: "
                + EIUV + " EIUV1: " + EIUV1 + "\nNMI: " + IUV / HUV);
        ami = (IUV - EIUV) / (Math.max(HU, HV) - EIUV);
    }

    @SuppressWarnings("unused")
    private double factorials(int max) {
        double fact = 1.0;
        for (int i = 1; i <= max; i++) {
            fact *= i;
        }
        return fact;
    }

    private void computeAri() {
        BitSet bsg = new BitSet(groundTruth.size() * groundTruth.size());
        List<Integer> cidgList = new ArrayList<>(OrderInformation.sortByKey(groundTruth).values());
        for (int i = 0; i < cidgList.size(); i++) {
            int cLeft = groundTruth.get(i);
            for (int j = i + 1; j < groundTruth.size(); j++) {
                if (cLeft == cidgList.get(j)) {
                    bsg.set(i * groundTruth.size() + j);
                } else {
                    bsg.clear(i * groundTruth.size() + j);
                }
            }
        }
        BitSet csg = new BitSet(clusters.size() * clusters.size());
        List<Integer> cidcList = new ArrayList<>(OrderInformation.sortByKey(clusters).values());
        for (int i = 0; i < clusters.size(); i++) {
            int cLeft = clusters.get(i);
            for (int j = i + 1; j < clusters.size(); j++) {
                if (cLeft == cidcList.get(j)) {
                    csg.set(i * clusters.size() + j);
                } else {
                    csg.clear(i * clusters.size() + j);
                }
            }
        }
        int PP = 0;
        int PN = 0;
        int NP = 0;
        int NN = 0;
        for (int i = 0; i < clusters.size(); i++) {
            for (int j = i + 1; j < clusters.size(); j++) {
                int offset = i * clusters.size() + j;
                if (bsg.get(offset)) {
                    if (csg.get(offset)) {
                        PP++;
                    } else {
                        PN++;
                    }
                } else {
                    if (csg.get(offset)) {
                        NP++;
                    } else {
                        NN++;
                    }
                }
            }
        }
        ari = (2 * ((double) NN * PP - NP * PN))
                / (((double) NN + NP) * ((double) NP + PP) + ((double) NN + PN) * ((double) PN + PP));
    }

    private void initialize() {
        this.commons = new int[groundTruth.size()][];
        for (int i = 0; i < groundTruth.size(); i++) {
            commons[i] = new int[clusters.size()];
        }
        this.pp = new double[groundTruth.size()][];
        for (int i = 0; i < groundTruth.size(); i++) {
            pp[i] = new double[clusters.size()];
        }
        this.pn = new double[groundTruth.size()][];
        for (int i = 0; i < groundTruth.size(); i++) {
            pn[i] = new double[clusters.size()];
        }
        this.np = new double[groundTruth.size()][];
        for (int i = 0; i < groundTruth.size(); i++) {
            np[i] = new double[clusters.size()];
        }
        this.nn = new double[groundTruth.size()][];
        for (int i = 0; i < groundTruth.size(); i++) {
            nn[i] = new double[clusters.size()];
        }
    }

    private void computeAcc() {
        testingMap = extract(clusters);
        groundTruthMap = extract(groundTruth);
        Map<Integer, Integer> gSizes = new HashMap<>();

        for (Entry<Integer, Set<Integer>> lp : groundTruthMap.entrySet()) {
            gSizes.put(lp.getKey(), lp.getValue().size());
            total += lp.getValue().size();
        }
        sequence = OrderInformation.keySortByValue(gSizes, false);

        Map<Integer, Integer> lSizes = new HashMap<>();
        for (Entry<Integer, Set<Integer>> lp : testingMap.entrySet()) {
            lSizes.put(lp.getKey(), lp.getValue().size());
        }
        labelSeq = OrderInformation.keySortByValue(lSizes, false);

        //for (Entry<Integer, Set<Integer>> lp : groundTrueMap.entrySet()) {
        for (int gId : sequence) {
            Set<Integer> gSet = groundTruthMap.get(gId);
            int maxOverlap = -1;
            int gCluster = -1;
            for (Entry<Integer, Set<Integer>> gp : testingMap.entrySet()) {
                if (lTogMap.containsValue(gp.getKey())) {
                    continue;
                }
                int overlap = 0;
                for (Integer lo : gSet) {
                    if (gp.getValue().contains(lo)) {
                        overlap++;
                    }
                }
                if (overlap > maxOverlap) {
                    maxOverlap = overlap;
                    gCluster = gp.getKey();
                }
            }
            lTogMap.put(gId, gCluster);
            lTogTrue.put(gId, maxOverlap);
        }
    }

    public double getAcc() {
        int rights = 0;
        for (Integer right : lTogTrue.values()) {
            rights += right;
        }
        acc = (double) rights / total;
        System.out.println(
                "Right: " + rights + " total: " + total + " maxLabel: " + testingMap.get(labelSeq.get(0)).size());
        return acc;
    }

    public double getAmi() {
        return ami;
    }

    public double getAri() {
        return ari;
    }

    private Map<Integer, Set<Integer>> extract(Map<Integer, Integer> labels) {
        Map<Integer, Set<Integer>> gMap = new HashMap<>();
        for (Entry<Integer, Integer> pair : labels.entrySet()) {
            if (gMap.containsKey(pair.getValue())) {
                gMap.get(pair.getValue()).add(pair.getKey());
            } else {
                Set<Integer> elements = new HashSet<>();
                elements.add(pair.getKey());
                gMap.put(pair.getValue(), elements);
            }
        }
        return gMap;
    }
}
