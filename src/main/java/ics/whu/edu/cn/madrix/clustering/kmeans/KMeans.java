/**
 *
 */
package ics.whu.edu.cn.madrix.clustering.kmeans;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import ics.whu.edu.cn.madrix.clustering.density.AbstractDensityClustering;
import ics.whu.edu.cn.madrix.clustering.density.IClustering;
import ics.whu.edu.cn.madrix.common.exceptions.MadrixException;

/**
 * @author Administrator
 *
 */
public class KMeans implements IClustering {
    private static final int DEFAULT_TIMEOUT = 100;

    private int round = 0;

    protected int type = -1;

    protected double[][] data;

    protected int[] ouputLabels;

    protected int dim;

    @SuppressWarnings("unused")
    private final int K;

    private final Map<Integer, Set<Integer>> clusters;

    private final double[][] dists;

    private final double[][] kernels;

    private double globalMaxDist = Double.MIN_VALUE;

    private final double epsilon;

    public KMeans(double[][] data, int K, double epsilon, int type) throws MadrixException {
        this.type = type;
        this.data = data;
        this.epsilon = epsilon;
        this.K = K;
        this.clusters = new HashMap<>();
        this.kernels = new double[K][];
        for (int i = 0; i < K; i++) {
            kernels[i] = new double[data[0].length];
            int kidx = /*data.length / (2 * K) +*/ data.length * i / K;
            //int kidx = (int) (Math.random() * (data.length - 1));
            while (this.clusters.containsKey(kidx)) {
                kidx = (int) (Math.random() * (data.length - 1));
            }
            clusters.put(i, new HashSet<Integer>());
            for (int j = 0; j < data[kidx].length; j++) {
                kernels[i][j] = data[kidx][j];
            }
        }
        dists = new double[data.length][];
        for (int i = 0; i < data.length; i++) {
            dists[i] = new double[data.length];
            for (int j = 0; j < data.length; j++) {
                dists[i][j] = AbstractDensityClustering.distance(data[i], data[j], type);
                if (dists[i][j] > globalMaxDist) {
                    globalMaxDist = dists[i][j];
                }
            }
        }
    }

    @Override
    public int[] export() {
        int[] labels = new int[data.length];
        for (Entry<Integer, Set<Integer>> pair : clusters.entrySet()) {
            for (int oid : pair.getValue()) {
                labels[oid] = pair.getKey();
            }
        }
        return labels;
    }

    @Override
    public void action() throws MadrixException {
        double avgEps = Double.MAX_VALUE;
        do {
            toClusters();
            avgEps = updateKernels();
        } while (round++ < DEFAULT_TIMEOUT && avgEps > epsilon * globalMaxDist);
        System.out.println("Round: " + round + " avgEps: " + avgEps);
    }

    private void toClusters() throws MadrixException {
        for (int i = 0; i < data.length; i++) {
            double minDist = Double.MAX_VALUE;
            int belongTo = -1;
            for (int k = 0; k < kernels.length; k++) {
                double dist = AbstractDensityClustering.distance(data[i], kernels[k], type);
                if (dist < minDist) {
                    minDist = dist;
                    belongTo = k;
                }
            }
            clusters.get(belongTo).add(i);
        }
    }

    private double updateKernels() throws MadrixException {
        double epsilon = .0;
        for (Entry<Integer, Set<Integer>> pair : clusters.entrySet()) {
            int kid = pair.getKey();
            System.out.print(kid + "\t");
            for (int j = 0; j < kernels[kid].length; j++) {
                System.out.print(kernels[kid][j] + "\t");
            }
            System.out.println();
            double[] oldPos = Arrays.copyOf(kernels[kid], kernels[kid].length);
            for (int i = 0; i < kernels[kid].length; i++) {
                kernels[kid][i] = .0;
            }
            for (int oid : pair.getValue()) {
                for (int i = 0; i < kernels[kid].length; i++) {
                    kernels[kid][i] += data[oid][i];
                }
            }
            for (int i = 0; i < kernels[kid].length; i++) {
                kernels[kid][i] /= pair.getValue().size();
            }
            epsilon += AbstractDensityClustering.distance(oldPos, kernels[kid], type);
        }
        System.out.println();
        return epsilon /= clusters.size();
    }
}
