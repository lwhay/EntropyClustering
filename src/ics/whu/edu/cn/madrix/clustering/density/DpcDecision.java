/**
 * 
 */
package ics.whu.edu.cn.madrix.clustering.density;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ics.whu.edu.cn.madrix.stream.utils.OrderInformation;

/**
 * @author Administrator
 *
 */
public class DpcDecision extends AbstractDensityClustering implements IClustering {
    private double[][] dist;

    private int K;

    private Map<Integer, Integer> dependency = new HashMap<>();

    private List<Integer> sequence = new ArrayList<>();

    List<Map.Entry<Integer, Integer>> density = new ArrayList<>();

    List<Map.Entry<Integer, Double>> decision = new ArrayList<>();

    public DpcDecision(double[][] data, int K, double cutoff, boolean isDTW) {
        this.isDTW = isDTW;
        this.data = data;
        this.K = K;
        this.ouputLabels = new int[data.length];
        this.rankingdist = new double[data.length * data.length];
        this.dist = new double[data.length][];
        this.cutoff = cutoff;
        int sidx = 0;
        for (double[] s : data) {
            int tidx = 0;
            dist[sidx] = new double[data.length];
            for (double[] t : data) {
                dist[sidx][tidx] = distance(s, t, isDTW);
                rankingdist[sidx * data.length + tidx] = dist[sidx][tidx];
                tidx++;
            }
            sidx++;
        }
        Arrays.sort(rankingdist);
        radius = rankingdist[rankingdist.length - 1] * cutoff;
    }

    @Override
    public void action() {
        Map<Integer, Integer> denest = new HashMap<>();
        //Map<Integer, List<Integer>> distSort = new HashMap<>();
        for (int i = 0; i < data.length; i++) {
            int cnt = 0;
            //Map<Integer, Double> dists = new HashMap<>();
            //int oid = 0;
            for (int j = 0; j < data.length; j++) {
                double[] pt = data[j];
                //dists.put(oid++, distance(pt, data[i]));
                if (distance(pt, data[i], isDTW) <= radius) {
                    cnt++;
                }
            }
            //List<Integer> sortedDists = OrderInformation.keySortByValue(dists, true);
            //distSort.put(i, new ArrayList<Integer>(sortedDists));
            denest.put(i, cnt);
        }
        density = OrderInformation.sortByValues(denest, false);
        sequence = OrderInformation.keySortByValue(denest, false);
        //for (Entry<Integer, Integer> pair : density) {
        //System.out.println(data[pair.getKey()][0] + "\t" + data[pair.getKey()][1] + "\t" + pair.getValue());
        //}
        Map<Integer, Double> mindist = new HashMap<>();
        int idx = 0;
        for (int oid : sequence) {
            if (idx == 0) {
                dependency.put(oid, -1);
                mindist.put(oid, rankingdist[rankingdist.length - 1]);
            } else {
                double dt = Double.MAX_VALUE;
                int nid = -1;
                for (int i = 0; i < idx; i++) {
                    double cdt = distance(data[sequence.get(i)], data[oid], isDTW);
                    if (cdt < dt) {
                        dt = cdt;
                        nid = sequence.get(i);
                    }
                }
                dependency.put(oid, nid);
                mindist.put(oid, dt);
            }
            idx++;
        }
        decision = OrderInformation.sortByValues(mindist, false);
    }

    public int[] peaks() {
        int[] peaks = new int[K];
        int idx = 0;
        for (Map.Entry<Integer, Double> pair : decision) {
            if (idx >= K) {
                break;
            }
            peaks[idx++] = pair.getKey();
        }
        return peaks;
    }

    @Override
    public int[] export() {
        Map<Integer, Integer> results = new HashMap<>();
        int[] pks = peaks();
        for (int i = 0; i < pks.length; i++) {
            results.put(pks[i], i + 1);
        }
        for (int oid : sequence) {
            //System.out.println(oid + "\t" + dependency.get(oid));
            if (results.containsKey(dependency.get(oid))) {
                if (!results.containsKey(oid)) {
                    results.put(oid, results.get(dependency.get(oid)));
                }
            } /*else {
                System.out.println("^-^");
                //results.put(oid, results.get(dep.getValue()));
              }*/
        }
        for (int i = 0; i < data.length; i++) {
            ouputLabels[i] = results.get(i);
        }
        return ouputLabels;
    }
}
