/**
 *
 */
package ics.whu.edu.cn.madrix.clustering.connect;

import ics.whu.edu.cn.madrix.clustering.density.AbstractDensityClustering;
import ics.whu.edu.cn.madrix.clustering.density.IClustering;
import ics.whu.edu.cn.madrix.common.exceptions.MadrixException;
import ics.whu.edu.cn.madrix.stream.utils.OrderInformation;

import java.util.*;

/**
 * @author Administrator
 *
 */
public class Transitive extends AbstractDensityClustering implements IClustering {
    private Map<Integer, List<Integer>> distSortedNbs = new HashMap<>();

    private Map<Integer, Set<Integer>> directTo = new HashMap<>();

    private final double[][] dists;

    private Map<Integer, Integer> labels = new HashMap<>();

    private Map<Integer, Integer> nbCnts = new HashMap<>();

    public Transitive(double[][] data, double cutoff, int type) throws MadrixException {
        this.type = type;
        this.data = data;
        this.cutoff = cutoff;
        rankingdist = new double[data.length * data.length];
        dists = new double[data.length][];
        for (int i = 0; i < data.length; i++) {
            dists[i] = new double[data.length];
            for (int j = 0; j < data.length; j++) {
                dists[i][j] = distance(data[i], data[j], type);
                rankingdist[i * data.length + j] = dists[i][j];
            }
        }
        Arrays.sort(rankingdist);
        radius = rankingdist[(int) (rankingdist.length * cutoff)];
        for (int i = 0; i < data.length; i++) {
            Map<Integer, Double> distMap = new HashMap<>();
            for (int j = 0; j < data.length; j++) {
                distMap.put(j, dists[i][j]);
            }
            distSortedNbs.put(i, OrderInformation.keySortByValue(distMap, true));
            int nCnt = 0;
            for (int j = 0; j < data.length; j++) {
                if (dists[i][distSortedNbs.get(i).get(j)] > radius) {
                    break;
                } else {
                    if (!directTo.containsKey(i)) {
                        directTo.put(i, new HashSet<>());
                    }
                    directTo.get(i).add(distSortedNbs.get(i).get(j));
                }
                nCnt++;
            }
            nbCnts.put(i, nCnt);
        }
    }

    @Override
    public int[] export() {
        int[] results = new int[labels.size()];
        for (int i = 0; i < labels.size(); i++) {
            results[i] = labels.get(i);
        }
        return results;
    }

    private Set<Integer> reached = new HashSet<>();

    @Override
    public void action() {
        int cid = 1;
        for (int i = 0; i < data.length; i++) {
            Set<Integer> related = new HashSet<>();
            if (reached.contains(i)) {
                //labels.put(i, labels.get(derivedFrom.get(i)));
            } else {
                iterExp(i, cid++, related);
            }
        }

    }

    private void iterExp(int id, int cid, Set<Integer> related) {
        related.add(id);
        labels.put(id, cid);
        reached.add(id);
        for (Integer p : directTo.get(id)) {
            if (!reached.contains(p)) {
                iterExp(p, cid, related);
            }
        }
    }
}
