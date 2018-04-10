/**
 * 
 */
package ics.whu.edu.cn.madrix.clustering.density;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ics.whu.edu.cn.madrix.stream.utils.OrderInformation;

/**
 * @author Administrator
 *
 */
public class DBSCAN extends AbstractDensityClustering implements IClustering {
    private final int mpt;

    private Map<Integer, List<Integer>> distSortedNbs = new HashMap<>();

    private Map<Integer, Integer> derivedFrom = new HashMap<>();

    private final double[][] dists;

    private Map<Integer, Integer> labels = new HashMap<>();

    private Map<Integer, Integer> nbCnts = new HashMap<>();

    public DBSCAN(double[][] data, double cutoff, int mpt, boolean isDTW) {
        this.isDTW = isDTW;
        this.data = data;
        this.cutoff = cutoff;
        this.mpt = mpt;
        rankingdist = new double[data.length * data.length];
        dists = new double[data.length][];
        for (int i = 0; i < data.length; i++) {
            dists[i] = new double[data.length];
            for (int j = 0; j < data.length; j++) {
                dists[i][j] = distance(data[i], data[j], isDTW);
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
        int nbs = nbCnts.get(id);
        if (nbs < mpt) {
            return;
        }
        for (int i = 0; i < nbCnts.get(id); i++) {
            int nb = distSortedNbs.get(id).get(i);
            if (!related.contains(nb)) {
                iterExp(nb, cid, related);
                derivedFrom.put(nb, id);
            }
        }
    }
}
