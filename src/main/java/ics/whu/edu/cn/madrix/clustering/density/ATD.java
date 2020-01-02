/**
 *
 */
package ics.whu.edu.cn.madrix.clustering.density;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Set;

import ics.whu.edu.cn.madrix.stream.utils.OrderInformation;
import org.apache.commons.math3.special.Gamma;

/**
 * @author Administrator
 *
 */
public class ATD extends AbstractDensityClustering implements IClustering {
    private final boolean isDTW;

    private final int KT;

    static private final int EXPENSIONDEGREE = 2;

    private Map<Integer, List<Integer>> distSortedNbs = new HashMap<>();

    private double[][] dists;

    private int freshid = 0;

    //private int KKm = 0;
    private int DH = 0;

    private Map<Integer, Integer> Index = new HashMap<>();

    private Map<Integer, Double> initEntropies = new HashMap<>();

    private List<Integer> priority = new ArrayList<>();

    private final double BALL;

    public ATD(double[][] data, int K, int k, boolean isDTW) {
        this.isDTW = isDTW;
        this.KT = k;
        this.data = data;
        this.dim = isDTW ? 2 : data[0].length;
        this.BALL = Math.pow(Math.PI, dim / 2) / Gamma.gamma(dim / 2 + 1);
        rankingdist = new double[data.length * data.length];
        dists = new double[data.length][];
        for (int i = 0; i < data.length; i++) {
            dists[i] = new double[data.length];
            Map<Integer, Double> distMap = new HashMap<>();
            for (int j = 0; j < data.length; j++) {
                dists[i][j] = distance(data[i], data[j], isDTW);
                rankingdist[i * data.length + j] = dists[i][j];
                distMap.put(j, dists[i][j]);
            }
            distSortedNbs.put(i, OrderInformation.keySortByValue(distMap, true));
            Set<Integer> nn = new HashSet<>();
            for (int j = 0; j <= KT; j++) { //KT * DH
                nn.add(distSortedNbs.get(i).get(j));
            }
            initEntropies.put(i, entropy(nn));
        }
        priority = OrderInformation.keySortByValue(initEntropies, false);
        for (int i = 0; i < data.length; i++) {
            Index.put(priority.get(i).intValue(), 0);
        }
        Arrays.sort(rankingdist);
    }

    @Override
    public int[] export() {
        // TODO Auto-generated method stub
        return null;
    }

    private Set<Integer> reached = new HashSet<>();

    private Map<Integer, List<Integer>> derivedPath = new HashMap<>();

    private Map<Integer, Set<Integer>> pathExpansion = new HashMap<>();

    @Override
    public void action() {
        for (int i = 0; i < priority.size(); i++) {
            if (!reached.contains(priority.get(i))) {
                List<Integer> path = new ArrayList<>();
                path.add(priority.get(i));
                intendFreshKT(priority.get(i), path);
                derivedPath.put(priority.get(i), path);
            }
        }
    }

    private double iterateFrom(int oid, Set<Integer> pathNNs) {
        double nnAvg = nnMaxavg(oid);
        levelNNs(oid, oid, nnAvg, pathNNs, EXPENSIONDEGREE);
        return nnAvg;
    }

    private void levelNNs(int ooid, int oid, double nnAvg, Set<Integer> pathNNs, int nnLevel) {
        Set<Integer> curNNs = new HashSet<>();
        Set<Integer> MaxavgDistNN = maxavgNN(oid, nnMaxavg(freshid)); // nnAvg);
        for (int noid : MaxavgDistNN) {
            if (!pathNNs.contains(noid) && !curNNs.contains(noid)) {
                //if (reachable(oid, noid, nnAvg)) {
                if ( /*initEntropies.get(noid) <= initEntropies.get(ooid) &&  */expansible(ooid, noid)) {
                    // && reachable(oid, noid, nnAvg) && expansible(ooid, noid)) {
                    pathNNs.add(noid);
                    curNNs.add(noid);
                    reached.add(noid);
                }
            }
        }
        /*if (--nnLevel > 0) {
            for (int noid : curNNs) {
                levelNNs(ooid, noid, nnAvg, pathNNs, nnLevel);
            }
        }*/
    }

    private void intendFreshKT(int coid, List<Integer> path) {
        freshid = coid;
        Set<Integer> pathNNs = new HashSet<>();
        pathNNs.add(coid);
        pathNNs.addAll(path);
        reached.add(coid);
        Set<Integer> oldPathNNs = new HashSet<>(pathNNs);
        iterateFrom(coid, pathNNs);
        Set<Integer> newPathNNs = new HashSet<>(pathNNs);
        newPathNNs.removeAll(oldPathNNs);
        newPathNNs.add(coid);
        pathExpansion.put(coid, newPathNNs);
        freshid = peakNNs(pathNNs);
        //int freshid = peakEntropyNNs(pathNNs, nnAvg);
        //int freshid = peakEntropyNNs(pathNNs, nnAvg);
        Integer[] temp = newPathNNs.toArray(new Integer[]{});
        if (freshid != coid) {
            path.add(freshid);
            for (int i = 0; i < newPathNNs.size(); i++) {
                if (Index.get(temp[i].intValue()) != 0) {
                    List<Integer> cindex = getKeyByValues(Index, Index.get(temp[i].intValue()), temp[i].intValue());
                    // List<Integer> cindex = getKeyByValue(Index , Index.get(temp[i].intValue()));
                    for (int j = 0; j < cindex.size(); j++) {
                        Index.put(cindex.get(j).intValue(), freshid);
                    }
                    // Collection<Integer> c = Index.values();
                    // Integer[] cindex = (Integer[]) c.toArray();

                } else {
                    Index.put(temp[i].intValue(), freshid);
                }
            }
            intendFreshKT(freshid, path);
        } else {
            for (int i = 0; i < newPathNNs.size(); i++) {
                Index.put(temp[i].intValue(), freshid);
            }
        }
    }

    @SuppressWarnings("unused")
    private int peakEntropyNNs(Set<Integer> pathNNs, double nnAvg) {
        int peak = -1;
        double minEntropy = Double.MAX_VALUE;
        for (int oid : pathNNs) {
            Set<Integer> nns = maxavgNN(oid, nnAvg);
            double curEntropy = entropy(nns);
            if (curEntropy < minEntropy) {
                minEntropy = curEntropy;
                peak = oid;
            }
        }
        return peak;
    }

    public List<Integer> getKeyByValues(Map<Integer, Integer> map, Object value, int m) {
        List<Integer> keys = new ArrayList<Integer>();
        Iterator<Entry<Integer, Integer>> it = map.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<Integer, Integer> entry = (Entry<Integer, Integer>) it.next();
            Object obj = entry.getValue();
            if (obj != null && obj.equals(value) && initEntropies.get(m) <= initEntropies.get(obj)) {
                keys.add((Integer) entry.getKey());
            }
        }
        return keys;
    }

    public List<Integer> getKeyByValue(Map<Integer, Integer> map, Object value) {
        List<Integer> keys = new ArrayList<Integer>();
        Iterator<Entry<Integer, Integer>> it = map.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<Integer, Integer> entry = (Entry<Integer, Integer>) it.next();
            Object obj = entry.getValue();
            if (obj != null && obj.equals(value)) {
                keys.add((Integer) entry.getKey());
            }
        }
        return keys;
    }

    private int peakNNs(Set<Integer> pathNNs) {
        int peak = -1;
        double minEntropy = Double.MAX_VALUE;
        for (int oid : pathNNs) {
            double curEntropy = initEntropies.get(oid);
            if (curEntropy < minEntropy) {
                minEntropy = curEntropy;
                peak = oid;
            }
        }
        for (int oid : pathNNs) {
            System.out.print(oid + "\t");
        }
        System.out.println(peak);
        return peak;
    }

    @SuppressWarnings("unused")
    private int peakNNs(Set<Integer> pathNNs, double nnAvg) {
        int peak = -1;
        double minEntropy = Double.MAX_VALUE;
        for (int oid : pathNNs) {
            Set<Integer> nns = new HashSet<>();
            for (int noid : distSortedNbs.get(oid)) {
                if (dists[oid][noid] <= nnAvg) {
                    nns.add(noid);
                } else {
                    break;
                }
            }
            double curEntropy = entropy(nns);
            if (curEntropy < minEntropy) {
                minEntropy = curEntropy;
                peak = oid;
            }
        }
        return peak;
    }

    private double nnMaxavg(int center) {
        Set<Integer> nns = new HashSet<>();
        for (int i = 0; i < KT + 1; i++) {
            nns.add(distSortedNbs.get(center).get(i));
        }
        return maxAvgDist(nns);
    }

    @SuppressWarnings("unused")
    private Set<Integer> minavgNN(int center, double radius) {
        Set<Integer> bs = new HashSet<>();
        for (int nn : distSortedNbs.get(center)) {
            bs.add(nn);
            if (maxAvgDist(bs) > radius) {
                break;
            }
        }
        return bs;
    }

    private Set<Integer> maxavgNN(int center, double radius) {
        Set<Integer> bs = new HashSet<>();
        for (int nn : distSortedNbs.get(center)) {
            bs.add(nn);
            if (maxAvgDist(bs) > radius) {
                bs.remove(nn);
                break;
            }
        }
        return bs;
    }

    public Map<Integer, List<Integer>> getPathMap() {
        return derivedPath;
    }

    public Map<Integer, Set<Integer>> getExpansion() {
        return pathExpansion;
    }

    // For debug purpose.
    public Map<Integer, List<Integer>> getNbs() {
        return this.distSortedNbs;
    }

    public double maxAvgDist(Set<Integer> bs) {
        double sumMaxDist = .0;
        for (int bid : bs) {
            double maxDist = .0;
            for (int sid : bs) {
                if (dists[bid][sid] > maxDist) {
                    maxDist = dists[bid][sid];
                }
            }
            sumMaxDist += maxDist;
        }
        return sumMaxDist / bs.size();
    }

    public double entropy(Set<Integer> pts) {
        double volume = .0f;
        for (int me : pts) {
            double md = .0;
            for (int you : pts) {
                double d = distance(data[me], data[you], isDTW);
                if (d > md) {
                    md = d;
                }
            }
            double revolume = Math.pow(md, dim);
            if (revolume > Double.MAX_VALUE) {
                System.out.println(md + ":" + dim);
                revolume = Double.MAX_VALUE;
            }
            volume += 1 / revolume;
        }
        //double entropy = -Math.log(1 / (volume * (pts.size() - 1) * BALL / pts.size()));
        double entropy = -Math.log(volume / ((pts.size()) * BALL));
        return entropy;
    }

    // To be used for kNN-based entropy.
    public double entropy(Set<Integer> pts, int k) {
        double volume = .0f;
        for (int me : pts) {
            double md = Double.MAX_VALUE;
            FastBoundedPriorityQueue<Double> queue = new FastBoundedPriorityQueue<>(k, true);
            for (int you : pts) {
                double d = distance(data[me], data[you], isDTW);
                queue.add(d);
            }
            //md = queue.pop();
            md = queue.peek();
            volume += 1 / Math.pow(md, dim);
        }
        //double entropy = -Math.log(1 / (volume * (pts.size() - 1) * BALL / pts.size()));
        double entropy = -Math.log(volume / ((pts.size()) * BALL));
        return entropy;
    }

    public double entropySlow(Set<Integer> pts, int k) {
        double volume = .0f;
        for (int me : pts) {
            double md = Double.MAX_VALUE;
            PriorityQueue<Double> queue = new PriorityQueue<Double>(k, new Comparator<Double>() {
                public int compare(Double lhs, Double rhs) {
                    if (lhs > rhs)
                        return 1;
                    if (lhs.equals(rhs))
                        return 0;
                    return -1;
                }
            });
            for (int you : pts) {
                double d = distance(data[me], data[you], isDTW);
                queue.add(d);
            }
            int polled = 0;
            while (!queue.isEmpty() && polled++ < k) {
                md = queue.poll();
            }
            volume += 1 / Math.pow(md, dim);
        }
        //double entropy = -Math.log(1 / (volume * (pts.size() - 1) * BALL / pts.size()));
        double entropy = -Math.log(volume / ((pts.size()) * BALL));
        return entropy;
    }

    // May disable the previous three functions in release.

    private boolean expansible(int begin, int end) {
        double dt = nnMaxavg(end);
        Set<Integer> bs = maxavgNN(begin, dt);
        return (bs.size() > (int) (KT / 2));
    }

    // May be further optimized without computation of begin.
    @SuppressWarnings("unused")
    private boolean reachable(int begin, int end, double radius) {
        Set<Integer> bs = maxavgNN(begin, radius);

        Set<Integer> es = maxavgNN(end, radius);

        bs.retainAll(es);

        return (bs.size() > (int) (KT / 2));
    }

    public List<Integer> getPriority() {
        return priority;
    }

    public double[][] getdists() {
        return dists;
    }

    public Map<Integer, Integer> getIndex() {
        return Index;
    }

    public int getKT() {
        return KT;
    }

    public int getDH() {
        return DH;
    }

    public Map<Integer, Double> getEntropies() {
        return initEntropies;
    }

    private static void exprCircle() {
        long bt = System.currentTimeMillis();
        int cnt = 20;
        int begin = 0;
        double data[][] = new double[4 * cnt + begin][2];
        data[0][0] = .0f;
        data[0][1] = .0f;
        for (int i = 0; i < cnt; i++) {
            data[begin + i + cnt * 0][0] = -1 + (double) i / cnt;
            data[begin + i + cnt * 0][1] = Math.sqrt(1 - (data[begin + i + cnt * 0][0] * data[begin + i + cnt * 0][0]));
            data[begin + i + cnt * 1][0] = 1 - (double) i / cnt;
            data[begin + i + cnt * 1][1] = Math.sqrt(1 - (data[begin + i + cnt * 1][0] * data[begin + i + cnt * 1][0]));
            data[begin + i + cnt * 2][0] = 1 - (double) i / cnt;
            data[begin + i + cnt * 2][1] =
                    -Math.sqrt(1 - (data[begin + i + cnt * 2][0] * data[begin + i + cnt * 2][0]));
            data[begin + i + cnt * 3][0] = -1 + (double) i / cnt;
            data[begin + i + cnt * 3][1] =
                    -Math.sqrt(1 - (data[begin + i + cnt * 3][0] * data[begin + i + cnt * 3][0]));
        }
        long et = System.currentTimeMillis();
        System.out.println("Init: " + (et - bt));
        bt = System.currentTimeMillis();
        ATD atd = new ATD(data, 2, 2, false);
        et = System.currentTimeMillis();
        System.out.println("Construct: " + (et - bt));
        bt = System.currentTimeMillis();
        Set<Integer> set = new HashSet<>();
        for (int i = 0; i < 4 * cnt + begin; i++) {
            set.add(i);
        }
        et = System.currentTimeMillis();
        System.out.println("Set: " + (et - bt));
        bt = System.currentTimeMillis();
        System.out.println(atd.entropy(set, cnt * 4 + begin));
        et = System.currentTimeMillis();
        System.out.println("Entropy: " + (et - bt));
    }

    private static void exprRectangle() {
        long bt = System.currentTimeMillis();
        int card = 60;
        double data[][] = new double[card * card][2];
        for (int i = 0; i < card; i++) {
            for (int j = 0; j < card; j++) {
                data[i * card + j][0] = i;
                data[i * card + j][1] = j;
            }
        }
        long et = System.currentTimeMillis();
        System.out.println("Init: " + (et - bt));
        bt = System.currentTimeMillis();
        ATD atd = new ATD(data, 2, 2, false);
        et = System.currentTimeMillis();
        System.out.println("Construct: " + (et - bt));
        bt = System.currentTimeMillis();
        Set<Integer> set = new HashSet<>();
        for (int i = 0; i < card * card; i++) {
            set.add(i);
        }
        et = System.currentTimeMillis();
        System.out.println("Set: " + (et - bt));
        bt = System.currentTimeMillis();
        System.out.println(atd.entropy(set, 40));
        et = System.currentTimeMillis();
        System.out.println("Entropy: " + (et - bt));
    }

    public static void main(String[] args) {
        exprCircle();
        exprRectangle();
    }
}
