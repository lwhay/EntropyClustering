package ics.whu.edu.cn.madrix.clustering.density;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections15.list.SetUniqueList;
import org.apache.commons.math3.special.Gamma;

public class CTD extends AbstractDensityClustering {

    private double[][] dists;
    private static Map<Integer, List<Integer>> distSortedNbs = new HashMap<>();
    private Map<Integer, Double> initEntropies = new HashMap<>();
    @SuppressWarnings("unused")
    private double[][] data;
    private int KKm;
    private int KT;
    private int DH;
    private Map<Integer, Integer> Index = new HashMap<>();
    private Map<Integer, Integer> Index0 = new HashMap<>();
    private Map<Integer, Double> rhoRange = new HashMap<Integer, Double>();
    private List<Integer> priority = new ArrayList<Integer>();
    private int dim;
    private double BALL;
    private final ATD atd;

    public CTD(double[][] data, ATD atd, List<Integer> rhoK, int K) {
        this.atd = atd;
        this.data = data;
        this.dim = data[0].length;
        this.BALL = Math.pow(Math.PI, dim / 2) / Gamma.gamma(dim / 2 + 1);
        priority = atd.getPriority();
        initEntropies = atd.getEntropies();
        KT = atd.getKT();
        DH = KT;
        //Specially, we can use KKm = KT + 2 in addition to ATD's commented lines to boost D31 and beyond.
        KKm = 2 * KT + 1;
        Map<Integer, Integer> NX = CBridgeKernels(atd, rhoK, K);
        Index = NX;
        Collection<Integer> c = NX.values();
        Object[] cindex = c.toArray();
        for (Object nid : cindex) {
            System.out.print((Integer) nid + 1 + "\t");
        }
        System.out.print("\n");
    }

    public Map<Integer, Integer> CBridgeKernels(ATD atd, List<Integer> rhoK, int K) {
        // Bridge kernels
        Index = atd.getIndex();
        dists = atd.getdists();
        Map<Integer, Double> initEntropies = atd.getEntropies();
        distSortedNbs = atd.getNbs();
        List<Integer> rhoH = SetUniqueList.decorate(rhoK);
        List<Integer> rhom = new ArrayList<Integer>();
        for (int id : priority) {
            if (rhoH.contains(id))
                rhom.add(id);
        }
        Map<Integer, Integer> cindex = new HashMap<Integer, Integer>();
        for (int id : rhom) {
            cindex.put(id, atd.getKeyByValue(Index, id).size());
        }
        int Ko = (int) Math.pow(KKm, KT) < Index.size() - 1 ? (int) Math.pow(KKm, KT) : Index.size() - 1;
        for (int i = 0; i < rhom.size(); i++) {
            if (cindex.get(rhom.get(i)) <= KT) {
                int NNindex = 1;
                int nns = 1;
                boolean nno = false;
                int NNcandidate = distSortedNbs.get(rhom.get(i)).get(NNindex);
                List<Integer> tc = atd.getKeyByValue(Index, rhom.get(i));
                boolean lo = false;
                while (NNindex < Ko) {
                    if (expansible(rhom.get(i), NNcandidate) && Isexpansible(tc, NNcandidate)
                            && (tc.contains(NNcandidate)) == false) {// && reachable(rhom.get(i), NNcandidate, nnMaxavg(rhom.get(i)))) {          //initEntropies.get(NNcandidate) <= initEntropies.get(rhom.get(i)) &&     // && expansible(NNcandidate,rhom.get(i))           
                        lo = true;
                        break;
                    } else {
                        if (Isexpansible(tc, NNcandidate) && nno == false) {
                            nns = NNindex;
                            nno = true;
                        }
                        NNindex = NNindex + 1;
                        NNcandidate = distSortedNbs.get(rhom.get(i)).get(NNindex);
                    }
                }
                if (lo == false)
                    NNcandidate = distSortedNbs.get(rhom.get(i)).get(nns);
                for (int id : tc) {
                    Index.put(id, Index.get(NNcandidate));
                }
                cindex.put(rhom.get(i), cindex.get(rhom.get(i)) + tc.size()); //update
            }
        }
        Collection<Integer> c = Index.values();
        Object[] tmp = c.toArray();
        List<Integer> rhoT = new ArrayList<Integer>();
        for (Object nid : tmp) {
            rhoT.add((Integer) nid);
            System.out.print((Integer) nid + 1 + "\t");
        }
        System.out.print("\n");

        Index0 = Index;
        //Bridge kernels
        List<Integer> rhoms = SetUniqueList.decorate(rhoT);
        List<Integer> CNrho = new ArrayList<Integer>();
        List<Integer> BDrho = new ArrayList<Integer>();
        Map<Integer, Double> CNOrder = new HashMap<Integer, Double>();
        Map<Integer, Double> BridgeOrder = new HashMap<Integer, Double>();
        for (int od : rhoms) {
            List<Integer> tc = atd.getKeyByValue(Index, od);
            int Km = 0;
            for (int ii = 0; ii <= KT; ii++) {
                if (tc.contains(distSortedNbs.get(od).get(ii)))
                    Km = Km + 1;
            }
            if (Km == KT + 1 && tc.size() >= KKm) {
                CNrho.add(od);
                CNOrder.put(od, initEntropies.get(od));
            } else {
                BDrho.add(od);
                BridgeOrder.put(od, initEntropies.get(od));
            }
        }
        /*Map<Integer, Double> CNOrders = sortMapByValue(CNOrder);//sortMapByValue
        Map<Integer, Integer> Index0 = new HashMap<Integer, Integer>();
        if (CNOrder.size() > K * (3 * dim)) {
            int count = CNrho.size();
            List<Integer> CNkernels = new ArrayList<Integer>();
            while (count > K * (3 * dim)) {
        Index0 = Index;
        CNkernels = new ArrayList<Integer>();
        Iterator<Integer> iter = CNOrders.keySet().iterator();
        while (iter.hasNext()) {
            Integer i = iter.next();
            List<Integer> tc = atd.getKeyByValue(Index0, i);
            Integer newindex = ComputediffrhoM(i, Index0, CNrho);
            List<Integer> ts = atd.getKeyByValue(Index0, newindex);
            if (expansA(tc, ts)) { //2*(KT + 1)                   
                for (int ii = 0; ii < tc.size(); ii++) {
                    Index0.put(tc.get(ii), newindex);
                }
            } else
                CNkernels.add(i);
        }
        if (count == CNkernels.size())
            break;
        else
            count = CNkernels.size();
        KKm = KKm + 1;
            }
        
            Index = Index0;
            KKm = KKm - 1;
            CNrho = CNkernels;
        }*/

        c = Index.values();
        tmp = c.toArray();
        rhoT = new ArrayList<Integer>();
        for (Object nid : tmp) {
            rhoT.add((Integer) nid);
            System.out.print((Integer) nid + 1 + "\t");
        }
        System.out.print("\n");
        Map<Integer, Double> BridgeOrders = sortMapByValue(BridgeOrder);
        MergeBridgeKernels(BridgeOrders, BDrho, CNrho, K + 1, 1, true);

        c = Index.values();
        tmp = c.toArray();
        rhoT = new ArrayList<Integer>();
        for (Object nid : tmp) {
            rhoT.add((Integer) nid);
            System.out.print((Integer) nid + 1 + "\t");
        }
        System.out.print("\n");
        rhoms = SetUniqueList.decorate(rhoT);

        BridgeOrder = new HashMap<Integer, Double>();
        for (int od : rhoms) {
            BridgeOrder.put(od, initEntropies.get(od));
        }
        BridgeOrders = sortMapByValue(BridgeOrder);
        Iterator<Integer> iter = BridgeOrders.keySet().iterator();
        while (iter.hasNext()) {
            int i = iter.next();
            /*rhoRange.put(i, initEntropies.get(i));  //1
            List<Integer> tx = new ArrayList<Integer>(); //2
            for(int ii = 0; ii < KKm ; ii++) { //KKm // 2*(KT + 1)  * DH
                tx.add(distSortedNbs.get(i).get(ii));
            }
            rhoRange.put(i, entropy(tx,tx.size()));*/

            List<Integer> tc = atd.getKeyByValue(Index0, i); // 2s    
            rhoRange.put(i, getRentropy(tc));

            /*List<Integer> tc = atd.getKeyByValue(Index, i);
            Set<Integer> tx = new HashSet<Integer>(); //3
            for (int ii = 0; ii < tc.size(); ii++) {
                if (tc.contains(distSortedNbs.get(i).get(ii)))
                    tx.add(distSortedNbs.get(i).get(ii));
                if (tx.size() >= KKm * DH)
                    break;
            }
            rhoRange.put(i, entropy(tx, tx.size()));*/
        }

        return Index;
    } // rank index, push into the stack

    private Double getRentropy(List<Integer> tx) {
        int top = 0; // find the least entropy index
        double maxentropy = Double.MAX_VALUE;
        for (int id : tx) {
            if (initEntropies.get(id) < maxentropy) {
                maxentropy = initEntropies.get(id);
                top = id;
            }
        }
        List<Integer> ts = new ArrayList<Integer>(); //2
        for (int ii = 0; ii < KKm; ii++) { //KKm // 2*(KT + 1)  * DH
            ts.add(distSortedNbs.get(top).get(ii));
        }
        return entropy(ts, ts.size());
    }

    private void MergeBridgeKernels(Map<Integer, Double> BridgeOrders, List<Integer> BDrho, List<Integer> CNrho, int nn,
                                    int m, boolean b) {
        int count = 0;
        if (BDrho.size() > 0) {
            Map<Integer, Double> Bridgenew = new HashMap<Integer, Double>();
            if (b) {
                Iterator<Integer> iter = BridgeOrders.keySet().iterator();
                while (iter.hasNext()) {
                    Integer i = iter.next();
                    List<Integer> tc = atd.getKeyByValue(Index, i);
                    Integer newindex = findBKrelayCN(i, CNrho, nn, m);
                    if (newindex != null) {
                        for (int ii = 0; ii < tc.size(); ii++) {
                            Index.put(tc.get(ii), newindex);
                        }
                        count = count + 1;
                        BDrho.remove(i);
                    } else
                        Bridgenew.put(i, BridgeOrders.get(i));
                }
                Collection<Integer> c = Index.values();
                Object[] tmp = c.toArray();
                List<Integer> rhoT = new ArrayList<Integer>();
                for (Object nid : tmp) {
                    rhoT.add((Integer) nid);
                    System.out.print((Integer) nid + 1 + "\t");
                }
                System.out.print("\n");
            } else {
                // can not merge !!!!
                Iterator<Integer> iter = BridgeOrders.keySet().iterator();
                while (iter.hasNext()) {
                    Integer i = iter.next();
                    List<Integer> tc = atd.getKeyByValue(Index, i);
                    //                    Integer newindex = ComputeClassdist(i, CNrho);
                    Integer newindex = ComputeMinrhorelay(i, CNrho, nn + 1);
                    for (int ii = 0; ii < tc.size(); ii++) {
                        Index.put(tc.get(ii), newindex);
                    }
                    BDrho.remove(i);
                }
            }

            if (count == 0)
                b = false;
            MergeBridgeKernels(Bridgenew, BDrho, CNrho, nn, m + 1, b);
        }
        return;
    }

    public Integer ComputeMinrhorelay(Integer i, List<Integer> CNrho, int nn) {
        Integer newindex = null;
        List<Integer> tc = atd.getKeyByValue(Index, i);
        List<Integer> CNs = new ArrayList<Integer>();
        CNs.addAll(CNrho);
        while (true) {
            while (CNs.size() > 0) {
                Integer tmp = ComputediffrhoM(i, Index0, CNs);
                List<Integer> ts = atd.getKeyByValue(Index0, tmp);
                if (expansA(tc, ts, nn)) {
                    newindex = tmp;
                    break;
                } else
                    CNs.remove(tmp);
            }
            if (nn == KKm) {
                newindex = ComputediffrhoM(i, Index0, CNrho);
                break;
            }

            if (newindex != null)
                break;
            else
                nn = nn + 1;
        }
        return newindex;
    }

    private Integer findBKrelayCN(Integer i, List<Integer> CNrho, int nn, int m) {
        Integer newindex = null;
        List<Integer> tc = atd.getKeyByValue(Index, i);
        List<Integer> CNs = new ArrayList<Integer>();
        CNs.addAll(CNrho);
        while (m > 0) {
            Integer tmp = ComputediffrhoM(i, Index0, CNs);
            List<Integer> ts = atd.getKeyByValue(Index0, tmp);
            if (expansA(tc, ts, nn)) {
                newindex = tmp;
                break;
            } else
                CNs.remove(tmp);
            m = m - 1;
        }
        return newindex;
    }

    public Map<Integer, Double> sortMapByValue(Map<Integer, Double> oriMap) {
        if (oriMap == null || oriMap.isEmpty()) {
            return null;
        }
        Map<Integer, Double> sortedMap = new LinkedHashMap<Integer, Double>();
        List<Map.Entry<Integer, Double>> entryList = new ArrayList<Map.Entry<Integer, Double>>(oriMap.entrySet());
        Collections.sort(entryList, new MapValueComparator());

        Iterator<Map.Entry<Integer, Double>> iter = entryList.iterator();
        Map.Entry<Integer, Double> tmpEntry = null;
        while (iter.hasNext()) {
            tmpEntry = iter.next();
            sortedMap.put(tmpEntry.getKey(), tmpEntry.getValue());
        }
        return sortedMap;
    }

    private boolean Isexpansible(List<Integer> tc, int NN) {
        boolean re = false;
        for (int id : tc) {
            for (int i = 0; i < KKm; i++) { //KKm
                if (distSortedNbs.get(id).get(i).equals(NN)) {
                    re = true;
                    break;
                }
            }
            for (int i = 0; i < KKm; i++) {
                if (distSortedNbs.get(NN).get(i).equals(id)) {
                    re = true;
                    break;
                }
            }
        }
        return re;
    }

    public boolean expansA(List<Integer> t1, List<Integer> t2, int nn) {
        Set<Integer> ba = maxavgNNA(t1, nn);
        Set<Integer> bs = maxavgNNA(t2, nn);
        boolean re = true;
        Set<Integer> res = new HashSet<>();
        res.clear();
        res.addAll(ba);
        res.retainAll(bs);
        if (res.size() < 1)//(int) KT / 2)
            re = false;
        return re;
    }

    private Set<Integer> maxavgNNA(List<Integer> nns, int nn) {
        Set<Integer> bs = new HashSet<Integer>();
        for (Integer id : nns) {
            Set<Integer> ba = getKTNN(id, nn);
            for (int ob : ba) {
                if (expansible(id, ob))
                    bs.add(ob);
            }
        }
        return bs;
    }

    private Set<Integer> getKTNN(int center, int nn) {
        Set<Integer> nns = new HashSet<>();
        for (int i = 0; i < nn; i++) { //(1 + DH )* 
            if (expansible(center, distSortedNbs.get(center).get(i)))
                nns.add(distSortedNbs.get(center).get(i));
        }
        return nns;
    }

    public Integer ComputediffrhoM(Integer i, Map<Integer, Integer> Index, List<Integer> rhoms) {
        List<Integer> tmp = new ArrayList<Integer>();
        tmp.addAll(rhoms);
        tmp.remove(i);
        //        double[] distC = new double[rhoms.size()-2];
        Map<Integer, Double> distC = new HashMap<Integer, Double>();
        int f = 0;
        while (f < tmp.size()) {
            //  the nearest yellow cluster
            List<Integer> x1 = atd.getKeyByValue(Index, i);
            List<Integer> x2 = atd.getKeyByValue(Index, tmp.get(f));
            List<Integer> ts2 = new ArrayList<Integer>();
            ts2.addAll(x2);
            ts2.retainAll(rhoms); // intersection
            if (ts2.contains(tmp.get(f)) == false)
                ts2.add(tmp.get(f));
            int end = ComputeClassdist(i, ts2);

            x2 = atd.getKeyByValue(Index, end);

            int sm = x1.size() < x2.size() ? x1.size() : x2.size();
            //            sm = KKm > sm ? KKm : sm;

            if (sm == x1.size()) {
                x2 = reBuild(x2, tmp.get(f), sm);
            } else {
                x1 = reBuild(x1, i, sm);
            }

            double dd1 = ComputeInclassdist(x1);

            double dd2 = ComputeInclassdist(x2);

            double dx1 = ComputeBeclassdist(x1, x2);

            double dx2 = ComputeBeclassdist(x2, x1);

            //double fz = dx1 > dx2 ? dx1 : dx2; //(dx1 + dx2) / 2;//

            distC.put(tmp.get(f), (dx1 * dx2) / (dd1 * dd2)); //(fz) / Math.sqrt(dd1*dd2)) (dx1*dx2) / (dd1*dd2)
            f = f + 1;
        }
        return getMinValue(distC);
    }

    public int ComputeClassdist(int ix, List<Integer> list) {
        int rely = 0;
        double distC = Double.MAX_VALUE;
        List<Integer> x1 = atd.getKeyByValue(Index, ix);
        for (int i = 0; i < list.size(); i++) {
            double mindist = 0; // Accumulate
            List<Integer> x2 = atd.getKeyByValue(Index, list.get(i));
            for (int od : x1) {
                double mm = Double.MAX_VALUE;
                for (int sd : x2) {
                    if (dists[od][sd] < mm)
                        mm = dists[od][sd];
                }
                mindist = mindist + mm;
            }

            for (int od : x2) {
                double mm = Double.MAX_VALUE;
                for (int sd : x1) {
                    if (dists[od][sd] < mm)
                        mm = dists[od][sd];
                }
                mindist = mindist + mm;
            }

            if (distC > (mindist / (x1.size() + x2.size())) && mindist / (x1.size() + x2.size()) > 0) { //distC mean
                distC = (mindist / (x1.size() + x2.size()));
                rely = list.get(i);
            }
        }
        return rely;
    }

    public List<Integer> reBuild(List<Integer> x2, Integer integer, int sm) {
        List<Integer> re = new ArrayList<Integer>();
        for (int id : distSortedNbs.get(integer)) {
            if (x2.contains(id))
                re.add(id);
            if (re.size() == sm)
                break;
        }
        return re;
    }

    public Integer getMinValue(Map<Integer, Double> map) {
        if (map == null)
            return null;
        Collection<Double> c = map.values();
        Object[] obj = c.toArray();
        Arrays.sort(obj);
        int reIndex = 0;
        for (Integer getKey : map.keySet()) {
            if (map.get(getKey).equals(obj[0])) {
                reIndex = getKey;
                break;
            }
        }
        return reIndex;
    }

    public double ComputeBeclassdist(List<Integer> x1, List<Integer> x2) {
        double ds1 = 0;
        for (int i = 0; i < x1.size(); i++) {
            double minds1 = Double.MAX_VALUE;
            for (int j = 0; j < x2.size(); j++) {
                double dist = dists[x1.get(i)][x2.get(j)];
                if (minds1 > dist && dist > 0)
                    minds1 = dist;
            }
            ds1 = ds1 + minds1;
        }
        return ds1 / x1.size();
    }

    public double ComputeInclassdist(List<Integer> x1) {
        double dd1 = 0;
        for (int i = 0; i < x1.size(); i++) {
            double mindd1 = Double.MAX_VALUE;
            for (int j = 0; j < x1.size(); j++) {
                double dist = dists[x1.get(i)][x1.get(j)];
                if (mindd1 > dist && dist > 0)
                    mindd1 = dist;
            }
            dd1 = dd1 + mindd1;
        }
        return dd1 / x1.size();
    }

    public boolean expansible(int begin, int end) {
        double dt = nnMaxavg(end);
        Set<Integer> bs = maxavgNN(begin, dt);
        return (bs.size() > (int) (KT / 2));
    }

    @SuppressWarnings("unused")
    private boolean reachable(int begin, int end, double radius) {
        Set<Integer> bs = maxavgNN(begin, radius);

        Set<Integer> es = maxavgNN(end, radius);

        bs.retainAll(es);

        return (bs.size() > (int) (KT / 2));
    }

    private double nnMaxavg(int center) {
        Set<Integer> nns = new HashSet<>();
        for (int i = 0; i < KT + 1; i++) {
            nns.add(distSortedNbs.get(center).get(i));
        }
        return maxAvgDist(nns);
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

    public double entropy(List<Integer> pts, int k) {
        double volume = .0f;
        for (int me : pts) {
            double md = Double.MAX_VALUE;
            FastBoundedPriorityQueue<Double> queue = new FastBoundedPriorityQueue<>(k, true);
            for (int you : pts) {
                double d = dists[me][you];
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

    public Map<Integer, Integer> getIndex() {
        return Index;
    }

    public Map<Integer, Double> getrhoRange() {
        return rhoRange;
    }

    public double[][] getdists() {
        return dists;
    }

    public int getKT() {
        return KT;
    }

    public int getKKm() {
        return KKm;
    }

    public List<Integer> getPriority() {
        return priority;
    }

    public Map<Integer, List<Integer>> getNbs() {
        return distSortedNbs;
    }

    public int getDH() {
        return DH;
    }

    public Map<Integer, Integer> getIndex0() {
        return Index0;
    }

}
