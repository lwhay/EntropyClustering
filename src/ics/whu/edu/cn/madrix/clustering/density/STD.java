package ics.whu.edu.cn.madrix.clustering.density;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections15.list.SetUniqueList;

public class STD {
    private Map<Integer, List<Integer>> distSortedNbs = new HashMap<>();
    private Map<Integer, Double> reStack = new HashMap<Integer, Double>();
    private Map<Integer, List<Integer>> rePath = new HashMap<Integer, List<Integer>>();
    private Map<Integer, Integer> Index = new HashMap<>();
    private Map<Integer, Integer> Index0 = new HashMap<>();
    //private int KT;
    private List<Integer> currentcenters = new ArrayList<Integer>();
    private List<Integer> rhoT = new ArrayList<Integer>();
    private List<Integer> rhoms = new ArrayList<Integer>();
    private static double[][] dists;
    private final CTD ctd;
    private final ATD atd;

    public STD(CTD ctd, ATD atd, int K, int DH) {
        this.ctd = ctd;
        this.atd = atd;
        dists = ctd.getdists();
        Index = ctd.getIndex();
        distSortedNbs = ctd.getNbs();
        Index0 = Index;
        Collection<Integer> c = Index0.values();
        Object[] tmp = c.toArray();
        rhoT = new ArrayList<Integer>();
        for (Object nid : tmp) {
            rhoT.add((Integer) nid);
        }
        rhoms = SetUniqueList.decorate(rhoT);
        Map<Integer, Double> rhoRange = ctd.sortMapByValue(ctd.getrhoRange());
        Iterator<Integer> iter = rhoRange.keySet().iterator();
        // the Global Max-- point 
        int s0 = iter.next();
        currentcenters.add(s0);
        List<Integer> tmppath = new ArrayList<Integer>();
        tmppath.add(s0);
        rePath.put(s0, tmppath);
        reStack.put(s0, Double.MAX_VALUE);
        Index.put(s0, s0);

        // 2th -- k
        int count = 1;
        while (iter.hasNext()) {
            int ix = iter.next();
            int l0 = ComputeClassdist(ix, currentcenters);
            List<Integer> path0 = new ArrayList<Integer>();
            path0.add(l0);
            rePath.put(ix, path0);
            reStack.put(ix, ComputeMindist(ix, l0));
            currentcenters.add(ix);
            count = count + 1;
            if (count == K)
                break;
        }

        // k --
        while (iter.hasNext()) {
            int ix = iter.next();
            int l1 = ComputeClassdist(ix, currentcenters);
            double currentdist = ComputeMindist(ix, l1);
            if (Ifpush(reStack.values(), currentdist)) {
                List<Integer> path1 = new ArrayList<Integer>();
                path1.add(l1);
                rePath.put(ix, path1);
                reStack.put(ix, currentdist);
                currentcenters.add(ix);
                int start = 0;
                for (int m = 0; m < currentcenters.size(); m++) {
                    if (currentcenters.get(m) == l1) {
                        start = m;
                        break;
                    }
                }
                for (int i = start + 1; i < currentcenters.size(); i++) {
                    List<Integer> newcenters = new ArrayList<Integer>();
                    for (int m = 0; m < currentcenters.size(); m++) {
                        boolean bo = ifcenterp(rePath, currentcenters, m, i);
                        if (m != i && bo) {
                            newcenters.add(currentcenters.get(m));
                        }
                    }
                    if (newcenters.size() > 0) {
                        int l2 = ComputeClassdist(currentcenters.get(i), newcenters);
                        List<Integer> path2 = rePath.get(currentcenters.get(i));
                        if (path2.get(path2.size() - 1) != l2) {
                            path2.add(l2);
                            int l2s = ComputeClassdist(currentcenters.get(i), path2);
                            path2.remove(path2.get(path2.size() - 1));
                            path2.add(l2s);
                        }
                        rePath.put(currentcenters.get(i), path2); // ????
                        reStack.put(currentcenters.get(i), ComputeMindist(currentcenters.get(i), l2));
                    }
                }
            } else {
                UpdateIndex(ix, l1);
                for (int i = 1; i < currentcenters.size(); i++) {
                    int curId = currentcenters.get(i);
                    List<Integer> path2 = rePath.get(curId);
                    if (path2.contains(l1)) {
                        List<Integer> tmpcenters = currentcenters;
                        for (int ii = 0; ii < tmpcenters.size(); ii++) {
                            if (tmpcenters.get(ii) == curId) {
                                tmpcenters.remove(tmpcenters.get(ii));
                                break;
                            }
                        }
                        while (true) {
                            int s1 = ComputeClassdist(curId, tmpcenters);
                            int addr = 0;
                            for (int ii = 0; ii < currentcenters.size(); ii++) {
                                if (currentcenters.get(ii) == s1) {
                                    addr = ii;
                                    break;
                                }
                            }
                            boolean bo = ifcenterp(rePath, currentcenters, addr, i);
                            if (bo == false) {
                                for (int ii = 0; ii < tmpcenters.size(); ii++) {
                                    if (tmpcenters.get(ii) == s1) {
                                        tmpcenters.remove(tmpcenters.get(ii));
                                        break;
                                    }
                                }
                            } else {
                                double xp = ComputeMindist(curId, s1);
                                if (xp < reStack.get(curId)) {
                                    reStack.put(curId, xp);
                                    if (path2.get(path2.size() - 1) != s1)
                                        path2.add(s1);
                                    rePath.put(curId, path2);
                                }
                                break;
                            }
                        }
                    }
                }
            }

            count = count + 1;
        }
    }

    private void UpdateIndex(int l0, int l1) {
        List<Integer> t1 = atd.getKeyByValue(Index, l0);
        List<Integer> ts = atd.getKeyByValue(Index, l1);
        if (expansA(t1, ts)) {
            for (int id : t1)
                Index.put(id, l1);
        } else {
            List<Integer> path1 = new ArrayList<Integer>();
            double currentdist = ComputeMindist(l0, l1);
            path1.add(l1);
            rePath.put(l0, path1);
            reStack.put(l0, currentdist);
            currentcenters.add(l0);
        }

    }

    public static boolean expansA(List<Integer> t1, List<Integer> t2) {
        boolean re = false;
        double dd1 = Double.MAX_VALUE;
        for (int i = 0; i < t1.size(); i++) {
            for (int j = 0; j < t2.size(); j++) {
                double dist = dists[t1.get(i)][t2.get(j)];
                if (dd1 > dist && dist > 0)
                    dd1 = dist;
            }
        }
        if (dd1 < nnMaxavg(t2))
            re = true;
        return re;
    }

    private static double nnMaxavg(List<Integer> t2) {
        double dd1 = Double.MAX_VALUE;
        for (int i = 0; i < t2.size(); i++) {
            for (int j = i + 1; j < t2.size(); j++) {
                double dist = dists[t2.get(i)][t2.get(j)];
                if (dd1 > dist && dist > 0)
                    dd1 = dist;
            }
        }
        return dd1;
    }

    private boolean ifcenterp(Map<Integer, List<Integer>> rePath, List<Integer> cur, int m, int i) {
        boolean re = true;
        m = cur.get(m);
        while (m != cur.get(0) && i < cur.size()) {
            if (rePath.get(m).contains(cur.get(i))) {
                re = false;
                break;
            } else {
                m = rePath.get(m).get(0); // search history
            }
        }
        return re;
    }

    @SuppressWarnings("unused")
    private Map<Integer, List<Integer>> getCenters(Map<Integer, List<Integer>> rePath2, Integer integer) {
        Iterator<Integer> iter = rePath2.keySet().iterator();
        List<Integer> re = new ArrayList<Integer>();
        Map<Integer, List<Integer>> rre = new HashMap<Integer, List<Integer>>();
        while (iter.hasNext()) {
            Integer ix = iter.next();
            List<Integer> le = rePath2.get(ix);
            re.add(le.get(le.size() - 1));
        }
        int addr = 0;
        for (int i = 0; i < re.size(); i++) {
            if (re.get(i) == integer) {
                addr = i;
                break;
            }
        }
        rre.put(addr, re);
        return rre;
    }

    private boolean Ifpush(Collection<Double> values, double currentdist) {
        boolean re = false;
        Object[] cindex = values.toArray();
        for (Object nid : cindex) {
            if ((double) nid < currentdist) {
                re = true;
                break;
            }
        }
        return re;
    }

    public Double ComputeMindist(int ix, int s1) {
        List<Integer> t1 = atd.getKeyByValue(Index, s1);
        List<Integer> ts1 = t1;
        ts1.retainAll(rhoms); // intersection
        int end = ComputeClassdist(ix, ts1);

        List<Integer> t2 = atd.getKeyByValue(Index, ix);
        List<Integer> ts2 = t2;
        ts2.retainAll(rhoms); // intersection
        int start = ComputeClassdist(end, ts2);

        double dist = ComputerhoKdist(start, end); //start , end
        return dist;
    }

    private double ComputerhoKdist(int start, int end) {
        List<Integer> xs1 = atd.getKeyByValue(Index0, start);
        List<Integer> xs2 = atd.getKeyByValue(Index0, end);

        int Num = xs1.size() < xs2.size() ? xs1.size() : xs2.size();//(4*(KT+1)); 
        List<Integer> x1 = new ArrayList<Integer>();
        int count = 0;
        for (int i = 0; i < xs1.size(); i++) {
            if (xs1.contains(distSortedNbs.get(start).get(i))) {
                x1.add(distSortedNbs.get(start).get(i));
                if (count < Num) //Math.floor(KKm / 2))
                    count = count + 1;
                else
                    break;
            }
        }
        List<Integer> x2 = new ArrayList<Integer>();
        count = 0;
        for (int i = 0; i < xs2.size(); i++) {
            if (xs2.contains(distSortedNbs.get(end).get(i))) {
                x2.add(distSortedNbs.get(end).get(i));
                if (count < Num)
                    count = count + 1;
                else
                    break;
            }
        }
        double dd1 = ctd.ComputeInclassdist(x1);

        double dd2 = ctd.ComputeInclassdist(x2);

        double dx1 = ctd.ComputeBeclassdist(x1, x2);

        double dx2 = ctd.ComputeBeclassdist(x2, x1);

        //double fz = dx1 > dx2 ? dx1 : dx2; //(dx1 + dx2) / 2;//

        double dist = (dx1 * dx2) / (dd1 * dd2); //(dx1*dx2) / (dd1*dd2); (fz) / Math.sqrt(dd1*dd2)

        return dist;
    }

    public int ComputeClassdist(Integer ix, List<Integer> list) {
        int rely = 0;
        List<Integer> tmp = new ArrayList<Integer>();
        tmp.addAll(list);
        tmp.remove(ix);
        double distC = Double.MAX_VALUE;
        List<Integer> x1 = atd.getKeyByValue(Index, ix);
        for (int i = 0; i < tmp.size(); i++) {
            double mindist = 0; // Accumulate
            List<Integer> x2 = atd.getKeyByValue(Index, tmp.get(i));
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

            if (distC > (mindist / (x1.size() + x2.size()))) { //distC mean
                distC = (mindist / (x1.size() + x2.size()));
                rely = tmp.get(i);
            }
        }
        return rely;
    }

    public Map<Integer, Double> getreStack() {
        return reStack;
    }

    public Map<Integer, Integer> getIndex() {
        return Index;
    }

    public Map<Integer, List<Integer>> getrePath() {
        return rePath;
    }

    public List<Integer> getrhoms() {
        return rhoms;
    }

}
