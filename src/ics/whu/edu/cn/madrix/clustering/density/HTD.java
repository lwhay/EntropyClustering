package ics.whu.edu.cn.madrix.clustering.density;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class HTD {

    private double[][] dists;
    private Map<Integer, Integer> Index = new HashMap<>();
    private Map<Integer, Integer> GKer = new LinkedHashMap<Integer, Integer>();
    private List<Integer> centers = new ArrayList<Integer>();
    private Map<Integer, List<Integer>> rePath = new HashMap<Integer, List<Integer>>();
    private final CTD ctd;
    private final ATD atd;

    public HTD(CTD ctd, STD std, ATD atd, Map<Integer, Double> resultK, List<Integer> GKernels,
            Map<Integer, List<Integer>> derivedPath) {
        this.ctd = ctd;
        this.atd = atd;
        List<Integer> priority = ctd.getPriority();
        Collections.reverse(priority);
        dists = ctd.getdists();
        rePath = std.getrePath();

        // stack index
        centers = getKeysinMap(resultK);
        Index = std.getIndex();

        for (int id : priority) {
            if (GKernels.contains(id))
                GKer.put(id, GKer.size());
            if (GKer.size() == GKernels.size())
                break;
        } // ranking kernels  

        Map<Integer, Integer> label = new LinkedHashMap<Integer, Integer>();
        for (Integer id : GKernels) {
            centers.remove(id);
            List<Integer> tc = atd.getKeyByValue(Index, id);
            for (int nd : tc)
                label.put(nd, GKer.get(id));
        }

        Map<Integer, Integer> cerl = new LinkedHashMap<Integer, Integer>();
        for (int id : priority) {
            if (centers.contains(id))
                cerl.put(id, 0);
            if (cerl.size() == centers.size())
                break;
        }

        //  after merge bridge kernels --index

        Collection<Integer> c = Index.values();
        Object[] tmp = c.toArray();
        List<Integer> rhoT = new ArrayList<Integer>();
        for (Object nid : tmp) {
            rhoT.add((Integer) nid);
        }
        //List<Integer> rhoms = SetUniqueList.decorate(rhoT);

        Iterator<Integer> iter = cerl.keySet().iterator(); // GKer  
        List<Integer> ls = new ArrayList<Integer>();
        List<Integer> ns = new ArrayList<Integer>();
        while (iter.hasNext()) {
            int i = iter.next();
            if (ns.contains(i) == false) {
                int is = rePath.get(i).get(rePath.get(i).size() - 1);
                if (GKernels.contains(is) == false) {
                    Integer str = ctd.ComputediffrhoM(i, Index, GKernels);
                    Integer end = ctd.ComputediffrhoM(is, Index, GKernels);
                    if (str.equals(end)) {
                        List<Integer> tc = atd.getKeyByValue(Index, i);
                        for (int id : tc) {
                            Index.put(id, str);
                            label.put(id, GKer.get(str));
                        }
                        tc = atd.getKeyByValue(Index, is);
                        for (int id : tc) {
                            Index.put(id, str);
                            label.put(id, GKer.get(str));
                        }
                        ns.add(i);
                        ns.add(is);
                    } else {
                        ls.add(i);
                    }
                } else {
                    ls.add(i);
                }
            }
        }

        //        Iterator<Integer> iter = GKer.keySet().iterator();   // GKer
        //        int count = 0;
        //        int i = iter.next();
        //        while(count < GKer.size()) {            
        //            Integer str = ComputediffrhoM(i,Index,centers);
        //            Integer end = ComputediffrhoM(str,Index,GKernels);
        //            if (end.equals(i)) {
        //                List<Integer> tc = ATD.getKeyByValue(Index , str);
        //                for (int id : tc) {
        //                    Index.put(id, i);
        //                    label.put(id, GKer.get(i));
        //                }
        //                centers.remove(str);
        //            }else {
        //                count = count + 1;
        //                i = iter.next();
        //            }
        //        }

        Collections.reverse(ls);
        for (int id : ls) {
            Integer str = ComputediffrhoM(id, Index, GKernels);
            List<Integer> tc = atd.getKeyByValue(Index, str);
            for (int od : tc) {
                Index.put(od, id);
                label.put(od, GKer.get(id));
            }
        }

        Collection<Integer> cs = Index.values();
        Object[] csindex = cs.toArray();
        List<Integer> rhosK = new ArrayList<Integer>();
        for (Object nid : csindex) {
            rhosK.add((Integer) nid);
            System.out.print((Integer) nid + 1 + "\t");
        }
        System.out.print("\n");

    }

    public Integer ComputediffrhoM(int i, Map<Integer, Integer> Index, List<Integer> tmp) {
        List<Integer> x1 = atd.getKeyByValue(Index, i);
        Map<Integer, Double> distC = new HashMap<Integer, Double>();
        int f = 0;
        while (f < tmp.size()) {

            List<Integer> x2 = atd.getKeyByValue(Index, tmp.get(f));

            double dd1 = ComputeInclassdist(x1);

            double dd2 = ComputeInclassdist(x2);

            double dx1 = ComputeBeclassdist(x1, x2);

            double dx2 = ComputeBeclassdist(x2, x1);

            //double fz = dx1 > dx2 ? dx1 : dx2; //(dx1 + dx2) / 2;//

            distC.put(tmp.get(f), (dx1 * dx2) / (dd1 * dd2)); //(fz) / Math.sqrt(dd1*dd2)) (dx1*dx2) / (dd1*dd2)
            f = f + 1;
        }
        return ctd.getMinValue(distC);
    }

    private double ComputeBeclassdist(List<Integer> x1, List<Integer> x2) {
        double ds1 = 0;
        for (int i = 0; i < x1.size(); i++) {
            double maxds1 = 0;
            for (int j = 0; j < x2.size(); j++) {
                double dist = dists[x1.get(i)][x2.get(j)];
                if (maxds1 < dist && dist > 0)
                    maxds1 = dist;
            }
            ds1 = ds1 + maxds1;
        }
        return ds1 / x1.size();
    }

    private double ComputeInclassdist(List<Integer> x1) {
        double dd1 = 0;
        for (int i = 0; i < x1.size(); i++) {
            double maxdd1 = 0;
            for (int j = 0; j < x1.size(); j++) {
                double dist = dists[x1.get(i)][x1.get(j)];
                if (maxdd1 < dist && dist > 0)
                    maxdd1 = dist;
            }
            dd1 = dd1 + maxdd1;
        }
        return dd1 / x1.size();
    }

    @SuppressWarnings("unused")
    private List<Integer> Rankintializ(List<Integer> rhoms, List<Integer> priority) {
        List<Integer> re = new ArrayList<Integer>();
        for (int id : priority) {
            if (rhoms.contains(id))
                re.add(id);
        }
        Collections.reverse(re);
        return re;
    }

    private <T> List<Integer> getKeysinMap(Map<Integer, T> resultK) {
        List<Integer> re = new ArrayList<Integer>();
        List<Entry<Integer, T>> entryList = new ArrayList<Map.Entry<Integer, T>>(resultK.entrySet());
        Iterator<Entry<Integer, T>> iter = entryList.iterator();
        Entry<Integer, T> tmpEntry = null;
        while (iter.hasNext()) {
            tmpEntry = iter.next();
            re.add(tmpEntry.getKey());
        }
        return re;
    }

}
