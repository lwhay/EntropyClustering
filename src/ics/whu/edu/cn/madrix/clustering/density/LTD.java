package ics.whu.edu.cn.madrix.clustering.density;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

public class LTD implements IClustering {

    //private double[][] dists;
    private int KT;
    private final CTD ctd;
    private final STD std;
    private final ATD atd;
    private final Map<Integer, Double> resultK;
    private final int Ker;
    private int KKm;
    //private int DH;
    private Map<Integer, Integer> Index = new HashMap<>();
    private Map<Integer, Integer> Index0 = new HashMap<>();
    private Map<Integer, Integer> Ind = new HashMap<>();
    private Map<Integer, List<Integer>> distSortedNbs = new HashMap<>();
    private Map<Integer, Double> initEntropies = new HashMap<>();
    private Map<Integer, Integer> GKer = new LinkedHashMap<Integer, Integer>();
    private List<Integer> rhoms = new ArrayList<Integer>();
    private Map<Integer, List<Integer>> rePath = new HashMap<Integer, List<Integer>>();
    private List<Integer> centers = new ArrayList<Integer>();
    private Map<Integer, Integer> labels = new LinkedHashMap<Integer, Integer>();

    public LTD(CTD ctd, STD std, ATD atd, Map<Integer, Double> resultK, int Ker) {
        this.ctd = ctd;
        this.std = std;
        this.atd = atd;
        this.resultK = resultK;
        this.Ker = Ker;
        rePath = std.getrePath();
        //dists = ctd.getdists();
        centers = getKeysinMap(resultK);
        Ind = ctd.getIndex();
        rhoms = std.getrhoms();
        Index0 = ctd.getIndex0();
        Index = std.getIndex();
        KT = Integer.valueOf(ctd.getKT());
        //KKm = ctd.getKKm();
        KKm = ctd.getKKm();
        distSortedNbs = ctd.getNbs();
        initEntropies = atd.getEntropies();
    }

    @SuppressWarnings("unused")
    private boolean IsexpansibleT(int ix, int s1) {
        List<Integer> t1 = atd.getKeyByValue(Ind, s1);
        List<Integer> ts1 = t1;
        ts1.retainAll(rhoms); // intersection
        int end = std.ComputeClassdist(ix, ts1);

        List<Integer> t2 = atd.getKeyByValue(Ind, ix);
        List<Integer> ts2 = t2;
        ts2.retainAll(rhoms); // intersection
        int start = std.ComputeClassdist(end, ts2);

        return Isexpansible(start, end); //start , end       
    }

    private boolean Isexpansible(Integer i, int is) {
        boolean re = false;
        List<Integer> t1 = atd.getKeyByValue(Index0, i);
        List<Integer> t2 = atd.getKeyByValue(Index0, is);
        if (expansA(i, t1, t2)) {
            re = true;
        }
        return re;
    }

    private int relyIndex(Integer i, List<Integer> gKernels) {
        int is = rePath.get(i).get(rePath.get(i).size() - 1);
        int is0 = 0;
        while (gKernels.contains(is) == false) {
            if (is0 != rePath.get(is).get(rePath.get(is).size() - 1)) {
                is0 = is;
                is = rePath.get(is).get(rePath.get(is).size() - 1);
            } else {
                is = is0;
                break;
            }
        }
        return is;
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

    public boolean expansA(int begin, List<Integer> t1, List<Integer> t2) {
        Set<Integer> ba = maxavgNNA(t1);
        Set<Integer> bs = maxavgNNA(t2);
        boolean re = true;
        Set<Integer> res = new HashSet<>();
        res.clear();
        res.addAll(ba);
        res.retainAll(bs);
        if (res.size() <= KT)//(int) KT / 2)
            re = false;
        return re;
    }

    private Set<Integer> maxavgNNA(List<Integer> nns) {
        Set<Integer> bs = new HashSet<Integer>();
        for (Integer id : nns) {
            Set<Integer> ba = getKTNN(id);
            for (int ob : ba) {
                if (ctd.expansible(id, ob))
                    bs.add(ob);
            }
        }
        return bs;
    }

    private Set<Integer> getKTNN(int center) {
        Set<Integer> nns = new HashSet<>();
        for (int i = 0; i < KKm; i++) { //(1 + DH )* 
            if (ctd.expansible(center, distSortedNbs.get(center).get(i)))
                nns.add(distSortedNbs.get(center).get(i));
        }
        return nns;
    }

    @Override
    public void action() {
        Iterator<Integer> iter = resultK.keySet().iterator();
        int count = 1;
        List<Integer> GKernels = new ArrayList<Integer>();
        while (iter.hasNext()) {
            int xx = iter.next();
            GKernels.add(xx);
            //System.out.print(xx + 1 + "\t");
            if (count == Ker)
                break;
            count = count + 1;
        }

        System.out.print("initialize Kernel centers :");
        for (int id : GKernels)
            System.out.print(id + 1 + "\t");
        System.out.print("\n");

        List<Integer> priority = ctd.getPriority();
        Collections.reverse(priority);

        for (int id : priority) {
            if (GKernels.contains(id))
                GKer.put(id, GKer.size());
            if (GKer.size() == GKernels.size())
                break;
        } // ranking kernels  

        for (Integer id : GKernels) {
            centers.remove(id);
        }

        Map<Integer, Double> cerl = new LinkedHashMap<Integer, Double>();
        for (int id : priority) {
            if (centers.contains(id))
                cerl.put(id, initEntropies.get(id));
            if (cerl.size() == centers.size())
                break;
        }

        if (centers.size() > 0) {
            cerl = ctd.sortMapByValue(cerl);
            Iterator<Integer> iters = cerl.keySet().iterator(); // GKer  
            while (iters.hasNext()) {
                Integer i = iters.next();
                Integer is = relyIndex(i, GKernels);
                if (Isexpansible(i, is) && (GKernels.contains(is) || centers.contains(is))) {
                    List<Integer> tc = atd.getKeyByValue(Index, i);
                    for (int id : tc) {
                        Index.put(id, is);
                    }
                    centers.remove(i);
                }
            }

            Map<Integer, Double> cerls = new LinkedHashMap<Integer, Double>();
            for (int id : centers)
                cerls.put(id, cerl.get(id));

            if (cerls.size() > 0) {
                cerls = ctd.sortMapByValue(cerls);
                iters = cerls.keySet().iterator();
                while (iters.hasNext()) {
                    HashSet<Integer> ft = new  HashSet<Integer>(); 
                    Integer i = iters.next();
                    ft.addAll(GKernels);
                    ft.addAll(centers);
                    ft.remove(i);
                    List<Integer> tmp = new ArrayList<>(ft);
                    while (tmp.size() > 0) {
                        Integer str = ctd.ComputeMinrhorelay(i, tmp, KT + 1); //ctd.ComputediffrhoM(i, Index0, tmp);
                        if (Isexpansible(i, str)) { //IsexpansibleT
                            List<Integer> tc = atd.getKeyByValue(Index, i);
                            for (int id : tc) {
                                Index.put(id, str);
                            }
                            centers.remove(i);
                            break;
                        } else
                            tmp.remove(str);
                    }
                    if (tmp.size() == 0) {
                        Map<Integer, Double> tt = new HashMap<Integer, Double>();
                        for (Integer id : GKernels) {                            
                            int l1 = std.ComputeClassdist(id, GKernels);
                            tt.put(id, std.ComputeMindist(id, l1));
                        }
                        Integer i1 = ctd.getMinValue(tt);
                        GKernels.remove(i1);
                        GKer.put(i, GKer.get(i1));
                        GKer.remove(i1);
                        GKernels.add(i);
                        ft = new  HashSet<Integer>();
                        ft.addAll(GKernels);
                        ft.addAll(centers);
                        ft.remove(i1);
                        tmp = new ArrayList<>(ft);
                        while (tmp.size() > 0) {
                            Integer str = ctd.ComputeMinrhorelay(i1, tmp, KT + 1); //ctd.ComputediffrhoM(i1, Index0, tmp);
                            if (Isexpansible(i1, str)) {
                                List<Integer> tc = atd.getKeyByValue(Index, i1);
                                for (int id : tc) {
                                    Index.put(id, Index.get(str));
                                }
                                break;
                            } else
                                tmp.remove(str);
                        }
                        if (tmp.size() == 0) {
                            ft = new  HashSet<Integer>();
                            ft.addAll(GKernels);
                            ft.addAll(centers);
                            ft.remove(i1);
                            tmp = new ArrayList<>(ft);
                            List<Integer> tc = atd.getKeyByValue(Index, i1);
                            Integer newindex = ctd.ComputeMinrhorelay(i1, tmp, KKm + 1);
                            for (int ii = 0; ii < tc.size(); ii++) {
                                Index.put(tc.get(ii), newindex);
                            }
                        }
                        System.out.print("update the Kernel centers :");
                        for (int id : GKernels)
                            System.out.print(id + 1 + "\t");
                        System.out.print("\n");
                    }
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

            for (Integer id : GKernels) {
                List<Integer> tc = atd.getKeyByValue(Index, id);
                for (int nd : tc) {
                    labels.put(nd, GKer.get(id));
                }
            }
        } else
            System.out.print("the number of local cores is less than Ker!!");
    }

    @Override
    public int[] export() {
        int[] results = new int[labels.size()];
        for (int i = 0; i < labels.size(); i++) {
            results[i] = labels.get(i);
        }
        return results;
    }
}
