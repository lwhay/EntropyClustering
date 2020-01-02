/**
 * 
 */
package ics.whu.edu.cn.madrix.clustering.density;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import ics.whu.edu.cn.madrix.basis.MadrixUtils;

/**
 * @author Administrator
 *
 */
public abstract class AbstractDensityClustering {
    protected boolean isDTW = false;

    protected double[][] data;

    protected double cutoff;

    protected double[] rankingdist;

    protected int[] ouputLabels;

    protected double radius;

    protected int dim;

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

    public static double distance(double[] s, double[] t, boolean isDTW) {
        double d = .0f;
        if (!isDTW) {
            for (int idx = 0; idx < s.length; idx++) {
                d += (s[idx] - t[idx]) * (s[idx] - t[idx]);
            }
            d = Math.sqrt(d);
        } else {
            d = MadrixUtils.vectorDTWDistance(s, t);
        }
        return d;
    }

    public double getMaximalDist() {
        return rankingdist[rankingdist.length - 1];
    }

    public double getCutoff() {
        return cutoff;
    }

    public double getRadius() {
        return radius;
    }
}
