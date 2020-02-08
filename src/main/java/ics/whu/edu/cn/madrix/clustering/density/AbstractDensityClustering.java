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
import ics.whu.edu.cn.madrix.common.exceptions.MadrixException;

/**
 * @author Administrator
 *
 */
public abstract class AbstractDensityClustering {
    protected static final int DEFAULT_NORM_POW = 2;

    protected int type = -1;

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

    public static double distance(double[] s, double[] t, int type) throws MadrixException {
        double d = .0f;
        switch (type) {
            case 0: {
                for (int idx = 0; idx < s.length; idx++) {
                    d += (s[idx] - t[idx]) * (s[idx] - t[idx]);
                }
                d = Math.sqrt(d);
                break;
            }
            case 1: {
                d = MadrixUtils.vectorDTWDistance(s, t);
                break;
            }
            case 2: {
                d = MadrixUtils.vectorPDistance(s, t, DEFAULT_NORM_POW);
            }
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
